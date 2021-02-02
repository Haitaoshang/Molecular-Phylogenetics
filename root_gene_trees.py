#!/usr/bin/env python2.7
#coding: utf-8

import ete3
import os
import re
import operator
import multiprocessing
import pickle as pkl
import pandas as pd
import numpy as np
import random
from commands import getoutput
import itertools
import seaborn as sns
import collections
import operator
from matplotlib import pyplot as plt
from time import time
import random
import plotly
import plotly.plotly as ptl
from plotly import graph_objs as go
import plotly.figure_factory as ff
from sklearn.neighbors import KernelDensity
ptl.sign_in('lthiberiol', 'm15ikp59lt')

class cd:
    """
    Context manager for changing the current working directory
    """
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

def jaccard(array1, array2):
    return 1 - len(array1.intersection(array2))/float(len(array1.union(array2)))

def root_like( ref_tree, tree2 ):
    tree_to_root = tree2.copy()
    for node in sorted( ref_tree.children, key=len ):
        if node.is_leaf():
            leaf = tree_to_root.get_leaves_by_name(node.name)[0]
            tree_to_root.set_outgroup(leaf)
            break
        else:
            is_it_monophyletic, clade_type, fucking_up = tree_to_root.check_monophyly(node.get_leaf_names(), 'name', unrooted=False)
            if is_it_monophyletic:
                equivalent = tree_to_root.get_common_ancestor(node.get_leaf_names())
                tree_to_root.set_outgroup(equivalent)
            else:
                gene_tree.set_outgroup(fucking_up.pop())
                equivalent = tree_to_root.get_common_ancestor(node.get_leaf_names())
                tree_to_root.set_outgroup(equivalent)
            break

    return tree_to_root

ncbi = ete3.NCBITaxa()

os.chdir('/work/Alphas_and_Cyanos')
species_tree = ete3.Tree('rooted_partitions-with_BB_support.treefile', format=0)
gene_tree    = ete3.Tree('ranger_input_trees/000070.tree',             format=0)

genomes = set()
for leaf in gene_tree.get_leaves():
    genome = re.match('GC[AF]\d+', leaf.name).group()
    leaf.add_feature('genome', genome)
    genomes.add(genome)

pruned_species_tree = species_tree.copy()
pruned_species_tree.prune(genomes)

child1           = set(pruned_species_tree.children[0].get_leaf_names())
child2           = set(pruned_species_tree.children[1].get_leaf_names())
gene_tree_leaves = set(gene_tree.get_leaves())

loop_count        = 0
possible_rootings = {}
for node in gene_tree.traverse():
    if node.is_leaf():
        continue

    bipartition1 = set([leaf.genome for leaf in node.get_leaves()])
    bipartition2 = genomes.difference(bipartition1)

    if not bipartition2:
        continue

    loop_count += 1
    node.add_feature('name', 'node_%i' %loop_count)

    jaccard_distance1 = jaccard(child1, bipartition1) + jaccard(child2, bipartition2)
    jaccard_distance2 = jaccard(child1, bipartition2) + jaccard(child2, bipartition1)
    if jaccard_distance1 < jaccard_distance2:
        possible_rootings[node.name] = jaccard_distance1
    else:
        possible_rootings[node.name] = jaccard_distance2

sorted_by_jaccard = sorted(possible_rootings.items(), key=operator.itemgetter(1))
optimal_root      = gene_tree.search_nodes(name=sorted_by_jaccard[0][0])[0]
gene_tree.set_outgroup(optimal_root)
gene_tree.write(outfile='test_optimal_rootings/000070.tax_optimal_root.tree', format=0)

########################################################################################################################
# Compare rooting with RANGER's                                                                                        #
########################################################################################################################

gene_tree_optRooted = root_like(ete3.Tree('test_optimal_rootings/000070.ranger_rooting.tree'), gene_tree)

random_node_sample        = random.sample(list(gene_tree_optRooted.traverse()), 100)
optimal_root_in_optRooted = gene_tree_optRooted.search_nodes(name=sorted_by_jaccard[0][0])[0]
ancestors                 = optimal_root_in_optRooted.get_ancestors()
ancestors.insert(0, optimal_root_in_optRooted)

with cd('test_optimal_rootings/random_rootings'):
    count = 0
    for node in random_node_sample:
        gene_tree_optRooted.set_outgroup(node)
        out    = open('random_gene_tree_rootings-%i.ranger_input' %count, 'wb')
        out.write(species_tree.write(format=9))
        out.write('\n')
        out.write(gene_tree_optRooted.write(format=9))
        out.close()
        count += 1

with cd('test_optimal_rootings'):
    count = 0
    for ancestor in ancestors:
        gene_tree_optRooted.set_outgroup(ancestor)
        print ancestor.name
        print gene_tree_optRooted.get_topology_id()
        out    = open('steps_to_optRoot-%i.ranger_input' %count, 'wb')
        out.write(species_tree.write(format=9))
        out.write('\n')
        out.write(gene_tree_optRooted.write(format=9))
        out.close()
        count += 1

reconciliation_costs = [int(cost) for cost in re.findall(
    'The minimum reconciliation cost is: (\d+)',
    getoutput('grep "The minimum reconciliation cost is" test_optimal_rootings/steps_to_optRoot-{0..13}.ranger_input.out'
))]

null_reconciliation_costs = [int(cost) for cost in re.findall(
    'The minimum reconciliation cost is: (\d+)',
    getoutput('grep "The minimum reconciliation cost is" test_optimal_rootings/random_rootings/random_gene_tree_rootings-*.ranger_input.out'
))]

mad_rooting_cost = int(re.match('The minimum reconciliation cost is: (\d+)', open('test_optimal_rootings/mad_rooting.ranger_input.out').readlines()[-3]).group(1))

g = sns.distplot(null_reconciliation_costs, kde_kws=dict(shade=True, color='black'), hist_kws=dict(color='#24476B', alpha=0.7))
sns.distplot(null_reconciliation_costs, hist=False, kde_kws=dict(color='white'), ax=g)
g.axes.plot((reconciliation_costs[0], reconciliation_costs[0]), g.axes.get_ybound(), color='#BF4046', alpha=1, label='Root matching species tree')
g.axes.plot((reconciliation_costs[-1], reconciliation_costs[-1]), g.axes.get_ybound(), color='#246B40', alpha=1, label='Most parsimonious DTL cost')
g.axes.plot((mad_rooting_cost, mad_rooting_cost), g.axes.get_ybound(), color='#613399', alpha=1, label='Root as predicted by MAD')
g.set_aspect(100)
g.figure.set_size_inches(w=20, h=10)
g.legend()
g.figure.savefig('test_optimal_rootings/null_costs.pdf', bbox_inches='tight', dpi=300)
g.clear()
plt.close()

upper_std = np.mean(null_reconciliation_costs)+np.std(null_reconciliation_costs)
lower_std = np.mean(null_reconciliation_costs)-np.std(null_reconciliation_costs)
plot = sns.barplot(x=range(1,len(reconciliation_costs)+1), y=reconciliation_costs, palette=sns.color_palette('seismic', n_colors=14, desat=None))
plot.axes.plot(plot.axes.get_xbound(), (upper_std, upper_std), color='#BF4046', linestyle='--', alpha=0.6, label='Standard deviantion of random root sample')
plot.axes.plot(plot.axes.get_xbound(), (lower_std, lower_std), color='#BF4046', linestyle='--', alpha=0.6)
plot.axes.plot(plot.axes.get_xbound(), (mad_rooting_cost, mad_rooting_cost), color='#613399', linestyle='--', alpha=0.6, label='Cos at root as predicted by MAD')
plot.set_ybound(lower=1550, upper=1570)
plot.axes.set_xlabel('From root matching species tree to most parsimonious')
plot.axes.set_ylabel('Reconciliation costs')
plot.figure.set_size_inches(w=15, h=10)
plot.legend()
plot.figure.savefig('test_optimal_rootings/reconciliation_costs.pdf', bbox_inches='tight', dpi=300)
plot.clear()
plt.close()

mad_costs = {}
for group in single_optimal_rooting:
    dtl = int(re.match('The minimum reconciliation cost is: (\d+)', open('comparing_rooting_methods/mad_reconciliations/%s/%s-MAD.ranger_out1' %(group,group)).readlines()[-3]).group(1))
    mad_costs[group] = dtl

ranger_optimal_costs = {}
for group in single_optimal_rooting:
    dtl = int(re.match('The minimum reconciliation cost is: (\d+)', open('reconciliations/{group}/{group}.reconciliation1'.format(group=group)).readlines()[-3]).group(1))
    ranger_optimal_costs[group] = dtl

delta_dtl = {group:mad_costs[group]-ranger_optimal_costs[group] for group in mad_costs.keys()}

fig, axs = plt.subplots(nrows=2)
sns.distplot(delta_dtl.values(), ax=axs[0])
axs[0].set_title("delta DTL (MAD's value is always higher)")
sns.distplot([float(delta_dtl[group])/ranger_optimal_costs[group] for group in delta_dtl.keys()], ax=axs[0])
axs[1].set_title("delta DTL as percentages of the optimal cost (MAD's value is always higher)")
fig.tight_layout()
fig.savefig('MAD_ranger_delta.pdf', dpi=600)

g = sns.distplot(delta_dtl.values())
g.axes.set_title("delta DTL of the optimal cost (MAD's value is always higher)")
g.figure.set_size_inches(w=20, h=10)
g.figure.savefig('MAD_ranger_delta.pdf', bbox_inches='tight', dpi=600)
g.clear()
plt.close()

g = sns.distplot([float(delta_dtl[group])/ranger_optimal_costs[group] for group in delta_dtl.keys()])
g.axes.set_title("delta DTL as percentages of the optimal cost (MAD's value is always higher)")
g.figure.set_size_inches(w=20, h=10)
g.figure.savefig('MAD_ranger_delta_fraction.pdf', bbox_inches='tight', dpi=600)
g.clear()
plt.close()

single_optimal_rooting = []
with cd('reconciliations/ranger_roots'):
    for group in os.listdir('.'):
        if not os.path.isdir(group) or not os.path.isfile('%s/%s.reconciliation1' %(group, group)):
            continue
        if os.path.isfile('%s/%s.optResolution1.ranger_input' %(group, group)) and not os.path.isfile('%s/%s.optResolution2.ranger_input' %(group, group)):
            single_optimal_rooting.append(group)

out = open('groups_with_single_optimal_DTL_root.list', 'wb')
out.write('\n'.join(single_optimal_rooting))
out.close()

transfers = [hgt.assess_transfers(transfer_file) for transfer_file in transfer_files]
results = {}
[results.update(transfer) for transfer in transfers]
out = open('aggregated/mad_index_transfers.pkl', 'w')
pkl.dump(results, out)
out.close()
transfers = [hgt.visualize_tree(*transfer, output_folder='/work/Alphas_and_Cyanos/index_transfer_trees/mad') for transfer in results.items()]


########################################################################################################################
def root_distances(group):
    ref_tree = ete3.Tree(open('/work/Alphas_and_Cyanos/comparing_rooting_methods/ranger/{group}/{group}.reconciliation1'.format(group=group)).readlines()[7], format=1)
    tree2    = ete3.Tree(open('/work/Alphas_and_Cyanos/comparing_rooting_methods/mad_reconciliations/{group}/{group}-MAD.ranger_out1'.format(group=group)).readlines()[7], format=1)
    for node in sorted( ref_tree.children, key=len ):
        if node.is_leaf():
            leaf = tree2.get_leaves_by_name(node.name)[0]
            root_distances = tree2.get_distance(leaf, topology_only=True)
            return (root_distances, root_distances/(2*len(tree2)-3))
        else:
            is_it_monophyletic, clade_type, fucking_up = tree2.check_monophyly(node.get_leaf_names(), 'name', unrooted=False)
            if is_it_monophyletic:
                equivalent = tree2.get_common_ancestor(node.get_leaf_names())
                root_distances = tree2.get_distance(equivalent, topology_only=True)
                return (root_distances, root_distances/(2*len(tree2)-3))
            else:
                continue
    return None

pool = multiprocessing.Pool(processes=4)
distances = pool.map(root_distances, single_optimal_rooting)
pool.close()
pool.join()
print 'yeah'

raw  = {}
norm = {}
for group, (i,j) in zip(single_optimal_rooting, distances):
    raw[group]  = i
    norm[group] = j

fig, axs = plt.subplots(nrows=2)
sns.distplot(raw.values(), ax=axs[0])
axs[0].set_title('Bipartitions between RANGER and MAD rootings')
sns.distplot(norm.values(), ax=axs[1])
axs[1].set_title('Normalized bipartitions between RANGER and MAD rootings')
fig.tight_layout()
fig.set_size_inches(15,10)
fig.savefig('distances_between_ranger_and_mad_root_positions1.pdf', dpi=600)
fig.clear()
plt.close()

fig, ax = plt.subplots()
ax.plot([raw[group]  for group in single_optimal_rooting],
        [norm[group] for group in single_optimal_rooting],
        'ko', alpha=0.6)
ax.set_xlabel('Bipartitions between RANGER and MAD rootings')
ax.set_ylabel('Normalized bipartitions between RANGER and MAD rootings')
fig.tight_layout()
fig.set_size_inches(10,10)
fig.savefig('yeah.pdf', dpi=600)
fig.clear()
plt.close()

sorted_by_root_distance = sorted(norm.items(), key=operator.itemgetter(1))
########################################################################################################################
########################################################################################################################
def parse_supported_transfers(handle, threshold=0.6):
    text                      = handle.read()
    number_of_reconciliations = int(re.match('Processed (\d+) files', text).group(1))
    if number_of_reconciliations != 20:
        return None
    transfers = re.findall('^.*, Transfers = [^0]\d+?\], \[Most Frequent mapping --> (n\S+), (\d+) times\], \[Most Frequent recipient --> (n\S+), (\d+) times\].', text, re.M)
    supported_pairs = []
    for donor, donor_support, recipient, recipient_support in transfers:
        if int(donor_support) < threshold*number_of_reconciliations or int(recipient_support) < threshold*number_of_reconciliations:
            continue
        supported_pairs.append((donor, recipient))

    return supported_pairs
########################################################################################################################
########################################################################################################################

os.chdir('/work/Alphas_and_Cyanos/comparing_rooting_methods')

def parse_supported_transfers(handle, threshold=0.9):
    text                      = handle.read()
    number_of_reconciliations = int(re.match('Processed (\d+) files', text).group(1))
    if number_of_reconciliations != 20:
        return None
    transfers = re.findall('^.*, Transfers = [^0]\d?\], \[Most Frequent mapping --> (\S+), (\d+) times\], \[Most Frequent recipient --> (\S+), (\d+) times\].', text, re.M)
    supported_pairs = []
    for donor, donor_support, recipient, recipient_support in transfers:
        if int(donor_support) < threshold*number_of_reconciliations or int(recipient_support) < threshold*number_of_reconciliations:
            continue
        supported_pairs.append((donor, recipient, int(donor_support), int(recipient_support)))

    return supported_pairs

with cd('/work/Alphas_and_Cyanos/comparing_rooting_methods/ranger/'):
    ranger_supported_pairs = {}
    for group in single_optimal_rooting:
        if not os.path.isdir(group) or not os.path.isfile('{group}/aggregated'.format(group=group)):
            os.mkdir(group)
            os.system('cp /work/Alphas_and_Cyanos/reconciliations/{group}/{group}.reconciliation* {group}/'.format(group=group))
            os.system('/work/ranger/CorePrograms/AggregateRanger_recipient {group}/{group}.reconciliation > {group}/aggregated'.format(group=group))

        handle = open('%s/aggregated' %group)
        ranger_supported_pairs[group] = parse_supported_transfers(handle, threshold=0.8)

with cd('/work/Alphas_and_Cyanos/comparing_rooting_methods/mad_reconciliations/'):
    mad_supported_pairs = {}
    for group in single_optimal_rooting:
        if not os.path.isdir(group) or not os.path.isfile('{group}/{group}-MAD.ranger_out1'.format(group=group)):
            continue
        if not os.path.isfile('{group}/aggregated'.format(group=group)):
            os.system('/work/ranger/CorePrograms/AggregateRanger_recipient {group}/{group}-MAD.ranger_out > {group}/aggregated'.format(group=group))

        handle = open('%s/aggregated' %group)
        mad_supported_pairs[group] = parse_supported_transfers(handle, threshold=0.8)

with cd('/work/Alphas_and_Cyanos/comparing_rooting_methods/random_root_reconciliations/'):
    random_root_supported_pairs = {}
    for group in single_optimal_rooting:
        if not os.path.isdir(group) or not os.path.isfile('{group}/{group}-random_root.ranger_out1'.format(group=group)):
            continue
        if not os.path.isfile('{group}/aggregated'.format(group=group)):
            os.system('/work/ranger/CorePrograms/AggregateRanger_recipient {group}/{group}-random_root.ranger_out > {group}/aggregated'.format(group=group))

        handle = open('%s/aggregated' %group)
        random_root_supported_pairs[group] = parse_supported_transfers(handle, threshold=0.8)

tmp    = []
[tmp.extend(value) for value in ranger_supported_pairs.values()]
ranger = collections.Counter(tmp)
tmp    = []
[tmp.extend(value) for value in mad_supported_pairs.values()]
mad = collections.Counter(tmp)
tmp    = []
[tmp.extend(value) for value in random_root_supported_pairs.values() if value]
random_roots= collections.Counter(tmp)

a = np.asarray(ranger.values())
b = np.asarray(mad.values())
c = np.asarray(random_roots.values())
for n in set().union(a,b):

    tmp_a = a[a>=n]
    tmp_b = b[b>=n]
    tmp_c = c[c>=n]

    fig, ax = plt.subplots()
#    sns.kdeplot(tmp_c, color='#246B40', bw=3, label='Random rooting', ax=ax)
    sns.kdeplot(tmp_a, color='#24476B', bw=2, label='Ranger', ax=ax)
    sns.kdeplot(tmp_b, color='#BF4046', bw=2, label='MAD', ax=ax)
    ax.set_xlim(left=0)
    fig.tight_layout()
    fig.set_size_inches(15,8)
    fig.savefig('yeah/yeah%i.png' %n)
    fig.clear()
    plt.close()

    tmp_a = a[a>=n]
    total = float(sum(tmp_a))
    tmp = list(set(a[a>=n]))
    tmp.sort()
    ax.plot(tmp, [tmp_a[tmp_a<=i].shape[0] for i in tmp], 'b-')

    tmp_b = b[b>=n]
    total = float(sum(tmp_b))
    tmp = list(set(b[b>=n]))
    tmp.sort()
    ax.plot(tmp, [tmp_b[tmp_b<=i].shape[0] for i in tmp], 'r-')

    ax.legend()
    fig.tight_layout()
    fig.savefig('yeah/yeah%i.png' %n)
    fig.clear()
    plt.close()

a = np.asarray(ranger.values())
b = np.asarray(mad.values())
c = np.asarray(random_roots.values())
fig, ax = plt.subplots()
#sns.kdeplot(random_roots.values(), color='#246B40', bw=1, label='Random rooting', ax=ax)
#sns.kdeplot(ranger.values(),       color='#24476B', bw=1, label='Ranger', ax=ax)
#sns.kdeplot(mad.values(),          color='#BF4046', bw=1, label='MAD', ax=ax)
sns.kdeplot(c[c>=np.percentile(c, 99)], cumulative=True, color='#246B40', bw=1, label='Random rooting', ax=ax)
sns.kdeplot(a[a>=np.percentile(c, 99)], cumulative=True, color='#24476B', bw=1, label='Ranger', ax=ax)
sns.kdeplot(b[b>=np.percentile(c, 99)], cumulative=True, color='#BF4046', bw=1, label='MAD', ax=ax)
fig.tight_layout()
fig.set_size_inches(15,8)
fig.savefig('donor_recipient_pair_frquency95.pdf', dpi=600)
fig.clear()
plt.close()

fig = ff.create_distplot([a[a>=np.percentile(a, 95)],b[b>=np.percentile(b, 95)]], ['ranger', 'mad'], bin_size=1, show_hist=False)
fig['layout'].update(title='Donor/recipient frequency distribution')
fig['layout'].update(xaxis=dict(title='# occurrences in gene families'))
fig['layout'].update(yaxis=dict(title='# of gene families'))
fig['layout'].update(hovermode='closest')
plotly.offline.plot(fig, filename='yeah.html')

fig, ax = plt.subplots()
sns.boxplot(data=[a[a>=np.percentile(a, 99)], b[b>=np.percentile(b, 99)], c[c>=np.percentile(c, 99)]], palette=['#24476B', '#BF4046', '#246B40'], ax=ax)
fig.tight_layout()
fig.set_size_inches(15,8)
fig.savefig('30/donor_recipient_pair_frquency95-boxplot.pdf', dpi=600)
fig.clear()
plt.close()

delta_ranger = [ranger[group]-random_roots[group] for group in random_roots.keys()]
delta_mad = [mad[group]-random_roots[group] for group in random_roots.keys()]
fig, ax = plt.subplots()
sns.kdeplot(delta_ranger, color='#24476B', bw=1, ax=ax, label='Ranger donor/recipient frequency - random rooting')
sns.kdeplot(delta_mad,    color='#BF4046', bw=1, ax=ax, label='MAD donor/recipient frequency - random rooting'   )
ax.set_xlabel('# of occurrence among gene families')
ax.set_ylabel('Ranger frequency - MAD frequency')
fig.tight_layout()
fig.set_size_inches(15,8)
fig.savefig('30/donor_recipient_pair_frquency-random_delta.pdf', dpi=600)
fig.clear()
plt.close()

pair_frequencies = list(set().union(*[a,b]))
pair_frequencies.sort()
delta = [a[a==frequency].shape[0]/float(a.shape[0]) - b[b==frequency].shape[0]/float(b.shape[0]) for frequency in pair_frequencies]
pair_frequencies.insert(0,0)
delta.insert(0,0)
delta=np.asarray(delta)
fig, ax = plt.subplots()
ax.plot(pair_frequencies, delta, 'k-', label='# of ranger pairs - # of MAD pairs')
ax.plot([],[],linewidth=5, label='Enriched in MAD', color='r',alpha=0.5)
ax.plot([],[],linewidth=5, label='Enriched in RANGER', color='b',alpha=0.5)
ax.fill_between(pair_frequencies, 0, delta, where=(delta <= 0), facecolor='r', alpha=0.5, interpolate=True)
ax.fill_between(pair_frequencies, 0, delta, where=(delta > 0), facecolor='b', alpha=0.5, interpolate=True)
ax.set_xlabel('# of occurrence among gene families')
ax.set_ylabel('Ranger frequency - MAD frequency')
fig.set_size_inches(15,8)
fig.tight_layout()
ax.legend()
fig.savefig('30/ranger_mad_delta.pdf', dpi=300,)
fig.clear()
plt.close()

for n in set().union(a,b):
    fig, ax = plt.subplots()

    tmp_a = a[a>=n]
    tmp_b = b[b>=n]
    sns.kdeplot(tmp_a, cumulative=True, color='#24476B', bw=2, label='Ranger', ax=ax)
    sns.kdeplot(tmp_b, cumulative=True, color='#BF4046', bw=2, label='MAD', ax=ax)
    ax.set_xlim(left=0)
    ax.legend(title='starting value is: %i' %n, loc=4)
    fig.tight_layout()
    fig.set_size_inches(15,8)
    fig.savefig('yeah_cdf/yeah%i.png' %n)
    fig.clear()
    plt.close()

    fig, ax = plt.subplots()
    tmp_a = a[a>=n]
    total = float(sum(tmp_a))
    tmp = list(set(a[a>=n]))
    tmp.sort()
    ax.plot(tmp, [tmp_a[tmp_a<=i].shape[0] for i in tmp], ls='-', color='#24476B', label='Ranger')

    tmp_b = b[b>=n]
    total = float(sum(tmp_b))
    tmp = list(set(b[b>=n]))
    tmp.sort()
    ax.plot(tmp, [tmp_b[tmp_b<=i].shape[0] for i in tmp], ls='-', color='#BF4046', label='MAD')

    ax.set_xlim(left=0)
    ax.legend(title='starting value is: %i' %n, loc=4)
    fig.set_size_inches(15,8)
    fig.tight_layout()
    fig.savefig('yeah/yeah%i.png' %n)
    fig.clear()
    plt.close()

fig, ax = plt.subplots(figsize=(30,10))
tmp_a = a[a>=8]
tmp = list(set(tmp_a))
tmp.sort()
ax.plot(tmp, [tmp_a[tmp_a<=i].shape[0] for i in tmp], 'b-')

tmp_b = b[b>=8]
tmp = list(set(tmp_b))
tmp.sort()
ax.plot(tmp, [tmp_b[tmp_b<=i].shape[0] for i in tmp], 'r-')

fig.tight_layout()
fig.savefig('../yeah.png')
fig.clear()
plt.close()

a_percentiles = {}
b_percentiles = {}
for n in range(90,101, 1):
    tmp_a = a[(a<np.percentile(a, n)) & (a>=np.percentile(a, n-1))]
    a_percentiles[n] = (tmp_a.mean(), tmp_a.std())
    tmp_b = b[(b<np.percentile(b, n)) & (b>=np.percentile(b, n-1))]
    b_percentiles[n] = (tmp_b.mean(), tmp_b.std())


cumulative_delta = [sum(delta[:i]) for i in range(len(delta))]
fig, ax = plt.subplots(figsize=(15,8))
sns.regplot(x=np.asarray(pair_frequencies[2:]), y=np.asarray(cumulative_delta[2:]), ax=ax, truncate=True, logistic=True,
            scatter_kws=dict(color='#24476B'),
            line_kws=dict(color='#BF4046'))
ax.set_xlabel('# of occurrence among gene families')
ax.set_ylabel('CUMULATIVE Ranger frequency - MAD frequency')
fig.tight_layout()
fig.savefig('30/ranger_mad_cumulative_delta.pdf', dpi=300,)
fig.clear()
plt.close()

max_value = np.max(mad.values())
delta_families = []
mad_delta_families = []
ranger_delta_families = []
for frequency in range(max_value+1):
    delta_families.append(ranger.values().count(frequency)/float(len(ranger)) - mad.values().count(frequency)/float(len(mad)))
    mad_delta_families.append(mad.values().count(frequency)/float(len(mad)) - random_roots.values().count(frequency)/float(len(random_roots)))
    ranger_delta_families.append(ranger.values().count(frequency)/float(len(ranger)) - random_roots.values().count(frequency)/float(len(random_roots)))

cumulative_delta        = [sum(   delta_families[    :i]) for i in range(len(       delta_families))]
cumulative_delta_mad    = [sum(   mad_delta_families[:i]) for i in range(len(   mad_delta_families))]
cumulative_delta_ranger = [sum(ranger_delta_families[:i]) for i in range(len(ranger_delta_families))]
fig, ax = plt.subplots(figsize=(30,10))
#ax.plot(mad_delta_families, 'r-', ranger_delta_families, 'b-', delta_families, 'g-', lw=0.5)
ax.plot(cumulative_delta_mad, 'r-', cumulative_delta_ranger, 'b-', cumulative_delta, 'g-', lw=0.5)
fig.tight_layout()
fig.savefig('../yeah.pdf', dpi=300,)
fig.clear()
plt.close()
