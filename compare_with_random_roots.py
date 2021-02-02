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
                tree_to_root.set_outgroup(fucking_up.pop())
                equivalent = tree_to_root.get_common_ancestor(node.get_leaf_names())
                tree_to_root.set_outgroup(equivalent)
            break

    return tree_to_root

ncbi = ete3.NCBITaxa()
os.chdir('/work/Alphas_and_Cyanos')
species_tree = ete3.Tree('rooted_partitions-with_BB_support.treefile', format=0)
species_tree = ete3.Tree('rooted_partitions-with_named_branches.treefile', format=1)

########################################################################################################################
#                                                                                                                      #
#                                                                                                                      #
########################################################################################################################

single_optimal_rooting = []
with cd('reconciliations/ranger_roots'):
    for group in os.listdir('.'):
        if not os.path.isdir(group) or not os.path.isfile('%s/%s.reconciliation1' %(group, group)):
            continue
        if os.path.isfile('%s/%s.optResolution1.ranger_input' %(group, group)) and not os.path.isfile('%s/%s.optResolution2.ranger_input' %(group, group)):
            single_optimal_rooting.append(group)

def parse_supported_transfers(handle, threshold=0.9, leaves_allowed=False):
    text                          = handle.read()
    number_of_reconciliations     = int(re.match('Processed (\d+) files', text).group(1))
    if number_of_reconciliations != 20:
        return None
    if leaves_allowed:
        transfers = re.findall('^.*, Transfers = [^0]\d?\], \[Most Frequent mapping --> (\S+), (\d+) times\], \[Most Frequent recipient --> (\S+), (\d+) times\].', text, re.M)
    else:
        transfers = re.findall('^.*, Transfers = [^0]\d+?\], \[Most Frequent mapping --> (n\S+), (\d+) times\], \[Most Frequent recipient --> (n\S+), (\d+) times\].', text, re.M)
    supported_pairs = []
    for donor, donor_support, recipient, recipient_support in transfers:
        if int(donor_support) < threshold*number_of_reconciliations or int(recipient_support) < threshold*number_of_reconciliations:
            continue
        supported_pairs.append((donor, recipient))

    return supported_pairs

def rename_branches(reconciliation_file, tree):
    branches = re.findall('^(m\d+) = LCA\[(\S+), (\S+)\]:', reconciliation_file, re.M)
    for name, leaf1, leaf2 in branches:
        node = tree.get_common_ancestor(leaf1, leaf2)
        if node.name:
            continue
        node.name = name
    return tree

def parse_aggregated(folder, threshold=0.0, leaves_allowed=False):
    with cd(folder):
        aggregated    = open('aggregated').read()
        reconciliation_file = re.search('%s(-\S+?\.ranger_out|\.reconciliation)1' %folder, '\t'.join(os.listdir('.'))).group()
        gene_tree     = {'named':ete3.Tree(linecache.getline(reconciliation_file, 8), format=1)}

    gene_tree['support'] = root_like(gene_tree['named'], ete3.Tree('/work/Alphas_and_Cyanos/ranger_input_trees/%s.tree' %folder))
    gene_tree            = rename_branches(aggregated, gene_tree['support'])

    num_replicates = float(re.match('Processed (\d+) files', aggregated).group(1))

    if not leaves_allowed:
        transfers = re.findall('^(m\d+) = .*, Transfers = [^0]\d+?\], \[Most Frequent mapping --> (n\d+), (\d+) times\], \[Most Frequent recipient --> (n\d+), (\d+) times\].', aggregated, re.M)
    else:
        transfers = re.findall('^(m\d+) = .*, Transfers = [^0]\d+?\], \[Most Frequent mapping --> (\S+), (\d+) times\], \[Most Frequent recipient --> (\S+), (\d+) times\].',   aggregated, re.M)

    supported_transfers = []
    for donor_map, donor, ranger_confidence_donor, recipient, ranger_confidence_recipient in transfers:
        if int(ranger_confidence_donor) < threshold*num_replicates or int(ranger_confidence_recipient) < threshold*num_replicates:
            continue
        supported_transfers.append((donor_map, donor, recipient))

    selected_transfers = []
    for donor_map_name, donor_name, recipient_name in supported_transfers:
        donor_map = gene_tree.search_nodes(name=donor_map_name)[0]
        if donor_map.support < 95:
            continue

        selected_transfers.append((donor_name, recipient_name))

    return selected_transfers

########################################################################################################################
#                                                                                                                      #
# Highly supported transfers                                                                                           #
########################################################################################################################

os.chdir('/work/Alphas_and_Cyanos/comparing_rooting_methods')

to_remove = []
with cd('/work/Alphas_and_Cyanos/comparing_rooting_methods/random_root_reconciliations/'):
    for group in single_optimal_rooting:
        if not os.path.isdir(group) or not os.path.isfile('{group}/{group}-random_root.ranger_out1'.format(group=group)):
            to_remove.append(group)
for group in to_remove:
    single_optimal_rooting.remove(group)

with cd('/work/Alphas_and_Cyanos/comparing_rooting_methods/random_root_reconciliations/'):
    pool = multiprocessing.Pool(processes=6)
    random_root_supported_pairs = pool.map(parse_aggregated, single_optimal_rooting)

with cd('/work/Alphas_and_Cyanos/comparing_rooting_methods/ranger/'):
    pool = multiprocessing.Pool(processes=6)
    ranger_supported_pairs = pool.map(parse_aggregated, single_optimal_rooting)

with cd('/work/Alphas_and_Cyanos/comparing_rooting_methods/mad_reconciliations/'):
    pool = multiprocessing.Pool(processes=6)
    mad_supported_pairs = pool.map(parse_aggregated, single_optimal_rooting)

shared_transfers = []
for a,b,c in zip(mad_supported_pairs, ranger_supported_pairs, random_root_supported_pairs):
    shared_transfers.extend(set(c).intersection(a,b))

shared_by_mad_ranger = []
for a,b in zip(mad_supported_pairs, ranger_supported_pairs):
    shared_by_mad_ranger.extend(set(a).intersection(b))

shared_by_mad_random = []
for a,c in zip(mad_supported_pairs, random_root_supported_pairs):
    shared_by_mad_random.extend(set(a).intersection(c))

shared_by_random_ranger = []
for b,c in zip(ranger_supported_pairs, random_root_supported_pairs):
    shared_by_random_ranger.extend(set(b).intersection(c))

#
# characterize ranger missed transfers
shared    = []
exclusive = []
mad    = []
ranger = []
for key, transfers in ranger_supported_pairs.items():
    mad_exclusive    = set(mad_supported_pairs[key]).difference(ranger_supported_pairs[key])
    ranger_exclusive = set(ranger_supported_pairs[key]).difference(mad_supported_pairs[key])

    if not mad_exclusive or not ranger_exclusive:
        continue

    for donor, recipient in mad_exclusive:
        mad.append(species_tree.get_distance(donor, recipient))
    for donor, recipient in ranger_exclusive:
        ranger.append(species_tree.get_distance(donor, recipient))

    shared_transfers    = set(transfers).intersection(set(mad_supported_pairs[key]).union(random_root_supported_pairs[key]))
    exclusive_transfers = shared_transfers.symmetric_difference(transfers)

    if not shared_transfers or not exclusive_transfers:
        continue

    if len(shared_transfers) > len(exclusive_transfers):
        for donor, recipient in sample(shared_transfers, len(exclusive_transfers)):
            shared.append(species_tree.get_distance(donor, recipient))
        for donor, recipient in exclusive_transfers:
            exclusive.append(species_tree.get_distance(donor, recipient))
    elif len(shared_transfers) < len(exclusive_transfers):
        for donor, recipient in shared_transfers:
            shared.append(species_tree.get_distance(donor, recipient))
        for donor, recipient in sample(exclusive_transfers, len(shared_transfers)):
            exclusive.append(species_tree.get_distance(donor, recipient))
    else:
        for donor, recipient in shared_transfers:
            shared.append(species_tree.get_distance(donor, recipient))
        for donor, recipient in exclusive_transfers:
            exclusive.append(species_tree.get_distance(donor, recipient))

fig, ax = plt.subplots()
sns.kdeplot(ranger,    ax=ax, label='Shared with other roots',     color='green', shade=True)
sns.kdeplot(mad, ax=ax, label='Not shared with other roots', color='red'  , shade=True)
fig.tight_layout()
fig.savefig('yeah.pdf', dpi=300)













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






































