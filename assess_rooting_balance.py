import ete3
import os
import re
import subprocess
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
import pickle as pkl
from itertools import combinations
import multiprocessing
import pandas as pd

os.chdir('/work/Alphas_and_Cyanos/reconciliations/ranger_roots')

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

########################################################################################################################
#                                                                                                                      #
# ASSESS BALANCE                                                                                                       #
########################################################################################################################
def assess_tree_balance_ranger(folder):
    ratios = []
    trees = []
    count = 1
    num_transfers = []
    num_duplications = []
    sub_count = 1
    while os.path.isfile('%s/%s.optResolution%i.ranger_input' %(folder, folder, count)):
#        trees.append(ete3.Tree(open('%s/%s.optResolution%i.ranger_input' %(folder, folder, count)).readlines()[1]))
        tree     = ete3.Tree(open('%s/%s.optResolution%i.ranger_input' %(folder, folder, count)).readlines()[1])
        children = sorted(tree.children, key=len)
        ratio    = len(children[0])/float(len(children[1]))
        while sub_count <= count * 20:
            reconciliation = open('%s/%s.reconciliation%i' %(folder, folder, sub_count)).read()
            ratios.append(ratio)
            num_transfers.append(   len(re.findall('Transfer, Mapping --> \S+, Recipient --> \S+$', reconciliation, re.M)))
            num_duplications.append(len(re.findall('Duplication, Mapping --> \S+$',                 reconciliation, re.M)))
            sub_count += 1

        count += 1

#    for tree in trees:
#        children = sorted(tree.children, key=len)
#        ratios.append(len(children[0])/float(len(children[1])))

    num_internal_nodes = len(re.findall('^m\d+ = LCA', reconciliation, re.M))
    dtl_cost = re.search('The minimum reconciliation cost is: (\d+)', reconciliation).group(1)

    return {'balance_ratios':ratios, 'num_transfers':num_transfers, 'num_duplications':num_duplications,
            'num_internal_nodes':num_internal_nodes, 'dtl_cost':int(dtl_cost)}
    return ratios, int(dtl_cost), len(tree)

def assess_tree_balance_mad(folder):
    if not os.path.isfile('%s/%s-MAD.ranger_input' %(folder, folder)):
        return []

    tree = ete3.Tree(open('%s/%s-MAD.ranger_input' %(folder, folder)).readlines()[1])
    children = sorted(tree.children, key=len)

    return len(children[0])/float(len(children[1]))

group_list = os.listdir('.')
to_remove  = []
for n in group_list:
    if not os.path.isfile('%s/%s.optResolution1.ranger_input' %(n, n)):
        to_remove.append(n)

for n in to_remove:
    group_list.remove(n)

pool = multiprocessing.Pool(processes=6)
tmp_balance = pool.map(assess_tree_balance_mad, group_list)
mad_balance = []
for n in tmp_balance:
    if n:
        mad_balance.append(n)

pool = multiprocessing.Pool(processes=4)
tmp_balance = pool.map(assess_tree_balance_ranger, group_list)

ranger_balance  = []
dtl_costs       = []
transfer_rate   = []
for group, data in zip(group_list, tmp_balance):
    ranger_balance.extend(data['balance_ratios'])

    tmp_rates = np.asarray(data['num_transfers'])/float(data['num_internal_nodes'])
    transfer_rate.extend(tmp_rates.tolist())

    tmp_dtl_ratio = data['dtl_cost']/float(data['num_internal_nodes'])
    dtl_costs.extend([tmp_dtl_ratio]*len(data['balance_ratios']))

chart = sns.jointplot(x=pd.Series(name='Balance ratio', data=ranger_balance),
                      y=pd.Series(name='Transfer rate', data=transfer_rate),
                      joint_kws=dict(alpha=0.3), marginal_kws=dict(kde=True), s=40, edgecolor='w')
chart.fig.tight_layout()
chart.savefig('../balance_VS_transfer_rate.pdf', dpi=300)

chart = sns.jointplot(x=pd.Series(name='Balance ratio', data=ranger_balance),
                      y=pd.Series(name='DTL cost', data=dtl_costs),
                      joint_kws=dict(alpha=0.3), marginal_kws=dict(kde=True), s=40, edgecolor='w', color='green')
chart.fig.tight_layout()
chart.savefig('../balance_VS_dtl.pdf', dpi=300)

for a, (n,m,k) in zip(group_list, tmp_balance):
    for i in n:
        ranger_balance.append(i)
        dtl_cost.append(m/float(k))

fig, ax = plt.subplots()
sns.kdeplot(compare_balance['with_transfer_on_root'], color='blue',     shade=True, ax=ax, label='Root labeled as transfer')
sns.kdeplot(compare_balance['without_transfer_on_root'], color='green', shade=True, ax=ax, label='Root labeled as speciation')
ax.set_xlabel('Normalized reconciliation costs')
ax.set_ylabel('Density')
fig.set_size_inches(15,8)
fig.tight_layout()
fig.savefig('../dtl_costs-transfer_on_root.pdf', dpi=300)

fig, axs = plt.subplots(nrows=2, sharex=True)
sns.kdeplot(ranger_balance, shade=True, color='blue', ax=axs[0])
axs[0].set_title('RANGER topologies tree balance')
axs[0].set_ylabel('Frequency')
sns.kdeplot(mad_balance, shade=True, color='green', ax=axs[1])
axs[1].set_title('MAD topologies tree balance')
axs[1].set_xlabel('root children balance distribution')
axs[1].set_ylabel('Frequency')
fig.tight_layout()
fig.savefig('../balance_distribution.pdf', dpi=300)

########################################################################################################################
#                                                                                                                      #
# Duplication Transfer ratio                                                                                           #
########################################################################################################################
def assess_balance_DT_ratio(folder):
    ratios    = []
    dt_ratios = []
    trees     = []
    count     = 1
    while os.path.isfile('%s/%s.optResolution%i.ranger_input' %(folder, folder, count)):
        trees.append(ete3.Tree(open('%s/%s.optResolution%i.ranger_input' %(folder, folder, count)).readlines()[1]))
        dtl_cost, d, t = re.match('The minimum reconciliation cost is: (\d+) \(Duplications: (\d+), Transfers: (\d+)', open('%s/%s.reconciliation%i' %(folder, folder, count)).readlines()[-3]).groups()
        dt_ratios.append(float(d)/float(t))
        count += 1

    for tree in trees:
        children = sorted(tree.children, key=len)
        ratios.append(len(children[0])/float(len(children[1])))

    return ratios, dt_ratios, int(dtl_cost), len(tree)

pool = multiprocessing.Pool(processes=8)
tmp_balance = pool.map(assess_balance_DT_ratio, group_list)
ranger_balance  = []
dtl_costs       = []
dt_ratios       = []
for group, (balance_ratio, dt_ratio, dtl_cost, num_leaves) in zip(group_list, tmp_balance):
    dt_ratios.append(np.mean(dt_ratio))
    for pos in range(len(balance_ratio)):
        ranger_balance.append(balance_ratio[pos])
        dt_ratios.append(dt_ratio[pos])
        dtl_costs.append(dtl_cost)

fig, ax = plt.subplots()
sns.kdeplot(dt_ratios, color='blue', shade=True, ax=ax)
ax.set_xlabel('Mean family DT ratio')
fig.set_size_inches(15,8)
fig.tight_layout()
fig.savefig('../tmp_ratio.pdf', dpi=300)

chart = sns.jointplot(x=pd.Series(name='Balance ratio', data=dtl_costs),
                      y=pd.Series(name='DT ratio', data=dt_ratios))
chart.fig.tight_layout()
chart.savefig('../tmp_jointplot.pdf', dpi=300)
########################################################################################################################
#                                                                                                                      #
# ASSESS ALPHA AND CYANO SEPARATION                                                                                    #
########################################################################################################################
os.chdir('/work/Alphas_and_Cyanos/reconciliations')

reference_tree = ete3.Tree('/work/Alphas_and_Cyanos/rooted_partitions-with_named_branches.treefile', format=1)

ncbi = ete3.NCBITaxa()
header = 'assembly_accession bioproject biosample wgs_master refseq_category taxid species_taxid organism_name infraspecific_name isolate version_status assembly_level release_type genome_rep seq_rel_date asm_name submitter gbrs_paired_asm paired_asm_comp ftp_path excluded_from_refseq relation_to_type_material'.split()
genbank_summary                     = pd.read_table('/work/assembly_summary_genbank.txt', comment='#', header=None, names=header, dtype={'taxid':str, 'infraspecific_name':str})
genbank_summary['refseq_category']  = genbank_summary['refseq_category'].str.lower()
genbank_summary['assembly_level']   = genbank_summary['assembly_level'].str.lower()
genbank_summary['genome_rep']       = genbank_summary['genome_rep'].str.lower()
genbank_summary.set_index('assembly_accession', inplace=True)
genbank_summary.index               = [re.sub('\.\d+$', '', index).replace('_', '') for index in genbank_summary.index]
genbank_summary                     = genbank_summary.reindex(labels=reference_tree.get_leaf_names(), axis='index')

lineages = pd.DataFrame.from_dict({taxid:{j: i for i, j in ncbi.get_rank(ncbi.get_lineage(taxid)).items()} for taxid in genbank_summary.taxid}).T

def assess_separation(folder):
    if not os.path.isfile('%s/%s-MAD.ranger_input' %(folder, folder)):
        return 0

    tree   = ete3.Tree(open('%s/%s-MAD.ranger_input' %(folder, folder)).readlines()[1])
    taxids = genbank_summary.loc[set([leaf.split('_')[0] for leaf in tree.get_leaf_names()]), 'taxid']
    phyla  = lineages.loc[taxids, 'phylum'].values.astype('int').astype('str').tolist()

    if not len(set(phyla)) == 2 or not set(['1117', '1224']).issubset(phyla):
        return 0

    child1_phyla = set(lineages.loc[genbank_summary.loc[set([leaf.split('_')[0] for leaf in tree.children[0].get_leaf_names()]), 'taxid'],
                                    'phylum'].values.astype('int').astype('str').tolist())
    child2_phyla = set(lineages.loc[genbank_summary.loc[set([leaf.split('_')[0] for leaf in tree.children[1].get_leaf_names()]), 'taxid'],
                                    'phylum'].values.astype('int').astype('str').tolist())

    if len(child1_phyla) == len(child2_phyla) == 1:
        return 1
    else:
        return 0

pool = multiprocessing.Pool(processes=8)
tmp_separation = pool.map(assess_separation, os.listdir('.'))

def assess_separation_ranger(folder):
    trees = []
    count = 1
    while os.path.isfile('%s/%s.optResolution%i.ranger_input' %(folder, folder, count)):
        trees.append(ete3.Tree(open('%s/%s.optResolution%i.ranger_input' %(folder, folder, count)).readlines()[1]))
        count += 1

    if not trees:
        return 0

    for tree in trees:
        taxids = genbank_summary.loc[set([leaf.split('_')[0] for leaf in tree.get_leaf_names()]), 'taxid']
        phyla  = lineages.loc[taxids, 'phylum'].values.astype('int').astype('str').tolist()

        if not len(set(phyla)) == 2 or not set(['1117', '1224']).issubset(phyla):
            continue

        child1_phyla = set(lineages.loc[genbank_summary.loc[set([leaf.split('_')[0] for leaf in tree.children[0].get_leaf_names()]), 'taxid'],
                                    'phylum'].values.astype('int').astype('str').tolist())
        child2_phyla = set(lineages.loc[genbank_summary.loc[set([leaf.split('_')[0] for leaf in tree.children[1].get_leaf_names()]), 'taxid'],
                                    'phylum'].values.astype('int').astype('str').tolist())

        if len(child1_phyla) == len(child2_phyla) == 1:
            return 1

    return 0

pool = multiprocessing.Pool(processes=8)
ranger_separation = pool.map(assess_separation_ranger, os.listdir('.'))
