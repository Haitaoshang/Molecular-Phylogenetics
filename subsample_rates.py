#!/usr/bin/env python2
#coding: utf-8

import ete3
import os
import pandas as pd
from Bio import SeqIO, SearchIO, AlignIO, Align, Alphabet
import subprocess
import numpy as np
from scipy.spatial.distance import squareform
from matplotlib import pyplot as plt
import seaborn as sns

os.chdir('/work/site_rate')
ncbi = ete3.NCBITaxa()
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

#
# already done!
#
header = '''assembly_accession bioproject biosample wgs_master refseq_category taxid species_taxid organism_name 
            infraspecific_name isolate version_status assembly_level release_type genome_rep seq_rel_date asm_name 
            submitter gbrs_paired_asm paired_asm_comp ftp_path excluded_from_refseq relation_to_type_material'''.split()
assembly_summary_outgroup                    = pd.read_table('outgroups.tab', comment='#', header=None, names=header,
                                                             dtype={'taxid':str, 'infraspecific_name':str})
assembly_summary_outgroup['refseq_category'] = assembly_summary_outgroup['refseq_category'].str.lower()
assembly_summary_outgroup['assembly_level']  = assembly_summary_outgroup['assembly_level'].str.lower()
assembly_summary_outgroup['genome_rep']      = assembly_summary_outgroup['genome_rep'].str.lower()
assembly_summary_outgroup.set_index('assembly_accession', inplace=True)

sampled_genomes = pd.read_table('sampled_genomes.tab', index_col=0)

tree = ete3.Tree('engaging.treefile', format=1)
if tree.check_monophyly(assembly_summary_outgroup.index, 'name'):
    outgroup = tree.get_common_ancestor(assembly_summary_outgroup.index.tolist())
tree.set_outgroup(outgroup)

sampled_genomes = sampled_genomes.reindex(tree.get_leaf_names())

lineages = {}
for index, row in sampled_genomes.iterrows():
    tmp_lineage = {j: int(i) for i, j in ncbi.get_rank(ncbi.get_lineage(row.taxid)).items()}
    lineages[index] = tmp_lineage

lineages = pd.DataFrame.from_dict(lineages, dtype=int).T
lineages.drop('no rank', axis=1, inplace=True)
lineages.dropna(how='any', subset=['genus'], inplace=True)

archaea_groups = {28889: u'Crenarchaeota', 28890: u'Euryarchaeota', 651137: u'Thaumarchaeota'}

########################################################################################################################
# parse site rates
#
########################################################################################################################
alignment = AlignIO.read('ribosomal_concat.aln', 'fasta')

subprocess.call(['iqtree', '-s', 'ribosomal_concat.aln',
                 '-spp', 'rate_partitions', '-nt', '3', '-redo',
                 '-safe', '-wsr', '-n', '0', '-te', 'partition_file.txt.treefile'])

#
# GAMMA distribution
#
rates = pd.read_table('rate_partitions.rate', comment='#')
simulated_partitions = {}
for partition_number in rates.Part.unique():
    partition                              = rates[rates.Part == partition_number]
    simulated_partitions[partition_number] = {category:{} for category in partition.Cat.unique()}

    for category in partition.Cat.unique():
        sites                                            = partition[partition.Cat == category]
        simulated_partitions[partition_number][category] = {block.name:[] for block in alignment}
        for sequence in alignment:
            simulated_partitions[partition_number][category][sequence.name].append(''.join([sequence[position] for position in sites.index]))

full_sequences = {}
for category in range(13):
    full_sequences[category] = {}
    for sequence in alignment:
        full_sequences[category][sequence.name] = ''

for partition_number, categories in simulated_partitions.items():
    print partition_number
    partition = rates[rates.Part == partition_number]

    for category, sequences in categories.items():
        for header,sequence in sequences.items():
            if header in 'GCF_001315945.1 GCF_001316065.1'.split():
                continue
            full_sequence = ''
            while len(full_sequence) <= partition.shape[0]:
                full_sequence += sequence[0]
            full_sequences[category][header] += full_sequence[:partition.shape[0]]

sorted_taxa = []
for category in range(1, 13):
    out = open('rate_categories/%i.aln' %category, 'wb')
    for sequence in alignment:
        if sequence.name in 'GCF_001315945.1 GCF_001316065.1'.split():
            continue
        if category == 1:
            sorted_taxa.append(sequence.name)
        out.write('>%s\n%s\n' %(sequence.name, full_sequences[category][sequence.name]))
    out.close()

ml_distances = pd.read_table('partition_file.txt.mldist', index_col=0, header=None, skiprows=1, sep=' ')
ml_distances.drop(132, axis='columns', inplace=True)
ml_distances.columns = ml_distances.index

fig, ax = plt.subplots()
with cd('rate_categories'):
    for category in range(1, 13):
        subprocess.call(['distmat', '-sequence', '%i.aln' % category, '-protmethod', '0', '-outfile', '%i.distmat' % category])
        uncorrected_distances = pd.read_table('%i.distmat' % category, skiprows=7, header=None, index_col=-1)
        uncorrected_distances.drop([0, 132], axis='columns', inplace=True)
        uncorrected_distances.columns = uncorrected_distances.index = sorted_taxa

        lower_triangle_indexes = np.tril_indices(uncorrected_distances.shape[0], -1)
        uncorrected_distances.values[lower_triangle_indexes] = uncorrected_distances.T.values[lower_triangle_indexes]

        uncorrected_distances = uncorrected_distances.reindex(index=ml_distances.index,
                                                              columns=ml_distances.columns,
                                                              tolerance=0)

        plot = sns.scatterplot(squareform(ml_distances.values),
                        squareform(uncorrected_distances.values),
                        ax=ax, label='category %i' % category, alpha=0.7, s=5
                        )

ax.set_xlabel('ML pairwise distances from original alignment')
ax.set_ylabel('Uncorrected distances from rate-category sites')

fig.set_size_inches(15,15)
fig.tight_layout()
fig.savefig('saturation_test-combined.pdf', dpi=300)
plt.close()

colors = '#29bece #bcbc35 #7f7f7f #e17ac1 #8b564c #936abb #d42a2f #339f34 #fd7f28 #2678b2 #29bece #bcbc35 '.split()
colors.reverse()
fig, axs = plt.subplots(nrows=3, ncols=4, sharex=True, sharey=True)
with cd('rate_categories'):
    row = col = 0
    for category in range(1, 13):
        if category  in [1,2,3,4]:
            row = 0
        elif category in [5,6,7,8]:
            row = 1
        else:
            row = 2

        if category in [1,5,9]:
            col = 0
        elif category in [2,6,10]:
            col = 1
        elif category in [3,7,11]:
            col = 2
        else:
            col = 3

        #subprocess.call(['distmat', '-sequence', '%i.aln' % category, '-protmethod', '0', '-outfile', '%i.distmat' % category])
        uncorrected_distances = pd.read_table('%i.distmat' % category, skiprows=7, header=None, index_col=-1)
        uncorrected_distances.drop([0, 132], axis='columns', inplace=True)
        uncorrected_distances.columns = uncorrected_distances.index = sorted_taxa

        lower_triangle_indexes = np.tril_indices(uncorrected_distances.shape[0], -1)
        uncorrected_distances.values[lower_triangle_indexes] = uncorrected_distances.T.values[lower_triangle_indexes]

        uncorrected_distances = uncorrected_distances.reindex(index=ml_distances.index,
                                                              columns=ml_distances.columns,
                                                              tolerance=0)

        plot = sns.scatterplot(squareform(ml_distances.values),
                               squareform(uncorrected_distances.values),
                               ax=axs[row, col], label='category %i' % category, alpha=0.5,
                               color=colors[category-1]
                               )

commom_area = fig.add_subplot(111, frameon=False)
commom_area.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
commom_area.set_xlabel('ML pairwise distances from original alignment')
commom_area.set_ylabel('Uncorrected distances from rate-category sites')

fig.set_size_inches(18,14)
fig.tight_layout()
fig.savefig('saturation_test.pdf', dpi=300)
plt.close()
