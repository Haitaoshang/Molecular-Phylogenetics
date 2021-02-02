import subprocess
import os
import ete3
import random
import numpy as np
from Bio import SeqIO
import pandas as pd
import multiprocessing
from itertools import product


#
# initial definitions
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


os.chdir('/work/site_rate/sequence_simulation')
random.seed(12345)
num_replicates      = 100
sequence_length     = 1000
num_threads         = 3
num_rate_categories = 12

indelible_conf = '''\
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
//  INDELible V1.03 control file - site-rate project                               //
//                                                                                 //
//      Automaticaly generated, more information:                                  //
//          github.com/lthiberiol/fournierLab/tree/master/site_rates               //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////

[TYPE]  AMINOACID 1

[SETTINGS]
    [randomseed]   12345

[MODEL] model1  [submodel]  LG         //   LG
                [rates]     0 1.3 20    //   pinv=0, alpha=1.3, discrete gamma categories=20
                [statefreq]
                    0.083   0.078   0.035   0.059   // list of 20 numbers
                    0.006   0.031   0.082   0.080   // A R N D C
                    0.018   0.063   0.070   0.078   // Q E G H I
                    0.025   0.025   0.046   0.050   // L K M F P
                    0.054   0.008   0.025   0.084   // S T W Y V

[TREE] T1  {t1}
[TREE] T2  {t2}
[TREE] T3  {t3}

//[PARTITIONS] partition1 [T1 model1 {length}]   // tree T1, model model1, root length 1000
[PARTITIONS] partition2 [T2 model1 {length}]   // tree T2, model model1, root length 1000
[PARTITIONS] partition3 [T3 model1 {length}]   // tree T3, model model1, root length 1000

[EVOLVE]    //partition1  {num_replicates}   {t1_name}
            partition2  {num_replicates}   {t2_name}
            partition3  {num_replicates}   {t3_name}
'''

partition_name = 'base_tree'
trees                 = {}
trees[partition_name] = ete3.Tree('../rate_categories/species.tree')

branch_lengths = np.asarray([node.dist for node in trees['base_tree'].traverse() if not node.is_root()])
branch_lengths.sort()
trees['short2long'] = trees['base_tree'].copy()
for position, node in enumerate(trees['short2long'].traverse()):
    if node.is_root():
        continue
    node.dist = branch_lengths[position-1]

reverse_order       = list(reversed(branch_lengths))
trees['long2short'] = trees['base_tree'].copy()
for position, node in enumerate(trees['long2short'].traverse()):
    if node.is_root():
        continue
    node.dist = reverse_order[position-1]

out = open('control.txt', 'w')
out.write(indelible_conf.format(t1=trees['base_tree'].write(format=5),  t1_name='base_tree',
                                t2=trees['short2long'].write(format=5), t2_name='short2long',
                                t3=trees['long2short'].write(format=5), t3_name='long2short',
                                num_replicates=num_replicates, length=sequence_length))
out.close()

subprocess.call(['/work/site_rate/indelible/INDELibleV1.03/bin/indelible_1.03_OSX_intel'])

for partition_name, tree in trees.items():
    fasta = open('%s.fas' % partition_name).read().strip()

    if not os.path.isdir(partition_name):
        os.mkdir(partition_name)
    else:
        os.system('rm %s/*' % partition_name)

    tree.write(outfile='%s/reference.tre' % partition_name, format=5)

    for count, block in enumerate(fasta.split('\n     \n')):
        out = open('%s/%i.fas' % (partition_name, count+1), 'w')
        out.write(block)
        out.close()

#
# classify sites into rate-categories
#
def write_rates((partition_name, replicate_number)):
    subprocess.call(['iqtree', '-s', '%s/%i.fas' % (partition_name, replicate_number), '-m', 'LG+G12', '-redo',
                     '-safe', '-wsr', '-nt', '1', '-n', '0', '-pre', '%s/%i' % (partition_name, replicate_number),
                     '-te', '%s/reference.tre' % partition_name, '-quiet'])


pool = multiprocessing.Pool(processes=num_threads)
pool.map(write_rates, product(trees.keys(), range(1, num_replicates+1)))

#
# parse rates classification
#
for partition_name in trees.keys():
    print partition_name

    with cd(partition_name):
        if not os.path.isdir('categories'):
            os.mkdir('categories')
        else:
            os.system('rm -r categories/*')

        for replicate in range(1, num_replicates+1):
            alignment = list(SeqIO.parse('%i.fas' % replicate, 'fasta'))

            ratesG = pd.read_table('%i.rate' % replicate, comment='#')
            for category in ratesG.Cat.unique():
                site_df        = ratesG[ratesG.Cat == category]
                category_aln   = {sequence.name:'' for sequence in alignment}
                for sequence in alignment:
                    category_aln[sequence.name] = ''.join([sequence[position] for position in site_df.index])

                out = open('categories/%i.%i.aln' % (replicate, category), 'w')
                for header, sequence in category_aln.items():
                    full_sequence = ''
                    while len(full_sequence) <= sequence_length:
                        full_sequence += sequence
                    out.write('>%s\n%s\n' % (header, full_sequence[:1000]))
                out.close()


def run_bootstrap((replicate_number, category)):
    subprocess.call(['iqtree', '-s', '%i.%i.aln' % (replicate_number, category), '-m', 'LG+G1', '-redo',
                     '-safe', '-nt', '1', '-pre', '%i.%i' % (replicate_number, category),
                     '-alrt', '1000', '-keep-ident', '-quiet', '-te', '../reference.tre'])


for partition_name in trees.keys():
    with cd('%s/categories' % partition_name):
        pool = multiprocessing.Pool(processes=num_threads)
        pool.map(run_bootstrap, product(range(1, num_replicates+1), range(1, num_rate_categories+1)))
