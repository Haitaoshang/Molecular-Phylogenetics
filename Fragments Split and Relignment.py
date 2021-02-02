from Bio import SeqIO
import os
from commands import getoutput
import subprocess
import itertools

os.chdir('/work/partial_transfers')

fragments = {}
for sequence in SeqIO.parse('Val_Leu.fasta', 'fasta'):
    fragments[sequence.id] = []
    for position in range(0, len(sequence), 5):
        tmp_fragment = sequence[position:position+10]
        tmp_fragment.description = '(%i:%i)' %(position, position+10)
        fragments[sequence.id].append(tmp_fragment)

product_count = 0
for product in itertools.product(*fragments.values()):
    SeqIO.write(product, 'fragment_combinations/%i.faa' %product_count, format='fasta' )
    product_count += 1
    break