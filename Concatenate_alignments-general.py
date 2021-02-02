#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from Bio import AlignIO, SeqIO, Align, Alphabet
import pandas as pd
import os, re, sys
from copy import deepcopy

aln_extension = '.aln'
get_ipython().run_line_magic('cd', '/base/alignment/folder/')


# In[ ]:


aln_alphabet = Alphabet.Gapped(Alphabet.IUPAC.ambiguous_dna)

genomes    = {}
for aln in os.listdir('.'):
    #
    # check if files have the desired extension
    if not aln.endswith(aln_extension):
        continue
    alignment    = AlignIO.read(aln, 'fasta')
    genomes[aln] = set()
    for entry in alignment:
        genome = entry.name

        if genome in genomes[aln]:
            #print ('\t**Error, duplicated genome in %s: %s' %(aln, genome))
            sys.exit('\t**Error, duplicated genome in %s: %s' %(aln, genome))
            break
        genomes[aln].add(genome)


# In[ ]:


genome_union  = set.union(*genomes.values())
missing_genes = {} # just to keep track of the number of missing marker genes in each genome
concatenation = {}
for genome in genome_union:
    missing_genes[genome]             = 0
    concatenation[genome]             = Align.SeqRecord( Align.Seq('', aln_alphabet) )
    concatenation[genome].name        = genome
    concatenation[genome].id          = genome
    concatenation[genome].description = genome


# In[ ]:


#
# fill the handles with the marker sequences from each genome
total_genes      = 0.0 # keep track of the number of genes added to the concatenation
current_position = 1
partitions       = open('concatenated_partitions', 'w')
for aln in os.listdir('.'):
    if not aln.endswith(aln_extension):
        continue
    
    print(aln.replace(aln_extension, ''))
    tmp_aln      = AlignIO.read(aln, 'fasta' )
    aln_length   = tmp_aln.get_alignment_length() # get the expected size of the alignment so you can 
                                                  #   compare if all have the same size
    total_genes += aln_length

    for entry in tmp_aln:
        # if this alignment has a different size from the rest, something is reaaaaaly wrong!
        if len(entry) != aln_length:
            sys.exit('\t**Error, block "%s" has a different length than the rest of the MSA: %s' % 
                     (entry.name, aln))

        genome = entry.name
        concatenation[genome] += deepcopy(entry.seq)

    partitions.write('LG, %s = %i-%i\n' %
                     (aln.replace(aln_extension, ''), current_position, current_position+aln_length-1) )
    current_position += aln_length

    #
    # add gaps for those genomes missing this gene (same size as the expected alignment)
    for genome in genome_union.difference(genomes[aln]):
        concatenation[genome] += Align.Seq( '-' * aln_length, aln_alphabet )
        missing_genes[genome] += aln_length
partitions.close()


# In[ ]:


#
# remove genomes missing more than 20% of the marker genes
counter = 0
for genome, num_missing_genes in missing_genes.items():
    if num_missing_genes/total_genes > 0.8:
        print('\t\t**%s: missing %.2f%% from concatenated alignment!' %(genome, (num_missing_genes/total_genes)*100))
        counter += 1

print('%i genomes missing more than 10%%' %counter)


# In[ ]:


AlignIO.write(Align.MultipleSeqAlignment(concatenation.values() ), 'concatenated_alignment.aln', 'fasta')

