#!/usr/bin/env python
# coding: utf-8

from Bio import AlignIO, SeqIO, Align, Alphabet
import pandas as pd
import os, re, sys
from copy import deepcopy

present_genes = []
for aln in os.listdir('.'):
    if not aln.endswith('.aln') or 'Mito_' in aln:
        continue
    present_genes.append(aln.replace('.aln', ''))

for aln in os.listdir('/work/Cyano_Clock_alignments/'):
    if not aln.startswith('Analysis14') or not aln.endswith('.fas.fas'):
        continue
    
    gene_name = re.search('Analysis14_(.*).fas.fas', aln).group(1)
    if gene_name in present_genes:
        continue
    
    os.system('cp %s/%s ./%s.aln' % ('/work/Cyano_Clock_alignments',
                                 aln,
                                 gene_name)
             )


aln_alphabet = Alphabet.Gapped(Alphabet.IUPAC.ambiguous_dna)

genomes    = {}
for aln in os.listdir('.'):
    if not aln.endswith('.aln') or 'Mito_' in aln:
        continue
    alignment    = AlignIO.read(aln, 'fasta')
    genomes[aln] = set()
    for entry in alignment:
        if '[' in entry.description and ']' in entry.description:
            genome = re.search('\[(.*)\]$', entry.description, re.M).group(1).replace(' ', '_')
            genome += '_mitochondria'
        else:
            genome = entry.name

        if genome in genomes[aln]:
            #print ('\t**Error, duplicated genome in %s: %s' %(aln, genome))
            sys.exit('\t**Error, duplicated genome in %s: %s' %(aln, genome))

        genomes[aln].add(genome)


genome_union  = set.union(*genomes.values())
missing_genes = {} # just to keep track of the number of missing marker genes in each genome
concatenation = {}
for genome in genome_union:
    missing_genes[genome]             = 0
    concatenation[genome]             = Align.SeqRecord( Align.Seq('', aln_alphabet) )
    concatenation[genome].name        = genome
    concatenation[genome].id          = genome
    concatenation[genome].description = genome


# fill the handles with the marker sequences from each genome
total_genes      = 0.0 # keep track of the number of genes added to the concatenation
current_position = 1
partitions       = open('concatenated_partitions', 'w')
for aln in os.listdir('.'):
    if not aln.endswith('.aln') or 'Mito_' in aln:
        continue
    
    print(aln.replace('.aln', ''))
    tmp_aln      = AlignIO.read(aln, 'fasta' )
    aln_length   = tmp_aln.get_alignment_length() # get the expected size of the alignment so you can compare if all have the same size
    total_genes += aln_length

    for entry in tmp_aln:
        # if this alignment has a different size from the rest, something is reaaaaaly wrong!
        if len(entry) != aln_length:
            sys.exit('\t**Error, block "%s" has a different length than the rest of the MSA: %s' %(entry.name, aln))

        if '[' in entry.description and ']' in entry.description:
            genome = re.search('\[(.*)\]$', entry.description, re.M).group(1).replace(' ', '_')
            genome += '_mitochondria'
        else:
            genome = entry.name
        concatenation[genome] += deepcopy(entry.seq)

    partitions.write('LG, %s = %i-%i\n' %(aln.replace('.aln', ''), current_position, current_position+aln_length-1) )
    current_position += aln_length

    #
    # add gaps for those genomes missing this gene (same size as the expected alignment)
    for genome in genome_union.difference(genomes[aln]):
        concatenation[genome] += Align.Seq( '-' * aln_length, aln_alphabet )
        missing_genes[genome] += aln_length
partitions.close()

# remove genomes missing more than 20% of the marker genes
counter = 0
for genome, num_missing_genes in missing_genes.items():
    if num_missing_genes/total_genes > 0.5:
        print('\t\t**%s: missing %.2f%% from concatenated alignment!' %(genome, (num_missing_genes/total_genes)*100))
        counter += 1

print('%i genomes missing more than 10%%' %counter)


AlignIO.write(Align.MultipleSeqAlignment(concatenation.values() ), 'concatenated_alignment.aln', 'fasta')

