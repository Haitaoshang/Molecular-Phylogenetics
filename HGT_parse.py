#!/usr/bin/env python

import ete3
import os
import re
import pickle as pkl
import pandas as pd
import numpy as np
import random
from commands import getoutput
import multiprocessing

########################################################################################################################
#
#   Source:
# https://stackoverflow.com/questions/3288595/multiprocessing-how-to-use-pool-map-on-a-function-defined-in-a-class
#
def fun(f, q_in, q_out):
    while True:
        i, x = q_in.get()
        if i is None:
            break
        q_out.put((i, f(x)))

def parmap(f, X, nprocs=multiprocessing.cpu_count()):
    '''Parallel mapping for "self" functions!
    parmap(<function to be prallelized>, <list of input arguments>, <num of threads>)

    '''
    q_in = multiprocessing.Queue(1)
    q_out = multiprocessing.Queue()

    proc = [multiprocessing.Process(target=fun, args=(f, q_in, q_out))
            for _ in range(nprocs)]
    for p in proc:
        p.daemon = True
        p.start()

    sent = [q_in.put((i, x)) for i, x in enumerate(X)]
    [q_in.put((None, None)) for _ in range(nprocs)]
    res = [q_out.get() for _ in range(len(sent))]

    [p.join() for p in proc]

    return [x for i, x in sorted(res)]
########################################################################################################################

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

class hgt:

    def __init__(self,
                 reference_tree_file='/work/Alphas_and_Cyanos/rooted_partitions-with_named_branches.treefile',
                 assembly_summary='/work/Alphas_and_Cyanos/assembly_summary_genbank.txt',
                 output_folder='.',
                 tree_folder='/work/Alphas_and_Cyanos/ranger_input_trees',
                 reconciliations_folder='/work/Alphas_and_Cyanos/reconciliations',
                 reconciliation_sufix='.reconciliation',):

        self.ncbi                   = ete3.NCBITaxa()
        self.named_reference_tree   = ete3.Tree(reference_tree_file, format=1)
        self.output_folder          = output_folder.strip()
        self.tree_folder            = tree_folder.strip()
        self.reconciliations_folder = reconciliations_folder.strip()
        self.reconciliation_sufix   = reconciliation_sufix.strip()

        header = 'assembly_accession bioproject biosample wgs_master refseq_category taxid species_taxid organism_name infraspecific_name isolate version_status assembly_level release_type genome_rep seq_rel_date asm_name submitter gbrs_paired_asm paired_asm_comp ftp_path excluded_from_refseq relation_to_type_material'.split()
        genbank_summary                     = pd.read_table(assembly_summary, comment='#', header=None, names=header, dtype={'taxid':str, 'infraspecific_name':str})
        genbank_summary['refseq_category']  = genbank_summary['refseq_category'].str.lower()
        genbank_summary['assembly_level']   = genbank_summary['assembly_level'].str.lower()
        genbank_summary['genome_rep']       = genbank_summary['genome_rep'].str.lower()
        genbank_summary.set_index('assembly_accession', inplace=True)
        genbank_summary.index               = [re.sub('\.\d+$', '', index).replace('_', '') for index in genbank_summary.index]
        self.assembly_summary               = genbank_summary.reindex(self.named_reference_tree.get_leaf_names())

    def rename_branches(self, reconciliation_file):
        reconciliation = open(reconciliation_file).read()
        tree = ete3.Tree(reconciliation.split('\n')[7], format=1)
        branches = re.findall('^(m\d+) = LCA\[(\S+), (\S+)\]', reconciliation, re.M)
        for node in tree.traverse():
            if node.is_leaf():
                continue
            else:
                node.name = ''
        for name, leaf1, leaf2 in branches:
            node = tree.get_common_ancestor(leaf1, leaf2)
            if node.name:
                print node.get_leaves_by_name()
                break
            node.name = name
        return tree

    def assess_transfers(self, transfer_file, reconciliation_sufix=None):
        if reconciliation_sufix is None:
            reconciliation_sufix = self.reconciliation_sufix

        if not transfer_file.endswith('.pkl'):
            return None

        group         = os.path.basename(transfer_file).replace('.pkl', '')
        transfer_data = pkl.load(open(transfer_file))
        mapped_nodes  = {}

        support_tree = ete3.Tree('%s/%s.tree' %(self.tree_folder, group))

        with cd('%s/%s' %(self.reconciliations_folder, group)):
            reconciliation_count = 1
            descendants = {}
            for donor, recipient, reticulation, topology_id in transfer_data:
                recipient_branch       = self.named_reference_tree.search_nodes(name=recipient)[0]
                descendants[recipient] = [descendant.name for descendant in recipient_branch.get_descendants() if not descendant.is_leaf()]
                descendants[recipient].append(recipient)

            while True:
                tmp_tree = self.rename_branches('{group}{sufix}{counter}'.format(group=group, sufix=reconciliation_sufix, counter=reconciliation_count))
                if tmp_tree.children[0].is_leaf():
                    support_tree.set_outgroup(tmp_tree.children[0].name)
                elif tmp_tree.children[1].is_leaf():
                    support_tree.set_outgroup(tmp_tree.children[1].name)
                else:
                    is_it_monophyletic, clade_type, fucking_up = support_tree.check_monophyly(tmp_tree.children[0].get_leaf_names(), 'name')
                    if is_it_monophyletic:
                        support_tree.set_outgroup(support_tree.get_common_ancestor(tmp_tree.children[0].get_leaf_names()))
                    else:
                        support_tree.set_outgroup(support_tree.get_common_ancestor(tmp_tree.children[1].get_leaf_names()))

                if support_tree.get_topology_id() != tmp_tree.get_topology_id():
                    return {group:[]}
                    #raise ValueError('%s support and named trees not matching: reconciliation_count=%i' %(group, reconciliation_count))

                mapped_donors       = {topology_id:'' for donor, recipient, reticulation, topology_id in transfer_data}
                map_count           = 0
                for node in tmp_tree.traverse():
                    if not node.is_leaf() and node.get_topology_id() in mapped_donors:
                        if support_tree.get_common_ancestor(node.get_leaf_names()).support >= 95:
                            mapped_donors[node.get_topology_id()] = node.name
                        map_count += 1
                    if map_count == len(mapped_donors):
                        break

                for donor, recipient, reticulation, topology_id in transfer_data:
                    index = tuple([donor, recipient, reticulation, topology_id])
                    if not mapped_donors[topology_id]:
                        continue

                    mapped_donor = tmp_tree.search_nodes(name=mapped_donors[topology_id])[0]

                    mapping_child1 = getoutput('grep -E "^%s = LCA\[\S+, \S+\]: \S+, Mapping --> (%s)" %s.reconciliation{%i..%i}' %(mapped_donor.children[0].name, '|'.join(descendants[recipient]), group, reconciliation_count, reconciliation_count+19))
                    mapping_child2 = getoutput('grep -E "^%s = LCA\[\S+, \S+\]: \S+, Mapping --> (%s)" %s.reconciliation{%i..%i}' %(mapped_donor.children[1].name, '|'.join(descendants[recipient]), group, reconciliation_count, reconciliation_count+19))

                    mapping = mapping_child1
                    child   = mapped_donor.children[0]
                    if len(mapping_child1.split('\n')) < len(mapping_child2.split('\n')):
                        mapping = mapping_child2
                        child   = mapped_donor.children[1]

                    if len(mapping.split('\n')) > 20:
                        raise ValueError('something is wrong with group %s, between %s and %s' %(group, donor, recipient))

                    if index not in mapped_nodes:
                        mapped_nodes[index] = []
                    mapped_nodes[index].append((topology_id, child.get_topology_id()))

                reconciliation_count += 20
                if not os.path.isfile('{group}{sufix}{counter}'.format(group=group, sufix=reconciliation_sufix, counter=reconciliation_count)):
                    break

        selected_transfers = []
        for index, mappings in mapped_nodes.items():
            if len(set(mappings)) != 1:
                break
                raise ValueError('Multiple reconciliations mappings for group %s, between %s and %s' %(group, index[0], index[1]))
            else:
                selected_transfers.append({'donor':index[0], 'recipient':index[1], 'reticulation':index[2], 'donor_map':mappings[0][0], 'recipient_map':mappings[0][1]})

        return {group:selected_transfers}

    def visualize_tree(self, group, transfers, output_folder=None):
        if output_folder is None:
            output_folder = self.output_folder

        tmp_names = {}

        gene_tree = ete3.Tree('{tree_folder}/{group}.tree'.format(tree_folder=self.tree_folder, group=group), format=0)
        for leaf in gene_tree.get_leaves():
            genome_name, gene_name = leaf.name.split('_')
            leaf.add_feature('genome', genome_name)

        count = 0
        for mappings in transfers:
            is_it_monophyletic, clade_type, fucking_up = gene_tree.check_monophyly(mappings['reticulation'].get_leaf_names(), 'name', unrooted=False)
            if is_it_monophyletic:
                donor_map = gene_tree.get_common_ancestor(mappings['reticulation'].get_leaf_names())
            else:
                gene_tree.set_outgroup(fucking_up.pop())
                donor_map = gene_tree.get_common_ancestor(mappings['reticulation'].get_leaf_names())

            if donor_map.get_topology_id() != mappings['donor_map']:
                raise ValueError('Topology ids not matching, {group} between {donor} and {recipient}' .format(group=group, donor=mappings['donor'], recipient=mappings['recipient']))

            for node in donor_map.children:
                if not node.is_leaf() and node.get_topology_id() == mappings['recipient_map']:
                    recipient_map = node
                    break

            #
            # set donor name
            tmp_id                            = str(random.random())
            tmp_names[tmp_id]                 = '[&support=%.2f,hgt_role="donor%i",hgt_count=%i]' %(donor_map.support, count, count)
            donor_map.name = tmp_id

            #
            # set recipient name
            tmp_id                                 = str(random.random())
            tmp_names[tmp_id]                      = '[&support=%.2f,hgt_role="recipient%i",hgt_count=%i]' %(recipient_map.support, count, count)
            recipient_map.name  = tmp_id

            count += 1
            del(donor_map, recipient_map)

        out  = open('%s/%s.Figtree.tree' %(output_folder, group), 'wb')
        out.write("#NEXUS\nbegin taxa;\n\tdimensions ntax=%i;\n\ttaxlabels\n" %len(gene_tree))
        for node in gene_tree.traverse():
            if node.is_leaf():
                taxid         = self.assembly_summary.loc[node.genome, 'taxid']
                lineage       = {j: i for i, j in self.ncbi.get_rank(self.ncbi.get_lineage(taxid)).items()}
                lineage_names = self.ncbi.get_taxid_translator(lineage.values())

                out.write('\t%s ' %(node.name))
                comment = ['source_name="%s"' %self.ncbi.get_taxid_translator([taxid]).itervalues().next()]
                for rank in ['phylum', 'class', 'order', 'family', 'genus']:
                    if rank in lineage:
                        comment.append('tax_%s="%s"' %(rank, lineage_names[lineage[rank]]))
                out.write('[&%s]\n' %' '.join(comment))

            else:
                if node.support and not node.name:
                    node.name = '[&support=%.2f]' %node.support

        rooted_tree = ete3.Tree(open('{folder}/{group}/{group}{sufix}1'.format(folder=self.reconciliations_folder, sufix=self.reconciliation_sufix, group=group)).readlines()[7], format=1)
        if rooted_tree.children[0].is_leaf():
            gene_tree.set_outgroup(rooted_tree.children[0].name)
        elif rooted_tree.children[1].is_leaf():
            gene_tree.set_outgroup(rooted_tree.children[1].name)
        else:
            is_it_monophyletic, clade_type, fucking_up = gene_tree.check_monophyly(rooted_tree.children[0].get_leaf_names(), 'name')
            if is_it_monophyletic:
                gene_tree.set_outgroup(gene_tree.get_common_ancestor(rooted_tree.children[0].get_leaf_names()))
            else:
                gene_tree.set_outgroup(gene_tree.get_common_ancestor(rooted_tree.children[1].get_leaf_names()))

        newick_text = gene_tree.write(format=1)
        newick_text = re.sub('_&support_(\d+\.\d\d)_', '[&support=\\1]', newick_text)
        for key, value in tmp_names.items():
            newick_text = newick_text.replace(key, value)
        out.write(';\nend;\n')
        out.write('begin trees;\n\ttree tree_1 = [&R] %s\nend;' %newick_text)
        out.close()
