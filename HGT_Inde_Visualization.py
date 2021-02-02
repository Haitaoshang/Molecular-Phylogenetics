#!/usr/bin/env python2.7
#coding: utf-8

import ete3
import os
import itertools
import re
import pandas as pd
from Bio import Phylo

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

def root_like( ref_tree, tree_to_root ):
    outgroup = ''
    for node in sorted( ref_tree.children, key=len ):
        if node.is_leaf():
            putative_outgroups = tree_to_root.get_leaves_by_name(node.name)
            for leaf in putative_outgroups:
                tree_to_root.set_outgroup( leaf )
                if tree_to_root.get_topology_id() == ref_tree.get_topology_id():
                    outgroup = leaf
                    break
        else:
            outgroup_members = []
            for leaf in node:
                outgroup_members.append( tree_to_root.get_leaves_by_name(leaf.name) )

            for outgroup_combination in itertools.product( *outgroup_members ):
                if tree_to_root.check_monophyly( [x.name for x in outgroup_combination], 'name' )[0]:
                    putative_outgroup = tree_to_root.get_common_ancestor( outgroup_combination )
                    tree_to_root.set_outgroup( putative_outgroup )
                    if tree_to_root.get_topology_id() == ref_tree.get_topology_id():
                        outgroup = putative_outgroup
                        break
        if outgroup:
            return tree_to_root

    if not outgroup:
        return None

os.chdir('/work/Alphas_and_Cyanos/groups_with_nested_transfers')
named_reference_tree = ete3.Tree('../rooted_partitions-with_named_branches.treefile', format=1)
final_transfers = {'003575': ['m83 = LCA[GCA000012905_ABA79637, GCA001459775_CUU43052]: Transfer, Mapping --> n262, Recipient --> n297'],
 '004119': ['m16 = LCA[GCA000019945_ACB79599, GCA900100665_SDG33339]: Transfer, Mapping --> n134, Recipient --> n293']}

genome_table       = pd.read_table('../selected_genomes.tab', index_col=31)
genome_table.index = [index.replace('_', '').split('.')[0] for index in genome_table.index]

reticulation_style = ete3.NodeStyle()
donor_style        = ete3.NodeStyle()
recipient_style    = ete3.NodeStyle()

#reticulation_style["fgcolor"]       = "#0f0f0f"
#reticulation_style["size"]          = 0
reticulation_style["vt_line_color"] = "#ff0000"
reticulation_style["hz_line_color"] = "#ff0000"
reticulation_style["vt_line_width"] = 5
reticulation_style["hz_line_width"] = 5
reticulation_style["vt_line_type"]  = 1
reticulation_style["hz_line_type"]  = 0

donor_style[    'bgcolor']          = 'LightSteelBlue'
recipient_style['bgcolor']          = 'DarkSeaGreen'

homologue_group = '004119'
named_gene_tree = ete3.Tree('%s.optimal_root.tree' %homologue_group, format=1)
gene_tree       = ete3.Tree('../alignments/%s.aln.treefile'      %homologue_group, format=0)
for leaf in gene_tree.get_leaves():
    leaf.name = re.sub('\.\d+', '', leaf.name.replace('GCA_', 'GCA'))

if set(named_gene_tree.get_leaf_names()).issubset(gene_tree.get_leaf_names()):
    gene_tree.prune(named_gene_tree.get_leaf_names())
gene_tree = root_like(named_gene_tree, gene_tree)

for node_name, descendant1, descendant2 in re.findall('^(m\d+)\s=\sLCA\[(\S+),\s(\S+)\]', open('../reconciliations/%s/aggregate.reconciliation' %homologue_group).read(), re.M):
    branch      = gene_tree.get_common_ancestor([descendant1, descendant2])
    branch.name = node_name

for transfer in final_transfers[homologue_group]:
    reticulation_name = transfer.split()[0]
    reticulation      = gene_tree.search_nodes(name=reticulation_name)[0]
    for node in reticulation.traverse():
        if node.is_leaf():
            node.add_feature('genome_name', node.name.split('_')[0])

    donor_name, recipient_name = re.search('Mapping --> (\S+), Recipient --> (\S+)$', transfer, re.M).groups()
    donor_branch               = named_reference_tree.search_nodes(name=donor_name    )[0]
    recipient_branch           = named_reference_tree.search_nodes(name=recipient_name)[0]

    donor_intersection     = [leaf.name for leaf in reticulation.get_leaves() if leaf.genome_name in donor_branch.get_leaf_names()    ]
    recipient_intersection = [leaf.name for leaf in reticulation.get_leaves() if leaf.genome_name in recipient_branch.get_leaf_names()]

    donor_branch_in_reticulation     = reticulation.get_common_ancestor(donor_intersection    )
    recipient_branch_in_reticulation = reticulation.get_common_ancestor(recipient_intersection)


    reticulation.set_style(reticulation_style)
    donor_branch_in_reticulation.set_style(        donor_style)
    recipient_branch_in_reticulation.set_style(recipient_style)

##color palette: blue, red, green , orange, purple, black
#color_palette   = '#2e2be2 #e22e2b #2be22e #f9a82c #8a2be2 #000000'
nameFaces       = {}
nameFaces[1117] = ete3.AttrFace("name", fsize=10, fgcolor='#000000')
nameFaces[1224] = ete3.AttrFace("name", fsize=10, fgcolor='#a20417')
def myLayout(node):
    if node.is_leaf():
        node.img_style['size'] = 0
        genome_name, gene_name = node.name.split('_')
        node.name = genome_table.loc[genome_name, 'organism_name']
        ete3.add_face_to_node(nameFaces[genome_table.loc[genome_name, 'phylum']], node, 2, aligned=True)
    elif node.is_root():
        node.img_style['size'] = 0
    else:
        if node.support >= 90:
            node.img_style['size'] = 5
            node.img_style['fgcolor'] = '#a20417'
        else:
            node.img_style['size'] = 0

treeStyle                    = ete3.TreeStyle()
treeStyle.layout_fn          = myLayout
treeStyle.show_leaf_name     = False
treeStyle.draw_guiding_lines = True
tree_to_plot                 = gene_tree.copy(method='deepcopy')

tree_to_plot.write(outfile='%s.tree' %homologue_group, format=1)
tree_plot = tree_to_plot.render('%s.pdf' %homologue_group, tree_style=treeStyle, dpi=600)
