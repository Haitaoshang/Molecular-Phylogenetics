#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import ete3
import re
import itertools
import multiprocessing
import random

import pandas as pd
import numpy  as np
import igraph as ig
import pickle as pkl

from scipy.spatial.distance import squareform, pdist
from scipy.stats            import mannwhitneyu
from collections            import Counter

ncbi = ete3.NCBITaxa()
get_ipython().run_line_magic('cd', '/work/eggNOG/')


# In[ ]:


sampled_genomes = pd.read_csv('../genomes.tab',
                              sep='\t',
                              index_col=0)


# In[ ]:


lineages = pd.DataFrame()
for taxid in sampled_genomes.species_taxid.unique():
    if pd.isna(taxid):
        continue
    lineages = lineages.append({tax_rank: tmp_taxid 
                                 for tmp_taxid, tax_rank in ncbi.get_rank(ncbi.get_lineage(taxid)).items()},
                                ignore_index=True)
lineages = lineages.reindex(columns=['class', 'family',  'genus', 'phylum',
                                     'order', 'species', 'superkingdom']).copy()
lineages = lineages.query('superkingdom == 2').copy()


# In[ ]:


working_groups  = pd.read_parquet('working_eggNOG_groups.parquet')
working_trees   = pd.read_parquet('working_eggNOG_trees.parquet' )
eggNOG_taxonomy = pd.read_parquet('eggNOG_taxonomy.parquet'      )


# In[ ]:





# In[ ]:


def get_pairwise_distances(group_id):
    
    tree = ete3.Tree(working_trees.loc[group_id, 'tree'])

    leaf_names = []
    for count, node in enumerate(tree.traverse()):
        if node.is_leaf():
            leaf_names.append(node.name)
        else:
            node.name = 'node_%i' % count
    leaf_names = np.array(leaf_names)

    nodes         = []
    children      = []
    branch_length = []
    for node in tree.traverse():
        if not node.is_leaf():
            for child in node.get_children():
                nodes.append(         node.name)
                children.append(     child.name)
                branch_length.append(child.dist)

    branch_length_df                  = pd.DataFrame()
    branch_length_df['node']          = nodes
    branch_length_df['child']         = children
    branch_length_df['branch_length'] = branch_length

    dag  = ig.Graph.TupleList(edges=branch_length_df[['node', 
                                                      'child', 
                                                      'branch_length']].itertuples(index=False), 
                                directed=False, 
                                weights=True)
    
    dist_matrix = pd.DataFrame(index  =leaf_names, 
                               columns=leaf_names, 
                               data   =np.array(dag.shortest_paths(source=leaf_names, 
                                                                   target=leaf_names, 
                                                                   weights='weight'))
                              )
    return(dist_matrix)


# In[ ]:


def create_taxa_graph(dist_matrix, phyla):
    triu_indices       = np.triu_indices_from(dist_matrix, k=1)
    
    edge_list                 = pd.DataFrame()
    edge_list['phylum1']      = phyla[triu_indices[0]]
    edge_list['phylum2']      = phyla[triu_indices[1]]
    edge_list['sequence1']    = dist_matrix.index[triu_indices[0]]
    edge_list['sequence2']    = dist_matrix.index[triu_indices[1]]
    edge_list['distance']     = dist_matrix.values[triu_indices]
    edge_list['inverse_dist'] = np.e**np.negative(edge_list.distance)

    graph  = ig.Graph.TupleList(edges=edge_list[['sequence1', 
                                                 'sequence2', 
                                                 'inverse_dist']].itertuples(index=False), 
                                directed=False, 
                                weights =True)
    
    return(edge_list, graph)


# In[ ]:


def cles(lessers, greaters):
    #
    # https://github.com/ajschumacher/cles/blob/master/cles.py
    #
    """Common-Language Effect Size
    Probability that a random draw from `greater` is in fact greater
    than a random draw from `lesser`.
    Args:
      lesser, greater: Iterables of comparables.
    """
    if len(lessers) == 0 and len(greaters) == 0:
        raise ValueError('At least one argument must be non-empty')
    # These values are a bit arbitrary, but make some sense.
    # (It might be appropriate to warn for these cases.)
    if len(lessers) == 0:
        return 1
    if len(greaters) == 0:
        return 0
    numerator = 0
    lessers, greaters = sorted(lessers), sorted(greaters)
    lesser_index = 0
    for greater in greaters:
        while lesser_index < len(lessers) and lessers[lesser_index] < greater:
            lesser_index += 1
        numerator += lesser_index  # the count less than the greater
    denominator = len(lessers) * len(greaters)
    return float(numerator) / denominator


# In[ ]:


def assess_cluster(reference_phylum, minimal_freq_phyla, cluster_edges, cluster_nodes):

    #
    # store distances between reference phylum and others
    cluster_dists = pd.DataFrame(columns=['phylum', 'median', 'distances'])

    #
    # traverse all phylum pairs containing the reference phylum
    for phylum in minimal_freq_phyla:
        if phylum == reference_phylum:
            continue

        inter_phyla = cluster_edges.loc[((cluster_edges.phylum1==reference_phylum) & (cluster_edges.phylum2==phylum)) |                                        ((cluster_edges.phylum2==reference_phylum) & (cluster_edges.phylum1==phylum))]
        
        #
        # create a quadratic matrix of pairwise distances between phyla
        #   rows and columns must be unique pairs of sequence names
        indices     = np.unique(inter_phyla[['sequence1', 'sequence2']])
        adjacencies = pd.DataFrame(index  =indices, 
                                   columns=indices,
                                   data   =0.0)
        
        #
        # add distances between sequences from the edge list to the quadratic matrix
        #   as the matrix is quadratic, add values to both directions
        indexer     = adjacencies.index.get_indexer
        adjacencies.values[indexer(inter_phyla.sequence1), indexer(inter_phyla.sequence2)] = inter_phyla.distance.values
        adjacencies.values[indexer(inter_phyla.sequence2), indexer(inter_phyla.sequence1)] = inter_phyla.distance.values

        #
        # sum the obtained distances into a single cell
        tmp_closest_to_phylum = adjacencies.loc[cluster_nodes.loc[cluster_nodes.phylum==reference_phylum, 'name'],
                                                cluster_nodes.loc[cluster_nodes.phylum==phylum,           'name']].sum()
        tmp_closest_to_phylum.sort_values(inplace=True)
        #
        # and grab the five sequences from <phylum> closest to all from <reference phylum>
        tmp_closest_to_phylum = tmp_closest_to_phylum.index[:5]

        try:
            #
            # get all inter-phyla distances between
            distances_to_reference_phylum = adjacencies.loc[
                #
                # all sequences from <reference phylum>
                cluster_nodes.loc[cluster_nodes.phylum==reference_phylum, 'name'],
                #
                # the 5 sequences from <phylum> closest to <reference phylum>
                tmp_closest_to_phylum
            ].values.flatten()
        except IndexError:
            continue        

        cluster_dists = cluster_dists.append(pd.Series(data =[phylum, 
                                                              np.median(distances_to_reference_phylum), 
                                                              distances_to_reference_phylum], 
                                                       index=['phylum', 'median', 'distances']),
                                             ignore_index=True)
    return(cluster_dists)


# In[ ]:


def get_phyla_evol_distances(group_id):    
    dist_matrix = get_pairwise_distances(group_id)

    taxids = [int(leaf.split('.')[0]) for leaf in dist_matrix.index]
    phyla  = eggNOG_taxonomy.loc[taxids, 'phylum'].values.astype(int)

    edge_list, graph  = create_taxa_graph(dist_matrix, phyla)

    random.seed(12345)
    clusters = graph.community_multilevel(weights='weight')

    node_data = pd.DataFrame(columns=['name', 'phylum', 'cluster'],
                             data   =zip(dist_matrix.index, 
                                         phyla, 
                                         clusters.membership)
                            )
    
    cluster_evol_relations = {}
    target_phyla = {1090,    # chlorobi
                    1117,    # cyanobacteria
                    1224,    # proteobacteria
                    200795,  # chloroflexi
                    976,     # bacteroidetes
                    1134404, # ignavibacteriae
                    1798710} # melainabacteria

    for cluster_num in set(clusters.membership):
        
        cluster_nodes      = node_data[node_data.cluster==cluster_num]
        minimal_freq_phyla = [phylum 
                              for phylum, frequency in Counter(cluster_nodes.phylum).items() 
                              if frequency>=5   # at least five sequences from a phylum
                              and phylum > 0]   # given pandas manipulation, unknown phyla are represented 
                                                #   as NAN, which when forced to be int are negative numbers
                                                #   that is why we ignore phylum taxids smaller than zero...
       
        #
        # if there are fewer than two phyla of interested within the tree cluster there is no reason to
        #   assess it...
        if len( target_phyla.intersection( minimal_freq_phyla ) ) < 2:
            continue
        
        cluster_evol_relations[cluster_num] = {}
        
        #
        # filter patristic distances from the whole tree to only those between sequences within the cluster
        cluster_edges = edge_list.loc[(edge_list.sequence1.isin(cluster_nodes.name)) &
                                      (edge_list.sequence2.isin(cluster_nodes.name)),
                                      ['phylum1', 'phylum2', 'sequence1', 'sequence2', 'distance']]

        #
        # filter again, leaving only between sequences whose phylum fulfills the minimal frequency
        cluster_edges      = cluster_edges[(cluster_edges.phylum1.isin(minimal_freq_phyla)) &                                           (cluster_edges.phylum2.isin(minimal_freq_phyla))]
        #
        # we will divide all pairwise distances by the cluster's mean
        normalizer         = np.median(cluster_edges.distance)
        #
        # remove all intra-phylum distances
        cluster_edges      = cluster_edges[cluster_edges.phylum1 != cluster_edges.phylum2] 

        #
        # assess all inter-phyla distances, using each phylum as a reference to itself
        #
        for ref_phylum in target_phyla.intersection(minimal_freq_phyla):
            cluster_dists = assess_cluster(ref_phylum, 
                                           minimal_freq_phyla, 
                                           cluster_edges,
                                           cluster_nodes)

            #
            # sort phyla by its avg. distance to <reference phylum>
            cluster_dists.sort_values('median', inplace=True)
            cluster_evol_relations[cluster_num][ref_phylum] = {
                'df'         :cluster_dists[['phylum', 'median']].copy(),
                'significant':False
            }
            if not cluster_dists.shape[0]:
                continue

            cluster_evol_relations[cluster_num][ref_phylum]['df']['median']   /= normalizer
            
            #
            # if there is only one <phylum> together with <reference phylum>, the proximity between
            #   both is automaticaly significant
            if cluster_dists.shape[0] == 1:
                cluster_evol_relations[cluster_num][ref_phylum]['significant'] = True
                continue

            try:
                #
                # test if distances from the closest phylum to <reference phylum> is significantly
                #   smaller than distances from the second closest phylum
                hypothesis = mannwhitneyu(cluster_dists.iloc[0, 2], 
                                          cluster_dists.iloc[1, 2], 
                                          alternative='less')
            except ValueError:
                continue
            else:
                #
                # both cles lines below should work identically... I am leaving the above because it is the one I
                #   used in the paper, no real reason...
                effect_size = hypothesis.statistic / (len(cluster_dists.iloc[0, 2])*len(cluster_dists.iloc[1, 2]))
#                 effect_size = 1-cles(cluster_dists.iloc[0, 2], cluster_dists.iloc[1, 2])

                if hypothesis.pvalue < 0.01 and effect_size < 0.2:
                    cluster_evol_relations[cluster_num][ref_phylum]['significant'] = True
    
    return(group_id, cluster_evol_relations)


# In[ ]:


# %%time
# get_phyla_evol_distances('COG0499')


# In[ ]:


get_ipython().run_cell_magic('time', '', 'pool    = multiprocessing.Pool(processes=10, maxtasksperchild=5)\nresults = pool.map_async(get_phyla_evol_distances, working_groups.group_id.values)\npool.close()\npool.join()')


# In[ ]:


with open('all_results.pkl', 'wb') as out:
    pkl.dump(results.get(), out)
del(results)

