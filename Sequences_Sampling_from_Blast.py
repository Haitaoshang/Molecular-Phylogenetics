#!/usr/bin/env python
# coding: utf-8

# # Sample phylum representative sequences from an all *VS.* all BLAST output
# 
# here we have a some generally recruited sequences which we previously scanned using HMMER, and now we have to reduce the number of sequences to a workable amount.
# 
# The pairwise matrix resulted from **BLAST** it is likely waaaay to big to work as a regular $m X m$ matrix, so we will use graph structures too simplify the analysis of the data.
# 
# > I have been working in a 100G cluster node and trying to use a regular $m X m$ matrices broke the job *every time*...
# 
# A Jupyter Notebook version of this code is available at my [github page](https://github.com/lthiberiol/jupyter_notebooks/blob/master/kelsey/sample_sequences_from_blast.ipynb). 

# In[1]:


import pandas as pd
import plotly
import chart_studio.plotly as ptl
import plotly.graph_objects as go
import networkx as nx
import multiprocessing
import subprocess
import ete3
import getpass

ncbi = ete3.NCBITaxa()
username = input('Plotly username:')
api_key  = getpass.getpass('Plotly API key:')
ptl.sign_in(username, api_key)


# ## Input data
# 
# the input is a regular all *VS.* all BLAST output in tabular format (either `-outfmt 7` or `-outfmt 8`)

# In[2]:


group = '1563_AAM72826.1-GCA_000006985.1'

blast_search = pd.read_csv('recruited_fastas/%s.bls' %group,
                           sep    ='\t',
                           comment='#',
                           header =None,
                           usecols=[0, 1, 10, 11],
                           names  =['query_acc',
                                    'subject_acc',
                                    'evalue',
                                    'bit_score'])


# ## General BLAST treatment
# 
# - remove BLAST alignments that aren't bi-directional, and remove multiple hits against the same subject (keeping only the 1st one)

# In[3]:


unique_queries = blast_search.query_acc.unique()
blast_search   = blast_search.query('subject_acc.isin(@unique_queries)').copy()

blast_search.drop_duplicates(subset=['query_acc', 'subject_acc'], inplace=True)


# \
# get bit scores for sequence alignments against itself (*important later to calculate normalized bit scores*)

# In[4]:


self_bit_scores = blast_search.loc[blast_search.query_acc==blast_search.subject_acc, 
                                   ['query_acc', 'bit_score']].copy()
self_bit_scores.set_index('query_acc', inplace=True)
self_bit_scores = self_bit_scores.bit_score


# \
# remove queries we could't find and alignment against itself (*shouldn't be frequent, but I have found a few cases...*)

# In[5]:


blast_search = blast_search.query("query_acc.isin(@self_bit_scores.index)").copy()


# ## Calculate *bit score Index*
# (it is just a normalization of bit score described [here](https://rdrr.io/cran/micropan/man/bDist.html))
# 
# first we have to add extra columns to our BLAST DataFrame:
# - one column for the query bit score against itself
# - another regarding the subject's bit score against itself

# In[6]:


blast_search['query_self_bscore'] = self_bit_scores[blast_search.query_acc].values
blast_search['sbjct_self_bscore'] = self_bit_scores[blast_search.subject_acc].values


# \
# now we add the **bit score Index** column:
# 
# > $Index_{BS}=\frac{2 BS_{ij}}{BS_{i}+BS_{j}}$\
# > where $BS_{i}$ is the query's self bit score, $BS_{j}$ is the subjects self bit score, and $BS_{ij}$ is the query *vs* the subject bit score

# In[7]:


blast_search['bscore_index'] = (blast_search.bit_score*2)/                               (blast_search.query_self_bscore+blast_search.sbjct_self_bscore)


# \
# simple cleanup, let's remove alignments for which we failed to calculate $Index_{BS}$, shouldn't happen but it will avoid errors latter

# In[8]:


blast_search.dropna(subset=['bscore_index'], inplace=True)


# ## Remove redundant sequences
# 
# create a sub-matrix containing only alignments between extremely similar sequences
# > I used a $Index_{BS}\geq0.95$, but feel free to adjust

# In[9]:


redundant_blast = blast_search[(blast_search.bscore_index > 0.95) &
                               (blast_search.query_acc != blast_search.subject_acc)]


# \
# use this sub-matrix to create a graph with edges only between extremely similar proteins
# 
# > The idea here is that this approach will divide our graph in multiple **connected components** each contaning very similar proteins not very likely to add diversity to our protein pool
# >
# > Doing this is more practical and straightforward than performing a regular community analysis

# In[10]:


tmp_graph = nx.from_pandas_edgelist(redundant_blast, 
                                    source='query_acc', 
                                    target='subject_acc', 
                                    edge_attr='bscore_index')


# \
# for each **connected component** we will only select the most connected node as the components representative,\
# ... and flag everything else for future removal

# In[11]:


to_delete = set()
for cc in nx.connected_components(tmp_graph):
    component_subgraph = tmp_graph.subgraph(cc)
    component_degree   = nx.degree_centrality(component_subgraph)
    #
    # sort connected components based on their degree
    sorted_nodes       = sorted(component_degree.items(),
                                key=lambda item: item[1], 
                                reverse=True)
    #
    # the first one should be the highest degree
    cc.remove(sorted_nodes[0][0])
    #
    # mark everyone else to say goodbye...
    to_delete.update(cc)


# \
# delete everything you flag'd to remove

# In[12]:


non_redundant_blast = blast_search.query('~query_acc.isin(@to_delete) and                                           ~subject_acc.isin(@to_delete)').copy()


# \
# and create the new graph, one **without redundant** connections

# In[13]:


graph = nx.from_pandas_edgelist(non_redundant_blast, 
                                source='query_acc', 
                                target='subject_acc', 
                                edge_attr='bscore_index')


# \
# proteins with hits against three or fewer other proteins are likely to only insert noise here,\
# ... we want to sample sequences representative of whole phyla not all the cases

# In[14]:


to_delete = []
for node in graph.nodes:
    if len(list(graph.neighbors(node)))<=3:
        to_delete.append(node)
graph.remove_nodes_from(to_delete)


# ## Add taxonomic information to nodes
# 
# using **blastdbcmd** we can query an indexed BLAST database for its sequences taxids
# 
# 1. let's create an input file for **blastdbcmd** with a list of all protein accessions

# In[15]:


with open('%s.accessions' % group, 'w') as out:
    out.write('\n'.join(graph.nodes))


# 
# 2. run **blastdbcmd**, retrieving pairs of accessions and taxids

# In[16]:


subprocess.call(['/home/thiberio/.conda/envs/py37/bin/blastdbcmd',
                 '-target_only',
                 '-db', 'recruited_fastas/%s' %group,
                 '-entry_batch', '%s.accessions' % group,
                 '-outfmt', '"%a %T"', # here we specify what kind of information we want from the db
                 '-out', '%s.taxIDs' % group])


# 3. remove quotes from the **blastdbcmd** output to simplify the pandas loading
#  1. delete intermediate files by changing the name of the new one (without quotes) to the old one (with quotes)

# In[17]:


#
# remove quotes from the blastdbcmd output to simplify pandas loading
#
with open('%s.taxIDs-no_quotes' % group, 'w') as stdout:
    subprocess.call(['sed', 's/"//g', '%s.taxIDs' % group],
                    stdout=stdout)
subprocess.call(['mv', 
                 '%s.taxIDs-no_quotes' % group, 
                 '%s.taxIDs' % group])


# 4. load the resulting file to a pandas DataFrame

# In[18]:


taxid_df = pd.read_csv('%s.taxIDs' % group, 
                       sep      =' ', 
                       index_col=0, 
                       header   =None, 
                       names    =['accession', 'taxid'])
taxids = taxid_df.taxid


# 5. obtaining taxids for each accession is just the beggining
#  1. using the taxids we will use ete3's **NCBITaxa** to obtain full taxonomic lineages for each taxid

# In[19]:


taxonomy_df   = pd.DataFrame()
missing_taxid = [] # if we can't obtain its lineage, mark it for deletion
for taxid in taxid_df.taxid.unique():
    try:
        #
        # NCBITaxa makes it very easy, damned are the days of parsing ncbi's taxonomy myself!
        lineage = {j:i
                   for i, j in ncbi.get_rank(
                       ncbi.get_lineage(taxid)
                   ).items()
                  }
    except ValueError:
        missing_taxid.append(taxid)
    else:
        lineage['taxid'] = taxid
        taxonomy_df = taxonomy_df.append(lineage, ignore_index=True)


# \
#  B. this will retrieve lots of ill-defined taxonomic ranks which will be un the name of **Unnamed**, let's keep just the classic ones...

# In[20]:


#
# create a DataFrame for posterior querying
#
columns_to_drop = []
for column in taxonomy_df.columns:
    if column not in ['class', 'species', 'superkingdom', 'genus',
                      'order', 'phylum',  'family',       'kingdom',
                      'taxid']:
        columns_to_drop.append(column)

taxonomy_df.drop(columns_to_drop, axis='columns', inplace=True)
taxonomy_df.set_index('taxid', inplace=True)


#  \
#  C. let's keep only bacterial lineages!
#  > Bacteria's taxID is 2!

# In[21]:


taxonomy_df = taxonomy_df.query('superkingdom==2')


# \
#  D. let's keep going and remove all lineages missing phylum information

# In[22]:


taxonomy_df.dropna(subset=['phylum'], inplace=True)


# \
#  E. now we have to trimm our taxid DataFrame and graph to match the changes we prunning we performed based on the taxonomic information!

# In[23]:


taxids = taxids[taxids.isin(taxonomy_df.index)]

graph.remove_nodes_from(set(graph.nodes).difference(taxids.index))


# ## Select sequence representatives for each phylum
# 
# using the taxonomic information we just retrieved we will create sub-graphs composed exclusively by sequences from the same phylum
# 
# **We are taking a few short-cuts here**
# 1. if a phylum contain 10 or less sequences we will ignore it
#  - it is neither a major phylum nor conserved in one
# 2. only going to sample 10% of sequences with the best **betweenness** from each phylum, but no more than 30 per phylum
# 3. if the graph is too big (i.e. more than 900 nodes) we will limit the calculation of the **betweenness centrality** to 100 random nodes per node

# In[24]:


def  select_phylum_representatives(phylum):

    #
    # select protein accessions from a single phylum, and create a phylum-specific sub-graph
    phylum_taxids     = taxonomy_df.query('phylum==@phylum').index
    phylum_accessions = taxids[taxids.isin(phylum_taxids)].index
    phylum_graph      = graph.subgraph(phylum_accessions)
    
    #
    # if there are less than 10 protein from this phylum, ignore it...
    if phylum_graph.order() < 10:
        return({phylum:[]})
        
    #
    # we will sample 10% of sequences each each phylum, but no more than 30 per phylum
    num_representatives = int(phylum_graph.order()/10) if phylum_graph.order() <= 300 else 30
    
    if phylum_graph.order() > 900:
    #
    # this one is tricky, to speeup betweenness calculations if the graph has more than 900 nodes
    #     we will limit betweenness calculation to random subsampling of 100 other nodes per node...
        phylum_betweenness = nx.betweenness_centrality(phylum_graph, 
                                                       k=100, 
                                                       weight='bscore_index')
    else:
        #
        # if the graph is smaller, just use everything
        phylum_betweenness = nx.betweenness_centrality(phylum_graph, 
                                                       weight='bscore_index')
    
    #
    # sort nodes based on their betweenness centrality and return the 10% top nodes
    #
    sorted_nodes       = sorted(phylum_betweenness.items(),
                                key=lambda item: item[1],
                                reverse=True)
    return({phylum:[item[0] for item in sorted_nodes[:num_representatives]]})


# \
# the most time consuming step here is calculating the **betweenness centrality**, feel free to change the centrality measure to degree. I like using **betweenness** because it better sample different areas of the graph!
# 
# the cluster I am currently running it in has 20 cores per node, adjust has you see fit

# In[25]:


#
# run it in parallel...
#
pool = multiprocessing.Pool(processes=20)
phyum_representatives = pool.map(select_phylum_representatives,
                                 taxonomy_df.phylum.unique())


# \
# here we are only going through the results from our paralellized function, removing mentions from bypassed phyla and throughing them into `dict` structures

# In[26]:


accession_phylum    = {}
all_representatives = {}
for tmp in phyum_representatives:
    #
    # if the no representatives are returned it is because there were too few protein from the phylum
    if tmp.values():
        all_representatives.update(tmp)
        #
        # translate phylum taxids to scientific names to ease visualization
        #
        phylum = list(tmp.keys())[0]
        for acc in list(tmp.values())[0]:
            #
            # here, this is the translating function
            accession_phylum[acc] = ncbi.translate_to_names([phylum])[0]


# \
# we are using **blastdbcmd** once more, this time to retrieve AA sequences for each selected representative
# 
# 1. first we create the input file, same as before
# 2. now we query the database, with one different from before: adding `%s` in the `-outfmt` parameter

# In[27]:


with open('%s.representative.accessions' % group, 'w') as out:
    out.write('\n'.join(accession_phylum.keys()))

subprocess.call(['/home/thiberio/.conda/envs/py37/bin/blastdbcmd',
                 '-target_only',
                 '-db', 'recruited_fastas/%s' % group,
                 '-entry_batch', '%s.representative.accessions' % group,
                 '-outfmt', '"%a %T %s"', # here me added the "%s" since we now want sequences
                 '-out', '%s.representatives.taxID_and_seqs' % group])


# \
# Here we will read the **blastdbcmd** output and write a **FASTA** file with all selected representatives

# In[28]:


with open('%s.representative.fasta' % group, 'w') as out,      open('%s.representatives.taxID_and_seqs' % group) as data:
    for line in data.readlines():
        line = line.strip('"').split()
        out.write('>%s|%s\n%s\n' % (line[0], taxonomy_df.loc[int(line[1]), 'phylum'], line[2]))


# 
# ### Done...
# ... almost!
# 
# the code above is the *important* part as it fulfil our objective, sampling phylum representatives for a gene family from a **BLAST** result.
# 
# ### But now we generate a visualization for our results
# 
# 1. create sub-graph containing only selected representatives
# 2. calculate betweenness centrality for our new sub-graph
# 3. compute 2D positions for our sub-graph
# 4. add 2D position and the estimated betweenness to `Plotly` traces
# 5. feed traces to a `Plotly.Figure` object
# 6. curate your results
# 7. enjoy!

# In[29]:


representative_graph = graph.subgraph(accession_phylum.keys())
#
# calcule some centrality measure, I chose betweenness (my favorite)
graph_betweenness = nx.betweenness_centrality(representative_graph, weight='bscore_index')
#
# calcula 2d positions for each node
node_positions = nx.spring_layout(representative_graph, weight='bscore_index')


# In[30]:


#
# let's generate an interactive visualiztion of our sampling
#     it makes easier to communicate our results and make them useful!
#
edge_x = []
edge_y = []
for edge in representative_graph.edges():
    x0, y0 = node_positions[edge[0]]
    x1, y1 = node_positions[edge[1]]
    edge_x.append(x0)
    edge_x.append(x1)
    edge_x.append(None)
    edge_y.append(y0)
    edge_y.append(y1)
    edge_y.append(None)

edge_trace = go.Scatter(
    x        =edge_x, 
    y        =edge_y,
    opacity  =0.5,
    name     ='Identities',
    line     =dict(width=0.5, 
                   color=('#888')),
    hoverinfo='none',
    mode     ='lines')

node_x      = []
node_y      = []
node_text   = []
node_color  = []
for node in representative_graph.nodes():
    x, y = node_positions[node]
    node_x.append(x)
    node_y.append(y)
    #
    # each node's name will be composed by its accession number and phylum name
    node_text.append('%s|%s' % (node, accession_phylum[node]))
    #
    # colors will come from the betweenness centrality
    node_color.append(graph_betweenness[node])
    
node_trace = go.Scatter(
    x        =node_x, 
    y        =node_y,
    text     =node_text,
    mode     ='markers',
    name     ='Proteins',
    hoverinfo='text',
    opacity  =0.7,
    marker   =dict(showscale   =True,
                   colorscale  ='RdBu',
                   reversescale=True,
                   symbol      ='circle',
                   color       =node_color,
                   size        =10,
                   line_width  =1,
                   line_color  ='white')
)


# In[31]:


fig = go.Figure(data=[edge_trace, node_trace], # here is important to add edges first so they don't cover nodes
                layout=go.Layout(template          ='simple_white',
                                 title             ='protein similarity network',
                                 titlefont_size    =16,
                                 showlegend        =True,
                                 legend_orientation='h',
                                 hovermode         ='closest',
                                 margin            =dict(b=20,
                                                         l=5,
                                                         r=5,
                                                         t=40),
                                 xaxis             =dict(showgrid     =False,
                                                         zeroline     =False,
                                                         showticklabels=False),
                                 yaxis             =dict(showgrid      =False, 
                                                         zeroline      =False, 
                                                         showticklabels=False))
                )

plotly.offline.iplot(fig, 'Identity network', config={'scrollZoom': True})

