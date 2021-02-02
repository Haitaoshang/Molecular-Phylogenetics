#!/usr/bin/env python
# coding: utf-8

# In[1]:


import multiprocessing
import pickle as pkl
import plotly
import plotly.plotly as ptl
from plotly import graph_objs as go
import pyparsing as pp
import subprocess
import os
import jdc

plotly_accession = open('/Users/thiberio/plotly_accession').read().split()
ptl.sign_in(plotly_accession[0], plotly_accession[1])

get_ipython().run_line_magic('run', 'base_functions.ipynb')
os.chdir('/work/Alphas_and_Cyanos')


# In[2]:


#os.chdir('/work/jupyter_notebooks/index_hgt')
#%run base_functions.ipynb
#os.chdir('/work/Alphas_and_Cyanos')


# In[3]:


parse_transfers = aggregate(
    reference_tree='rooted_partitions-with_named_branches.treefile',
    gene_tree_folder='/work/Alphas_and_Cyanos/ranger_input_trees-no_long_branches/',
    aggregate_folder='/work/Alphas_and_Cyanos/aggregated/mad_roots-stricter_branch_lengths/',
    reconciliation_folder='.')


# In[4]:


with cd('reconciliations/mad_roots-stricter_branch_lengths'):
    pool    = multiprocessing.Pool(processes=15)
    results = pool.map(parse_transfers.parse_aggregated, os.listdir('.'))
    pool.close()
    pool.join()

    transfers = {}
    for filtered in results:
        if  list(filtered.values()) != [None] and list(filtered.values())[0][0] != []:
            transfers.update(filtered)

out = open('aggregated/mad_transfers-test.pkl', 'wb')
pkl.dump(transfers, out)
out.close()


# In[5]:


get_ipython().run_line_magic('pinfo', 'parse_transfers.species_tree.write')


# In[6]:


out = open('aggregated/maxtic.constrains-test', 'w')
for group, (transfer_data, gene_tree) in transfers.items():
    for transfer in transfer_data:
        out.write('%s\t%s\n' % (transfer['donor'], transfer['recipient']))
out.close()

parse_transfers.species_tree.write(outfile='tmp_species.tre', format=1, dist_formatter='%.20f')

subprocess.call(['python',
                 '/work/ale/maxtic/MaxTiC.py',
                 'tmp_species.tre',
                 'aggregated/maxtic.constrains-test',
                 'ls=180',
                ])

maxtic = pd.read_table('aggregated/maxtic.constrains-test_MT_output_partial_order',
                       header=None,
                       names=['donor', 'recipient', 'weight', 'no_idea'],
                       sep=' ')


# In[7]:


maxtic_compatible_transfers = {}
for group, (transfer_data, gene_tree) in transfers.items():
    tmp_transfers = []
    for transfer in transfer_data:
        if maxtic[(maxtic.donor==transfer['donor']) & (maxtic.recipient==transfer['recipient'])].shape[0]:
            tmp_transfers.append(transfer)
    if tmp_transfers:
        maxtic_compatible_transfers[group] = [tmp_transfers, gene_tree]


# In[8]:


pool = multiprocessing.Pool(processes=18)
results = pool.map(parse_transfers.assess_dtl_dist,
                   list(maxtic_compatible_transfers.items()))


# In[9]:


donor_dtl_distances = {}
for element in results:
    donor_dtl_distances.update(element)


# In[10]:


reference_tree         = ete3.Tree('rooted_partitions-with_named_branches.treefile', format=1)
transfer_distances     = {}
donor_distance_to_root = {}
donor_complexity_ratio = {}
for group, (transfer_data, gene_tree) in maxtic_compatible_transfers.items():
    for transfer in transfer_data:
        pair         = frozenset([transfer['donor'], transfer['recipient']])
        donor_branch = reference_tree.search_nodes(name=transfer['donor']    )[0]
        
        if pair not in transfer_distances:
            recipient_branch = reference_tree.search_nodes(name=transfer['recipient'])[0]
            transfer_distances[pair] = donor_branch.get_distance(recipient_branch, topology_only=False)

        if transfer['donor'] not in donor_distance_to_root:
            tmp_dist = reference_tree.get_distance(
                transfer['donor'],
                topology_only=False
            )
            donor_distance_to_root[transfer['donor']] = tmp_dist
            donor_subtree_complexity = sum([node.dist
                                            for node in donor_branch.traverse()
                                            if node.name != transfer['donor']])
            tmp_dist = reference_tree.get_distance(
                transfer['donor'],
                topology_only=True
            )            
            donor_complexity_ratio[transfer['donor']] = tmp_dist/len(donor_branch)


# In[11]:


maxtic.loc[
            (maxtic.donor==maxtic_compatible_transfers[group][0][position]['donor']) &
            (maxtic.recipient==maxtic_compatible_transfers[group][0][position]['recipient']),
].shape[0]


# In[ ]:


maxtic.head()


# In[ ]:


tracer = {'color':[], 'x':[], 'y':[], 'text':[], 'marker_size':[]}
for group in donor_dtl_distances.keys():
    for position in range(len(donor_dtl_distances[group])):
        if not maxtic.loc[
            (maxtic.donor==maxtic_compatible_transfers[group][0][position]['donor']) &
            (maxtic.recipient==maxtic_compatible_transfers[group][0][position]['recipient']),
            ].shape[0]:
            continue

        tracer['x'    ].append(
            transfer_distances[frozenset(
                [maxtic_compatible_transfers[group][0][position]['donor'],
                 maxtic_compatible_transfers[group][0][position]['recipient']]
            )]
        )
        tracer['y'    ].append(
            donor_complexity_ratio[maxtic_compatible_transfers[group][0][position]['donor']]
        )
        tracer['text' ].append('%s-#%i' %(group, position))
        tracer['color'].append(donor_dtl_distances[group][position])
        
        transfer_count = maxtic.loc[
            (maxtic.donor==maxtic_compatible_transfers[group][0][position]['donor']) &
            (maxtic.recipient==maxtic_compatible_transfers[group][0][position]['recipient']),
            'weight'].squeeze()
        tracer['marker_size'].append(10+transfer_count*0.7)


color_range          = np.linspace(np.min(tracer['color']), np.max(tracer['color']), 100)
tracer['color_bins'] = np.digitize(tracer['color'], color_range)
tracer_df = pd.DataFrame.from_dict(tracer)

binned_df = tracer_df.groupby(by='color_bins')

bins        = []
for bin in binned_df.groups.keys():
    tmp_df = binned_df.get_group(bin)
    bins.append(
        go.Scatter(
            x=tmp_df.x.values,
            y=tmp_df.y.values,
            mode='markers',
            text=tmp_df.text.values,
            name=str(round(color_range[bin-1], 4)),
            hoverinfo='text',
            showlegend=False,
            marker=dict(
                size=tmp_df.marker_size.values,
                color=tmp_df.color.values,
                colorscale='RdBu',
                cmax=tracer_df.color.values.max(),
                cmin=tracer_df.color.values.min(),
                symbol='circle',
                opacity=.7,
            )
        )
    )

#
# source: https://plot.ly/python/sliders/
steps = [dict(label='All',
                method='restyle',
                args=[
                    'visible', [True] * (len(bins) + 1)
                ])
]
for i in range(len(bins)):
    step = dict(label=bins[i]['name'],
                method='restyle',
                args=[
#                    'visible', [False] * i + [True] * (len(bins) - i)
                    'visible', [False]  * (len(bins))
                ])
    step['args'][1].append(True)
    step['args'][1][i] = True
    steps.append(step)
slider = dict(steps=steps, currentvalue={'prefix':'Donor subtree DTL: '}, pad={'t':50})
bins.append(
    go.Scatter(
        x=[np.min(tracer['x']), np.max(tracer['x'])],
        y=[np.min(tracer['y']), np.max(tracer['y'])],
        showlegend=False,
        mode='markers',
        marker=dict(
            size=10, 
            color=[0.5], 
            colorscale='RdBu', 
            cmax=np.max(tracer['color']), 
            cmin=np.min(tracer['color']), 
            symbol='circle', 
            opacity=0,
            colorbar=dict(title='Donor subtree DTL cost')
        )
    )
)

layout    = go.Layout(
    title='Donor/Recipient subtree reconciliation costs',
    hovermode='closest',
    width=1200, height=1000,
    xaxis=dict(title='Donor-Recipient distance'),
    yaxis=dict(title='Donor branch distance to root and donor subtree ratio'),
    sliders=[slider])
fig       = go.Figure(data=bins, layout=layout)
plot      = plotly.offline.plot(fig, filename='./yeah.html', auto_open=False)

