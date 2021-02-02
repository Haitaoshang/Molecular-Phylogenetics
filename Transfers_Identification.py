#!/usr/bin/env python
# coding: utf-8

# In[2]:


from matplotlib import pyplot as plt
import seaborn as sns
import subprocess
import itertools
from copy import deepcopy

plotly_accession = open('/Users/plotly_accession').read().split()
ptl.sign_in(plotly_accession[0], plotly_accession[1])


# In[3]:


if not os.path.isdir('reconciliation_aggregates'):
    os.mkdir('reconciliation_aggregates')
for folder in os.listdir('ranger/'):
    if os.path.isfile('reconciliation_aggregates/%s' % folder):
        continue
    os.system('cp ranger/%s/aggregated reconciliation_aggregates/%s' % (folder, folder))


# In[4]:


if not os.path.isdir('renamed_gene_tree'):
    os.mkdir('renamed_gene_tree')

for treefile in os.listdir('gene_trees/'):
    if not treefile.endswith('.treefile'):
        continue
    if os.path.isfile('renamed_gene_tree/%s' % treefile):
        continue
    
    tmp_tree = ete3.Tree('gene_trees/%s' % treefile)
    for leaf in tmp_tree.get_leaves():
        if leaf.name.count('_') == 1:
            gene, genome = leaf.name.split('_')
        elif leaf.name.count('_') > 1 and re.search('GC[AF]_', leaf.name):
            gene, genome = re.search('^([^.]+).+?(GC[AF]_\d+)', leaf.name, re.M).groups()
        elif leaf.name.count('_') == 2 and re.search('_PRJ', leaf.name):
            gene, genome = re.search('^(.+)_(PRJ.+)$', leaf.name, re.M).groups()
        else:
            print(leaf.name)
        gene = gene.split('.')[0]
        leaf.name = '%s_%s' % (genome.replace('_', ''), gene.replace('_', ''))
    
    tmp_tree.write(outfile='renamed_gene_tree/%s' % treefile, format=0, dist_formatter='%.10f')


# In[5]:


ranger_parser = aggregate(reference_tree='species_tree-renamed',
                          gene_tree_folder='renamed_gene_tree/',
                          aggregate_folder='reconciliation_aggregates/',
                          reconciliation_folder='ranger/',
                          overall_tree_support_thresh=20,
                          leaves_allowed=False)


# In[73]:


if not os.path.isfile('transfers.tab'):
    transfers = []
    for count, group in enumerate(os.listdir('ranger')):
        if not os.path.getsize('ranger/%s/%s.output1' % (group, group)):
            continue
        transfers.append(ranger_parser.parse_aggregated(group))

    transfer_df = pd.concat([n[0] for n in transfers 
                             if type(n) is not dict], ignore_index=True, axis=0, sort=False)
    transfer_df.loc[:,['ranger_confidence',
                       'ranger_confidence_donor',
                       'ranger_confidence_recipient']] = transfer_df.loc[:,['ranger_confidence',
                                                                            'ranger_confidence_donor',
                                                                            'ranger_confidence_recipient']] / 50

    transfer_df.to_csv('transfers.tab', sep='\t')
else:
    transfer_df = pd.read_csv('transfers.tab', sep='\t', index_col=0)

print(transfer_df.shape)
transfer_df.head()


# In[7]:


fig, axs = plt.subplots(nrows=4, figsize=(15,10), dpi=300)
colors   = sns.color_palette(n_colors=4).as_hex()
for row, column in enumerate(['bipartition_support', 
                              'ranger_confidence',
                              'ranger_confidence_donor',
                              'ranger_confidence_recipient']):
    sns.kdeplot(transfer_df[column],
                shade=True, color=colors[row], ax=axs[row])


# In[77]:


transfer_df = transfer_df[(transfer_df.bipartition_support        >=95)]
print(transfer_df.shape, len(transfer_df.groupby('donor recipient'.split()).groups))
transfer_df.head()


# In[78]:


transfer_df = transfer_df[(transfer_df.ranger_confidence_donor    >=0.8)]
print(transfer_df.shape, len(transfer_df.groupby('donor recipient'.split()).groups))
transfer_df.head()


# In[79]:


transfer_df = transfer_df[(transfer_df.ranger_confidence_recipient    >=0.8)]
print(transfer_df.shape, len(transfer_df.groupby('donor recipient'.split()).groups))
transfer_df.head()


# In[9]:


ranger_parser.name_species_tree_nodes(
    reconciliation_file='ranger/4468_BAB72593.1-GCA_000009705.1/4468_BAB72593.1-GCA_000009705.1.output1'
)


# In[ ]:


ranger_parser.species_tree.write(outfile='species_tree_named_nodes',
                               format=1,
                               dist_formatter='%.10f')


# In[ ]:


out = open('maxtic.input', 'w')
for index, row in transfer_df[['donor', 'recipient']].iterrows():
    out.write('%s\n' % '\t'.join(row.tolist()))
out.close()


# In[ ]:


subprocess.call(['python2.7',
                 '/work/ale/maxtic/MaxTiC.py',
                 'species_tree_named_nodes',
                 'maxtic.input',
                 'ls=180'])


# In[10]:


donor_recipient_pairs = []
for line in open('maxtic.input_MT_output_partial_order').readlines():
    donor_recipient_pairs.append('-'.join(line.split()[:2]))

transfer_df['donor-recipient'] = transfer_df['donor']+'-'+transfer_df['recipient']
transfer_df = transfer_df[transfer_df['donor-recipient'].isin(donor_recipient_pairs)]
transfer_df.drop(labels='donor-recipient', axis=1, inplace=True)
transfer_df.info()


# In[11]:


extended_df = ranger_parser.assess_transfer_distance(transfer_df)
extended_df = ranger_parser.assess_dtl_cost(extended_df)

print(extended_df.shape)
extended_df.head()


# In[12]:


backup = extended_df.copy()


# In[55]:


extended_df = backup.copy()

extended_df['donor_ancestry'] = extended_df.apply(
    lambda row:\
    [ancestor.name for ancestor in next(
        ranger_parser.species_tree.iter_search_nodes(name=row.donor)).get_ancestors()],
    axis=1)
extended_df['recipient_ancestry'] = extended_df.apply(
    lambda row:\
    [ancestor.name for ancestor in next(
        ranger_parser.species_tree.iter_search_nodes(name=row.recipient)).get_ancestors()],
    axis=1)
extended_df.head()


# In[56]:


clusters = []
for name, row in extended_df.iterrows():
    matching_transfers = extended_df[(
                (extended_df.donor_ancestry.apply(lambda x: row.donor in x)) |
                (extended_df.donor == row.donor)
            ) &
            (
                (extended_df.recipient_ancestry.apply(lambda x: row.recipient in x)) |
                (extended_df.recipient == row.recipient)
            )].index
    existing_cluster = False
    for index, cluster in enumerate(clusters):
        if not cluster.isdisjoint(matching_transfers):
            existing_cluster = True
            clusters[index].update(matching_transfers)
            break
    if not existing_cluster:
        clusters.append(set(matching_transfers))


# In[153]:


fig, ax             = plt.subplots(figsize=(15, 5), dpi=300)
cluster_count       = 0
donor_distances     = []
recipient_distances = []
for cluster in clusters:
    if len(cluster) < 4:
        continue
    for pair in itertools.combinations(cluster, 2):
        donor1 = next(ranger_parser.species_tree.iter_search_nodes(name=extended_df.loc[pair[0], 'donor']))
        donor2 = next(ranger_parser.species_tree.iter_search_nodes(name=extended_df.loc[pair[1], 'donor']))
        donor_distances.append(
            ranger_parser.species_tree.get_distance(donor1, donor2, topology_only=True)
        )

        recipient1 = next(ranger_parser.species_tree.iter_search_nodes(name=extended_df.loc[pair[0], 'recipient']))
        recipient2 = next(ranger_parser.species_tree.iter_search_nodes(name=extended_df.loc[pair[1], 'recipient']))
        recipient_distances.append(
            ranger_parser.species_tree.get_distance(recipient1, recipient2, topology_only=True)
        )
    
sns.distplot(donor_distances, kde=False, ax=ax, label='Nodes between donors')
sns.distplot(recipient_distances, kde=False, ax=ax, label='Nodes between recipients')
ax.legend()


# In[154]:


fig, ax             = plt.subplots(figsize=(15, 5), dpi=300)
cluster_count       = 0
donor_distances     = []
recipient_distances = []
for cluster in clusters:
    if len(cluster) < 4:
        continue
    for pair in itertools.combinations(cluster, 2):
        donor1 = next(ranger_parser.species_tree.iter_search_nodes(name=extended_df.loc[pair[0], 'donor']))
        donor2 = next(ranger_parser.species_tree.iter_search_nodes(name=extended_df.loc[pair[1], 'donor']))
        donor_distances.append(
            ranger_parser.species_tree.get_distance(donor1, donor2, topology_only=False)
        )

        recipient1 = next(ranger_parser.species_tree.iter_search_nodes(name=extended_df.loc[pair[0], 'recipient']))
        recipient2 = next(ranger_parser.species_tree.iter_search_nodes(name=extended_df.loc[pair[1], 'recipient']))
        recipient_distances.append(
            ranger_parser.species_tree.get_distance(recipient1, recipient2, topology_only=False)
        )
    
sns.distplot(donor_distances, kde=False, ax=ax, label='Branch length between donors')
sns.distplot(recipient_distances, kde=False, ax=ax, label='Branch length between recipients')
ax.legend()


# In[283]:


count = 80
print(extended_df.loc[clusters[count], 'donor'].unique())
print()
print(extended_df.loc[clusters[count], 'recipient'].unique())


# In[284]:


clusters[80]


# In[282]:


#define constraints for each cluster
def myLayout(node):
    node.img_style["vt_line_color"] = "#000000"
    node.img_style["hz_line_color"] = "#000000"
    node.img_style["vt_line_width"] = 1
    node.img_style["hz_line_width"] = 1
    
    if node.is_leaf():
        node.img_style['size'] = 0
    else:
        if node.ranger_name in possible_donors and node.ranger_name in possible_recipients:
            node.img_style['fgcolor'] = 'purple'
        elif node.ranger_name in possible_donors:
            node.img_style['fgcolor'] = 'red'
        elif node.ranger_name in possible_recipients:
            node.img_style['fgcolor'] = 'green'
        else:
            node.img_style['size'] = 0

treeStyle                    = ete3.TreeStyle()
treeStyle.layout_fn          = myLayout
treeStyle.show_leaf_name     = False
treeStyle.branch_vertical_margin = 1.5

for count, cluster in enumerate(clusters):
    if len(cluster) < 4:
        continue
    tmp_tree = ranger_parser.species_tree.copy()
    
    possible_donors       = set(extended_df.loc[cluster, 'donor'    ].unique().tolist())
    donor_constrainer     = possible_donors.pop()
    constrainer_depth     = extended_df.loc[extended_df.donor==donor_constrainer,
                                            'donor_depth'].tolist()[0]
    for donor in possible_donors:
        tmp_depth = extended_df.loc[extended_df.donor==donor, 'donor_depth'].tolist()[0]
        if tmp_depth > constrainer_depth:
            constrainer_depth = tmp_depth
            donor_constrainer = donor

    possible_recipients   = set(extended_df.loc[cluster, 'recipient'].unique().tolist())
    recipient_constrainer = possible_recipients.pop()
    constrainer_depth     = extended_df.loc[extended_df.recipient==recipient_constrainer,
                                            'recipient_depth'].tolist()[0]
    for recipient in possible_recipients:
        tmp_depth = extended_df.loc[extended_df.recipient==recipient, 'recipient_depth'].tolist()[0]
        if tmp_depth < constrainer_depth:
            constrainer_depth = tmp_depth
            recipient_constrainer = recipient
    
    possible_donors       = set(extended_df.loc[cluster, 'donor'    ].unique().tolist())
    possible_recipients   = set(extended_df.loc[cluster, 'recipient'].unique().tolist())

    treeStyle                    = ete3.TreeStyle()
    treeStyle.layout_fn          = myLayout
    treeStyle.show_leaf_name     = False
    treeStyle.branch_vertical_margin = 1.5
    treeStyle.title.add_face(ete3.TextFace("cluster: %i" % count, fsize=10), column=0)
    tree_plot = tmp_tree.render('hgt_clusters/cluster_%i.png' % count,
                                tree_style=treeStyle, dpi=1200, w=2000, units='px')
    
    #print(count, donor_constrainer, recipient_constrainer)


# In[13]:


fig, axs = plt.subplots(nrows=5, figsize=(18,10), dpi=300)
colors   = sns.color_palette(n_colors=5).as_hex()
for row, column in enumerate(['donor_depth', 
                              'recipient_depth',
                              'donor_subtree_size',
                              'recipient_subtree_size',
                              'donor_dtl_size_ratio']):
    sns.kdeplot(extended_df.loc[(extended_df.donor_depth<=5)&(extended_df.recipient_depth<=5),
                                column],
                shade=True, color=colors[row], ax=axs[row])


# In[15]:


extended_df = ranger_parser.map_taxonomic_level(extended_df, taxa_table='../genomes.tab')
extended_df.head()


# In[16]:


phylum_transfers = extended_df.loc[extended_df.transfer_level == 'phylum']
kingdom_transfers = extended_df.loc[(extended_df.transfer_level == 'kingdom') |
                                    (extended_df.transfer_level == 'superkingdom')]


# In[17]:


gregs_table      = pd.read_excel('TimeCalibratingHGTList_7_10.xlsx')
useful_transfers = gregs_table.loc[gregs_table['Use?']=='Y', 'd/r pair'].tolist()
useful_transfers = extended_df.reindex(index=useful_transfers, copy=True)
useful_transfers.dropna(axis=0, how='all', inplace=True)


# In[69]:


get_ipython().run_cell_magic('add_to', 'aggregate', 'def visualize_in_gene_figtree(self, df, taxa_table=None):\n    ncbi     = ete3.NCBITaxa()\n\n    taxa_df = pd.read_csv(taxa_table, sep=\'\\t\')\n    taxa_df[\'Unnamed: 0\'] = taxa_df[\'Unnamed: 0\'].apply(lambda x: x.split(\'.\')[0])\n    taxa_df[\'accession\'] = taxa_df[\'accession\'].apply(lambda x: x.split(\'.\')[0])\n    taxa_df.set_index(\'Unnamed: 0\', inplace=True)\n    \n    folder = \'highlighted_gene_trees\'\n    if not os.path.isdir(folder):\n        os.mkdir(folder)\n    else:\n        os.system(\'rm -rf %s/*\' % folder)\n    \n    for group in df.family.unique():\n        group_df = df.query(\'family==@group\')\n        group_num   = int(group.split(\'_\')[0])\n        newick_text = open(\'gene_trees/%s.treefile.rooted\' % group).read()\n        gene_tree   = ete3.Tree(re.sub(\'\\[.*\\];$\', \';\', newick_text.strip(), flags=re.M))\n        for leaf in gene_tree.get_leaves():\n            if leaf.name.count(\'_\') == 1:\n                gene, genome = leaf.name.split(\'_\')\n            elif leaf.name.count(\'_\') > 1 and re.search(\'GC[AF]_\', leaf.name):\n                gene, genome = re.search(\'^([^.]+).+?(GC[AF]_\\d+)\', leaf.name, re.M).groups()\n            elif leaf.name.count(\'_\') == 2 and re.search(\'_PRJ\', leaf.name):\n                gene, genome = re.search(\'^(.+)_(PRJ.+)$\', leaf.name, re.M).groups()\n            else:\n                print(leaf.name)\n            gene = gene.split(\'.\')[0]\n            leaf.add_feature(\'true_name\', leaf.name)\n            leaf.add_feature(\'genome\', genome)\n            leaf.name = \'%s_%s\' % (genome.replace(\'_\', \'\'), gene.replace(\'_\', \'\'))\n        \n        tmp = self.name_branches_as_reconciliation(open(\'ranger/%s/%s.output1\' % (group, group)).read(),\n                                                   gene_tree)\n        out  = open(\'%s/group_%i-hgt.figTree\' % (folder, group_num), \'w\')\n        out.write("#NEXUS\\nbegin taxa;\\n\\tdimensions ntax=%i;\\n\\ttaxlabels\\n" %len(gene_tree))\n        branch_names = {}\n        for count, node in enumerate(tmp[0].traverse()):\n\n            if node.is_leaf():\n                if node.genome in taxa_df.index:\n                    node_name = node.genome\n                elif node.genome in taxa_df.accession.values:\n                    node_name = taxa_df.query(\'accession==@node.genome\').index[0]\n                else:\n                    out.write(\'\\t%s\\n\' %(node.true_name))\n                    continue\n\n                comment = [\'source_name="%s"\' % taxa_df.loc[node_name, \'Organism\']]\n                if pd.isnull(taxa_df.loc[node_name, \'taxid\']):\n                    out.write(\'\\t%s \' %(node.name))\n                else:\n                    taxid = taxa_df.loc[node_name, \'taxid\']\n                    lineage = {j:i\n                               for i, j in ncbi.get_rank(\n                                   ncbi.get_lineage(taxid)).items()\n                              }\n                    lineage_names = ncbi.get_taxid_translator(lineage.values())\n\n                    out.write(\'\\t%s \' % (node.name))\n                    for rank in [\'class\', \'phylum\', \'order\', \'family\']:\n                        if rank in lineage:\n                            comment.append(\'tax_%s="%s"\' % (rank, lineage_names[lineage[rank]]))\n                out.write(\'[&%s]\\n\' %\' \'.join(comment))\n\n            else:\n                branch_names[\'_branch_%i_\' % count] = \'&ranger_name=%s\' %(node.ranger_name)\n                as_donor     = {}\n                as_recipient = {}\n                for index, row in group_df.query(\'donor_map==@node.ranger_name\').iterrows():\n                    if not row.family.split(\'_\')[0] in as_donor:\n                        as_donor[row.family.split(\'_\')[0]] = \'\'\n                    as_donor[row.family.split(\'_\')[0]] += \'#%i\' % index\n                for index, row in group_df.query(\'recipient_map==@node.ranger_name\').iterrows():\n                    if not row.family.split(\'_\')[0] in as_recipient:\n                        as_recipient[row.family.split(\'_\')[0]] = \'\'\n                    as_recipient[row.family.split(\'_\')[0]] += \'#%i\' % index\n\n                for group, role in as_donor.items():\n                    branch_names[\'_branch_%i_\' % count] += \',role=donor%s\' % role\n                for group, role in as_recipient.items():\n                    if group in as_donor:\n                        branch_names[\'_branch_%i_\' % count] += \'/recipient%s\' % role\n                    else:\n                        branch_names[\'_branch_%i_\' % count] += \',role=recipient%s\' % role\n\n                node.name = \'_branch_%i_\' % count\n\n        newick_text = tmp[0].write(format=1, dist_formatter=\'%.10f\')\n        for key, value in branch_names.items():\n            newick_text = newick_text.replace(key, \'[%s]\' % value)\n        out.write(\';\\nend;\\n\')\n        out.write(\'begin trees;\\n\\ttree tree_1 = [&R] %s\\nend;\' %newick_text)\n        out.close()')


# In[71]:


len(transfer_df.family.unique())


# In[70]:


ranger_parser.visualize_in_gene_figtree(transfer_df, '../genomes.tab')


# In[151]:


ranger_parser.visualize_in_gene_figtree(useful_transfers, '../genomes.tab')


# In[ ]:


ranger_parser.interactive_dynamic_plot(
    extended_df)


# In[295]:


useful_transfers[useful_transfers.family.str.startswith('21_')].copy()


# In[305]:


ranger_parser.visualize_in_gene_figtree(useful_transfers[useful_transfers.family.str.startswith('21_')],
                                       '../genomes.tab')


# In[ ]:


get_ipython().run_cell_magic('add_to', 'aggregate', "def interactive_dynamic_plot(self, df):\n    tracers = []\n    max_x   = 0\n    max_y   = 0\n    rank_order = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']\n    rank_color = dict(zip(rank_order, np.linspace(0,1,15)[1::2]))\n    rank_color['superkingdom'] = rank_color['kingdom']\n\n    colors = cl.scales['9']['qual']['Paired']\n    colorscale = []\n    for pos, value in enumerate(np.linspace(0,1,8)):\n        colorscale.append([value, colors[pos]])\n        colorscale.append([value, colors[pos+1]])\n\n    for group in df.family.unique():\n        tracer = {'x':[], 'y':[], 'text':[], 'marker_color':[], 'marker_size':[]}\n        group_number = int(group.split('_')[0])\n        for index, row in df.query('family==@group').iterrows():\n            corroborated_by = extended_df.query('donor==@row.donor & recipient==@row.recipient').shape[0]\n            tracer['x'   ].append(row.donor_recipient_distance)\n            tracer['y'   ].append(row['donor_depth/size_ratio'])\n            tracer['text'].append('group_%i#%i' % (group_number, index))\n            tracer['marker_size'].append(10+corroborated_by*2)\n            tracer['marker_color'].append(rank_color[row.transfer_level]\n                                          if pd.notnull(row.transfer_level) else 1)\n    \n        if np.max(tracer['x']) > max_x:\n            max_x = np.max(tracer['x'])\n        if np.max(tracer['y']) > max_y:\n            max_y = np.max(tracer['y'])\n        tracer = go.Scatter(x=tracer['x'],\n                            y=tracer['y'],\n                            mode='markers',\n                            text=tracer['text'],\n                            name='group_%s' % group.split('_')[0],\n                            hoverinfo='text', showlegend=True,\n                            marker=dict(size=tracer['marker_size'],\n                                        color=tracer['marker_color'],\n                                        colorscale=colorscale,\n                                        cmax=1,\n                                        cmin=0,\n                                        symbol='circle',\n                                        opacity=0.7)\n                           )\n        tracers.append(tracer)\n\n    tracers = sorted(tracers, key = lambda x: int(x['name'].split('_')[1]))\n    \n    layout    = go.Layout(\n        title='Interactive index HGT candidates plot!',\n        hovermode='closest',\n        width=1500, height=1000,\n        xaxis=dict(title='Donor-Recipient distance', \n                   autorange=False, \n                   range=[0, max_x+max_x*0.01]),\n        yaxis=dict(title='Donor depth/size ratio', \n                   autorange=False, \n                   range=[0, max_y+max_y*0.01]),\n        updatemenus=[\n            {'buttons':[{'label':'Show all',\n                         'method':'restyle',\n                         'args': [ 'visible', True]},\n                        {'label':'Hide all',\n                         'method':'restyle',\n                         'args': [ 'visible', ['legendonly']*len(tracers)+[True]]}]}\n        ]\n    )\n    \n    tracers.append(go.Scatter(x=[max_x],\n                             y=[max_y],\n                             mode='markers',\n                             name='colorbar',\n                             showlegend=False,\n                             marker=dict(size=10,\n                                        color=0,\n                                        symbol='circle',\n                                        opacity=0.0,\n                                        colorscale=colorscale,\n                                        cmin=0,\n                                        cmax=1,\n                                        colorbar=dict(title='HGT within:',\n                                                      x=1.25,\n                                                      titleside = 'top',\n                                                      tickvals = np.linspace(0,1,15)[1::2],\n                                                      ticktext = rank_order,\n                                                      ticks = 'outside')\n                                        )\n                           )\n                  )\n    \n    fig       = go.Figure(data=tracers, layout=layout)\n    plot      = plotly.offline.plot(fig, filename='./test.html', auto_open=False)")


# In[ ]:


#ranger_parser.name_species_tree_nodes(
#    reconciliation_file='ranger/4468_BAB72593.1-GCA_000009705.1/4468_BAB72593.1-GCA_000009705.1.output1'
#)
ranger_parser.visualize_in_figtree(kingdom_transfers, taxa_table='../genomes.tab')

