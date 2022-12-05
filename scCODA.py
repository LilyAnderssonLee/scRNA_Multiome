#!/usr/bin/env python
# coding: utf-8

# **Current analysis was run under conda env: /proj/snic2021-23-715/private/Lili/analysis/conda_env/scRNA_MouseDevel_env**
# 
# **scCODA - Compositional analysis of single-cell data**
# 
# **citation: scCODA is a Bayesian model for compositional single-cell data analysis**
# 
# **Github: https://github.com/theislab/scCODA**
# 

# In[2]:


get_ipython().system('pip install sccoda')


# In[3]:


# Setup
import importlib
import warnings
import platform
warnings.filterwarnings("ignore")

import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt

from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization as viz

import sccoda.datasets as scd
import scanpy as sc


# In[4]:


sc.settings.verbosity = 3             
sc.settings.set_figure_params(dpi=150,figsize=(4,4))


# In[21]:


print(platform.platform())


# In[46]:


pip list


# # Developing mouse model

# In[5]:


mouse = sc.read_h5ad('../results/Mouse_harmony_rmBatchbetweenSamples_rmDoublets.h5ad')


# In[6]:


mouse


# In[8]:


sc.pl.umap(mouse,color=['Age','CellType_edited4'],wspace=0.3,size=2)


# In[7]:


mouse = mouse[mouse.obs.Age!='Adult']


# In[8]:


mouse


# ## boxplot & T test

# Prepare the cell counts data frame 

# ### Relative abundance

# In[37]:


n_cell_df = mouse.obs.groupby(["batch", "CellType_edited4"]).size().reset_index(name="n_cells")
n_cell_df['Age'] = [iterm.split('_')[4] for iterm in n_cell_df.batch]
n_cell_df.rename(columns={'CellType_edited4':'CellType'},inplace=True)

n_cell_df.head()


# In[38]:


n_cells_batch = mouse.obs.batch.value_counts().sort_index().tolist()
#five cell types within each batch
batch_cells = np.repeat(n_cells_batch,n_cell_df.batch.value_counts().sort_index().tolist())


# In[39]:


n_cell_df['% Cells'] = n_cell_df['n_cells']/batch_cells*100


# In[45]:


import natsort
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator


# In[46]:


fig, ax1 = plt.subplots(nrows=1,ncols=1,constrained_layout=True,figsize=(8,6))
ax = sns.boxplot(data = n_cell_df,x='CellType',y='% Cells',hue='Age',ax=ax1,width=0.5)
sns.stripplot(data = n_cell_df,x='CellType',y='% Cells',hue='Age',ax=ax1,dodge=True)
plt.setp(ax.get_xticklabels(), rotation=45)
plt.legend(bbox_to_anchor=[1.1,1])
plt.ylabel('% Cells')
plt.xlabel('Cell Type')

handles, labels = ax.get_legend_handles_labels()
# When creating the legend, only use the first three elements
l = plt.legend(handles[0:3], labels[0:3], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

#add annotation of T- test between P3 and P60
annot_df = n_cell_df[n_cell_df.CellType!='NeuroGlia']

hue_plot_params = {
    'data': annot_df,
    'x': 'CellType',
    'y': '% Cells',
    "hue": "Age"
}

pairs =[
    [('Epithelial','P3'),('Epithelial','P60')],
    [('Mesenchymal','P3'),('Mesenchymal','P60')],
    [('Endothelial','P3'),('Endothelial','P60')],
    [('Immune','P3'),('Immune','P60')],
]


annotator = Annotator(ax, pairs, **hue_plot_params)
annotator.configure(test="Mann-Whitney").apply_and_annotate()


# ### Normalization

# In[17]:


pd.DataFrame(mouse.obs.value_counts(['batch','CellType_edited4'],normalize=True)).reset_index()


# In[10]:


freq_cell_df = pd.DataFrame(mouse.obs.value_counts(['CellType_edited4','batch'],normalize=True)).reset_index()
freq_cell_df.rename(columns={0:'freq_CellType','CellType_edited4':'CellType'},inplace=True)
freq_cell_df['Age'] = [iterm.split('_')[4] for iterm in freq_cell_df.batch]
# sort the df
freq_cell_df.sort_values(['CellType','Age'],inplace=True,key=natsort.natsort_keygen())
freq_cell_df['-log2(freq_CellType)'] = -(np.log2(freq_cell_df.freq_CellType))
freq_cell_df['pct'] = freq_cell_df.freq_CellType*100
freq_cell_df['n_cells'] = n_cell_df.n_cells
freq_cell_df.head()


# In[47]:


fig, ax1 = plt.subplots(nrows=1,ncols=1,constrained_layout=True,figsize=(8,6))
ax = sns.boxplot(data = freq_cell_df,x='CellType',y='-log2(freq_CellType)',hue='Age',ax=ax1,width=0.5)
sns.stripplot(data = freq_cell_df,x='CellType',y='-log2(freq_CellType)',hue='Age',ax=ax1,dodge=True)
plt.setp(ax.get_xticklabels(), rotation=45)
plt.legend(bbox_to_anchor=[1.1,1])
plt.ylabel('-log2(pct)')
plt.xlabel('Cell Type')

handles, labels = ax.get_legend_handles_labels()
# When creating the legend, only use the first three elements
l = plt.legend(handles[0:3], labels[0:3], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

#add annotation of T- test between P3 and P60
annot_df = freq_cell_df[freq_cell_df.CellType!='NeuroGlia']

hue_plot_params = {
    'data': annot_df,
    'x': 'CellType',
    'y': '-log2(freq_CellType)',
    "hue": "Age"
}

pairs =[
    [('Epithelial','P3'),('Epithelial','P60')],
    [('Mesenchymal','P3'),('Mesenchymal','P60')],
    [('Endothelial','P3'),('Endothelial','P60')],
    [('Immune','P3'),('Immune','P60')],
]


annotator = Annotator(ax, pairs, **hue_plot_params)
annotator.configure(test="Mann-Whitney").apply_and_annotate()


# **Conclusion:**
# 
# None of cell types are significantly changed during the postnatal development. 

# ## scCODA

# In[11]:


cov_df = pd.DataFrame({"Age": ["P3", "P3", "P3",'P11','P11','P60','P60','P60']}, index=mouse.obs.batch.unique().tolist())
print(cov_df)


# In[12]:


data_scanpy_1 = dat.from_scanpy(
    mouse,
    cell_type_identifier="CellType_edited4",
    sample_identifier="batch",
    covariate_df=cov_df
)
print(data_scanpy_1)


# ## Grouped boxplots: relative abundance

# In[18]:


# Grouped boxplots. No facets, relative abundance, no dots.
viz.boxplots(
    data_scanpy_1,
    feature_name="Age",
    plot_facets=False,
    y_scale="relative",
    add_dots=True,
    figsize=(12,6),
    dpi=200,
)
plt.show()


# In[24]:


# Grouped boxplots. Facets, log scale, added dots and custom color palette.
viz.boxplots(
    data_scanpy_1,
    feature_name="Age",
    plot_facets=True,
    y_scale="log",
    add_dots=True,
    cmap="Reds",
    dpi=200
)
plt.show()


# ## Re-run Composition Analysis by different reference levels, automatically chosen by scCODA

# In[25]:


viz.rel_abundance_dispersion_plot(
    data=data_scanpy_1,
    abundant_threshold=0.9,
    figsize=(8,8),dpi=100
)
plt.show()


# In[26]:


model_all1 = mod.CompositionalAnalysis(data_scanpy_1, formula="Age", reference_cell_type="Endothelial")
all_results1 = model_all1.sample_hmc()
all_results1.summary()


# In[35]:


all_results1.set_fdr(est_fdr=0.3)
all_results1.summary()


# In[36]:


model_all2 = mod.CompositionalAnalysis(data_scanpy_1, formula="C(Age,Treatment('P3'))", reference_cell_type="Endothelial")
all_results2 = model_all2.sample_hmc()
all_results2.summary()


# In[39]:


all_results2.set_fdr(est_fdr=0.3)
all_results2.summary()


# **Endothelial is reference and FDR=0.3, Epithelial and Mesenchymal cells credible decreased in P60 Comparing to P3**

# # Zoom in Epithelial cells

# ## Top markers

# In[6]:


epith = sc.read_h5ad('../subcell_analysis/results/Epithelial_har_noDbl_2.h5ad')


# In[7]:


epith = epith[epith.obs.Age!='Adult']


# In[8]:


epith.obs.Epithelial_subcell.replace({'EP Proliferating':'Ciliogenesis+prolif','Ciliogenesis':'Ciliogenesis+prolif'},inplace=True)


# In[110]:


sc.pl.umap(epith,color=['Age','louvain_0.2','Epithelial_subcell'],wspace=0.5,size=3)


# In[13]:


epith.obs['Epithelial_subcell_edited'] = epith.obs['louvain_0.2']
epith.obs['Epithelial_subcell_edited'].replace({'0':'EP1+3','1':'EP2','2':'Ciliogenesis+prolif','3':'Gata3+'},inplace=True)


# In[44]:


sc.pl.umap(epith,color=['Age','Epithelial_subcell_edited'],wspace=0.5,size=3)


# In[45]:


sc.tl.rank_genes_groups(epith,'Epithelial_subcell_edited',method='wilcoxon',key_added='wilcoxon_Epithelial_subcell_edited')


# In[46]:


sc.tl.filter_rank_genes_groups(epith, groupby='Epithelial_subcell_edited', key = "wilcoxon_Epithelial_subcell_edited", 
                               key_added='wilcoxon_filtered_Epithelial_subcell_edited',min_in_group_fraction=0.1,
                               max_out_group_fraction=0.5,min_fold_change=0.25)


# ## Go analysis

# In[48]:


def Enrichr_Plot(geneset_name,enr_res,clt):
    df = enr_res.results[enr_res.results['Gene_set']==geneset_name]
    df_sig = df[df['Adjusted P-value']<0.05]
    x_values = df_sig['Adjusted P-value'].apply(lambda x: -np.log10(x)).to_list()
    y_values = df_sig['Term']
    
    if df_sig.shape[0]>20:
        ax = sns.barplot(x = x_values[:30], y =y_values[:30],color='#78c679')
        ax.set_title('Enrichment from %s were found in %s' %(geneset_name,clt))
        ax.set_xlabel('-log10(Adjusted P-value)')
    else:
        ax = sns.barplot(x = x_values, y =y_values,color='#78c679')
        ax.set_title('Enrichment from %s were found in %s' %(geneset_name,clt))
        ax.set_xlabel('-log10(Adjusted P-value)')


# ### EP2

# In[26]:


import gseapy
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# In[49]:


glist_cilio = sc.get.rank_genes_groups_df(epith, group='EP2', key='wilcoxon_filtered_Epithelial_subcell_edited')[['names','logfoldchanges','pvals_adj','scores']]
glist_cilio = glist_cilio[~glist_cilio['names'].isnull()]
enr_res_cilio = gseapy.enrichr(gene_list=glist_cilio['names'][glist_cilio['pvals_adj']<0.05].tolist(),
                         organism='mouse',
                         gene_sets=['GO_Biological_Process_2021','GO_Molecular_Function_2021'],
                         description='pathway',
                         cutoff = 0.5,
                         outdir ='GoTerm_subcell/Enrichr/Epithelial/EP2')
# check if enrichment were found in all geneset
for geneset_name in ['GO_Biological_Process_2021','GO_Molecular_Function_2021']:
    df = enr_res_cilio.results[enr_res_cilio.results['Gene_set']==geneset_name]
    df_sig = df[df['Adjusted P-value']<0.05]
    
    if df_sig.shape[0]==0:
        print('No enrichments from %s!' %geneset_name)
        
    else:
        print('Oh YEAH there are enrichments in %s!' %geneset_name)
    
enr_res_cilio_sig = enr_res_cilio.results[enr_res_cilio.results['Adjusted P-value']<0.05]
print(enr_res_cilio_sig.shape)


# In[50]:


# plot geneset with enrichments
for geneset_name in ['GO_Biological_Process_2021']:
    fig, ax = plt.subplots(1,1, figsize=(10, 8))
    Enrichr_Plot(geneset_name,enr_res_cilio,'EP2')


# ### EP1+3

# In[51]:


glist_cilio = sc.get.rank_genes_groups_df(epith, group='EP1+3', key='wilcoxon_filtered_Epithelial_subcell_edited')[['names','logfoldchanges','pvals_adj','scores']]
glist_cilio = glist_cilio[~glist_cilio['names'].isnull()]
enr_res_cilio = gseapy.enrichr(gene_list=glist_cilio['names'][glist_cilio['pvals_adj']<0.05].tolist(),
                         organism='mouse',
                         gene_sets=['GO_Biological_Process_2021','GO_Molecular_Function_2021'],
                         description='pathway',
                         cutoff = 0.5,
                         outdir ='GoTerm_subcell/Enrichr/Epithelial/EP13')
# check if enrichment were found in all geneset
for geneset_name in ['GO_Biological_Process_2021','GO_Molecular_Function_2021']:
    df = enr_res_cilio.results[enr_res_cilio.results['Gene_set']==geneset_name]
    df_sig = df[df['Adjusted P-value']<0.05]
    
    if df_sig.shape[0]==0:
        print('No enrichments from %s!' %geneset_name)
        
    else:
        print('Oh YEAH there are enrichments in %s!' %geneset_name)
    
enr_res_cilio_sig = enr_res_cilio.results[enr_res_cilio.results['Adjusted P-value']<0.05]
print(enr_res_cilio_sig.shape)


# In[52]:


# plot geneset with enrichments
for geneset_name in ['GO_Biological_Process_2021']:
    fig, ax = plt.subplots(1,1, figsize=(10, 8))
    Enrichr_Plot(geneset_name,enr_res_cilio,'EP1+3')


# ## Prepare input for scCODA

# In[32]:


mouse.obs.CellType_edited4 = mouse.obs.CellType_edited4.cat.add_categories(['0','1','2','3'])


# In[33]:


for cell in epith.obs.index:
    mouse.obs.loc[cell,'CellType_edited4'] = epith.obs.loc[cell,'louvain_0.2']


# In[34]:


mouse.obs.CellType_edited4 = mouse.obs.CellType_edited4.cat.remove_categories('Epithelial')


# In[35]:


mouse.obs.CellType_edited4.replace({'2':'Ciliogenesis+prolif','3':'Gata3+','1':'EP2','0':'EP1+3'},inplace=True)


# In[36]:


sc.pl.umap(mouse,color=['Age','CellType_edited4'],wspace=0.5,size=3)


# In[37]:


mouse.obs.CellType_edited4.value_counts()


# ## scCODA

# In[38]:


cov_df = pd.DataFrame({"Age": ["P3", "P3", "P3",'P11','P11','P60','P60','P60']}, index=mouse.obs.batch.unique().tolist())
print(cov_df)


# In[39]:


data_scanpy_1 = dat.from_scanpy(
    mouse,
    cell_type_identifier="CellType_edited4",
    sample_identifier="batch",
    covariate_df=cov_df
)
print(data_scanpy_1)


# ## Grouped boxplots: relative abundance

# In[40]:


# Grouped boxplots. No facets, relative abundance, no dots.
viz.boxplots(
    data_scanpy_1,
    feature_name="Age",
    plot_facets=False,
    y_scale="relative",
    add_dots=False,
    figsize=(12,6),
    dpi=200
)
plt.show()


# In[41]:


# Grouped boxplots. Facets, log scale, added dots and custom color palette.
viz.boxplots(
    data_scanpy_1,
    feature_name="Age",
    plot_facets=True,
    y_scale="log",
    add_dots=True,
    cmap="Reds",
    dpi=200
)
plt.show()


# ## Re-run Composition Analysis by different reference levels, automatically chosen by scCODA

# In[42]:


viz.rel_abundance_dispersion_plot(
    data=data_scanpy_1,
    abundant_threshold=0.9,
    figsize=(8,8),dpi=100
)
plt.show()


# In[43]:


model_all1 = mod.CompositionalAnalysis(data_scanpy_1, formula="Age", reference_cell_type="automatic")
all_results1 = model_all1.sample_hmc()
all_results1.summary()


# **Gata3+ is the reference cell, EP1+3 and Ciliogenesis credible decreased and EP2 increased in P60**

# In[36]:


model_all2 = mod.CompositionalAnalysis(data_scanpy_1, formula="C(Age,Treatment('P3'))", reference_cell_type="Endothelial")
all_results2 = model_all2.sample_hmc()
all_results2.summary()


# In[39]:


all_results2.set_fdr(est_fdr=0.3)
all_results2.summary()


# **Endothelial is reference and FDR=0.3, Epithelial and Mesenchymal cells credible decreased in P60 Comparing to P3**

# In[ ]:





# In[ ]:





# In[ ]:




