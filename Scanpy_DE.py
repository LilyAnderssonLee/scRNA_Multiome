#!/usr/bin/env python
# coding: utf-8

# In[2]:


#pip install rpy2
import rpy2.robjects as robjects
import rpy2.robjects.lib.ggplot2 as ggplot2
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')


# In[1]:


#load packages
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy 
import gc
import os

sc.settings.verbosity = 3             
sc.settings.set_figure_params(dpi=80)


# In[2]:


adata = sc.read_h5ad('../Scanpy/results/scanpy_harmony_corrected_DE.h5ad')


# # Epithelial cells

# In[9]:


epith = adata[adata.obs['Canonical_marker_CT']=='Epithelial']


# ## Wilcoxon
# Output genes were ranked by Z-scores.

# In[10]:


sc.tl.rank_genes_groups(epith, 'Condition', method = 'wilcoxon', key_added = "wilcoxon_Condition")
epith.write_h5ad('../Scanpy/GSEA/Epithelial_DE.h5ad')


# ### min_logfoldchange=1
#  Here logfoldchange is log2 fold change for each gene for each group.

# In[26]:


sc.pl.rank_genes_groups_stacked_violin(epith, n_genes=40, key='wilcoxon_Condition', groupby='Condition',show=True,figsize=(20,4),min_logfoldchange=1)


# ### mim_logfoldchange=NULL

# In[27]:


sc.pl.rank_genes_groups_stacked_violin(epith, n_genes=40, key='wilcoxon_Condition', groupby='Condition',show=True,figsize=(20,4))


# In[8]:


def rank_genes_groups_df(adata_file, wilcoxon_key):
    # create a data frame with columns from .uns['rank_genes_groups'] (eg. names, 
    # logfoldchanges, pvals). 
    # Ideally, the list of columns should be consistent between methods
    # but 'logreg' does not return logfoldchanges for example

    dd = []
    groupby = adata_file.uns[wilcoxon_key]['params']['groupby']
    for group in adata_file.obs[groupby].cat.categories:
        cols = []
        # inner loop to make data frame by concatenating the columns per group
        for col in adata_file.uns[wilcoxon_key].keys():
            if col != 'params':
                   cols.append(pd.DataFrame(adata_file.uns[wilcoxon_key][col][group], columns=[col]))
        
        df = pd.concat(cols,axis=1)
        df['group'] = group
        dd.append(df)

    # concatenate the individual group data frames into one long data frame
    rgg = pd.concat(dd)
    rgg['group'] = rgg['group'].astype('category')
    return rgg.set_index('group')



# In[9]:


epith=sc.read_h5ad('../Scanpy/GSEA/Epithelial_DE.h5ad')
df_rank_genes = rank_genes_groups_df(epith,'wilcoxon_Condition')
#only statistically significant genes by pvals_adj<0.05 and abs(logfoldchange)>1
df_rank_genes = df_rank_genes.loc[(df_rank_genes['pvals_adj']<0.05) & (df_rank_genes['logfoldchanges'].abs()>0.25)]
df_rank_genes.to_csv('markers/Epithelial/Epithelial_Scapy_wilcoxon.csv')


# # Immune cells

# In[3]:


immune = adata[adata.obs['Canonical_marker_CT']=='Immune']


# In[6]:


sc.tl.rank_genes_groups(immune, 'Condition', method = 'wilcoxon', key_added = "wilcoxon_Condition")
immune.write_h5ad('../Scanpy/GSEA/Immune_DE.h5ad')


# In[7]:


#min_logfold=1
sc.pl.rank_genes_groups_stacked_violin(immune, n_genes=40, key='wilcoxon_Condition', groupby='Condition',show=True,figsize=(20,4),min_logfoldchange=1)
#no filter for logfc
sc.pl.rank_genes_groups_stacked_violin(immune, n_genes=40, key='wilcoxon_Condition', groupby='Condition',show=True,figsize=(20,4))


# In[9]:


df_rank_genes = rank_genes_groups_df(immune,'wilcoxon_Condition')
#only statistically significant genes by pvals_adj<0.05 and abs(logfoldchange)>0.25
df_rank_genes = df_rank_genes.loc[(df_rank_genes['pvals_adj']<0.05) & (df_rank_genes['logfoldchanges'].abs()>0.25)]
df_rank_genes.to_csv('markers/Immune/Immune_Scapy_wilcoxon.csv')


# # Mesenchymal cells

# In[12]:


mesen = adata[adata.obs['Canonical_marker_CT']=='Mesenchymal']
sc.tl.rank_genes_groups(mesen, 'Condition', method = 'wilcoxon', key_added = "wilcoxon_Condition")
mesen.write_h5ad('../Scanpy/GSEA/Mesenchymal_DE.h5ad')

#min_logfold=1
sc.pl.rank_genes_groups_stacked_violin(mesen, n_genes=40, key='wilcoxon_Condition', groupby='Condition',show=True,figsize=(20,4),min_logfoldchange=1)
#no filter for logfc
sc.pl.rank_genes_groups_stacked_violin(mesen, n_genes=40, key='wilcoxon_Condition', groupby='Condition',show=True,figsize=(20,4))


# In[13]:


df_rank_genes = rank_genes_groups_df(mesen,'wilcoxon_Condition')
#only statistically significant genes by pvals_adj<0.05 and abs(logfoldchange)>0.25
df_rank_genes = df_rank_genes.loc[(df_rank_genes['pvals_adj']<0.05) & (df_rank_genes['logfoldchanges'].abs()>0.25)]
df_rank_genes.to_csv('markers/Mesenchymal/Mesenchymal_Scapy_wilcoxon.csv')


# # NeuronGlial cells

# In[16]:


neuroglial = adata[(adata.obs['Canonical_marker_CT']=='Oligos') | (adata.obs['Canonical_marker_CT']=='Neurons') | (adata.obs['Canonical_marker_CT']=='Astrocyte')]
sc.tl.rank_genes_groups(neuroglial, 'Condition', method = 'wilcoxon', key_added = "wilcoxon_Condition")
neuroglial.write_h5ad('../Scanpy/GSEA/NeuronGlials_DE.h5ad')

#min_logfold=1
sc.pl.rank_genes_groups_stacked_violin(neuroglial, n_genes=40, key='wilcoxon_Condition', groupby='Condition',show=True,figsize=(20,4),min_logfoldchange=1)
#no filter for logfc
sc.pl.rank_genes_groups_stacked_violin(neuroglial, n_genes=40, key='wilcoxon_Condition', groupby='Condition',show=True,figsize=(20,4))



# In[17]:


df_rank_genes = rank_genes_groups_df(neuroglial,'wilcoxon_Condition')
#only statistically significant genes by pvals_adj<0.05 and abs(logfoldchange)>0.25
df_rank_genes = df_rank_genes.loc[(df_rank_genes['pvals_adj']<0.05) & (df_rank_genes['logfoldchanges'].abs()>0.25)]
df_rank_genes.to_csv('markers/NeuronGlial/NeuronGlial_Scapy_wilcoxon.csv')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




