#!/usr/bin/env python
# coding: utf-8

#load packages
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import gc
import os
from matplotlib.pyplot import rc_context
import gseapy

sc.settings.verbosity = 3             
sc.settings.set_figure_params(dpi=100,figsize=(4,4))


# # Load processed adata


#all cell
adata_raw = sc.read_h5ad('../../results/Human_filtered_rawCounts.h5ad')
adata = sc.read_h5ad('../../results/Human_harmony_CellType.h5ad')
#subcell
epith = sc.read_h5ad('../results/Epithelial_har_noDbl_1.h5ad')
endoth = sc.read_h5ad('../results/Endothelial_har_noDbl_1.h5ad')
immune = sc.read_h5ad('../results/Immune_har_noDbl_1.h5ad')
mesen = sc.read_h5ad('../results/Mesenchymal_har_noDbl_1.h5ad')
neuroglia = sc.read_h5ad('../results/NeuroGlia_har_noDbl_1.h5ad')


adata.obs['CellType_intergrate_edited6'] = adata.obs['CellType_intergrate_edited5']
adata.obs['CellType_intergrate_edited6'].replace({'Epithelial-Immune':'Doublets','NeuroGlia-Epithelial':'Doublets','Epithelial-Immune-Neuron':'Doublets',
                                                 'Mesenchymal-Epithelial':'Doublets','Mesenchymal-Immune':'Doublets','Endothelial-Immune':'Doublets',
                                                  'NeuroGlia-Immune':'Doublets','Endothelial-Epithelial':'Doublets'},inplace=True)


# # Add NeuroGlia, Immune, and Mesenchymal info to adata


adata.obs['CellType_Subcell'] = adata.obs['CellType_intergrate_edited6']
#NeuroGlia
adata.obs['CellType_Subcell'] = adata.obs['CellType_Subcell'].cat.add_categories(['Oligo','Astrocyte','Neuron1','Neuron2','OPC'])
for cell in neuroglia.obs.index:
    adata.obs.loc[cell,'CellType_Subcell']=neuroglia.obs.loc[cell,'NeuroGlia_subcell']
#Immune
adata.obs['CellType_Subcell'] = adata.obs['CellType_Subcell'].cat.add_categories(['Macrophage','T cells','Epiplexus','NK cells','B cells','Neutrophil/Monocyte','Macrophage Proliferating'])
for cell in immune.obs.index:
    adata.obs.loc[cell,'CellType_Subcell']=immune.obs.loc[cell,'Immune_subcell']
#Mesenchymal
adata.obs['CellType_Subcell'] = adata.obs['CellType_Subcell'].cat.add_categories(['FB2','FB1','vSMC','Pericytes'])
for cell in mesen.obs.index:
    adata.obs.loc[cell,'CellType_Subcell']=mesen.obs.loc[cell,'Mesenchymal_subcells']


# In[34]:


adata.obs['CellType_Subcell'].replace({'FB1':'Fibroblast','FB2':'Fibroblast','T cells':'Lymphocyte','NK cells':'Lymphocyte','B cells':'Lymphocyte',
                                      'Neutrophil/Monocyte':'Macrophage','Neutrophil/Monocyte':'Macrophage','Macrophage Proliferating':'Macrophage',
                                      'Neuron1':'Neurons','Neuron2':'Neurons','Epiplexus':'Macrophage'},inplace=True)

adata.obs['CellType_Subcell']=adata.obs['CellType_Subcell'].cat.remove_categories(['Mesenchymal','NeuronGlial','Immune'])



adata.obs['CellType_Subcell'].cat.reorder_categories(['Epithelial', 'Endothelial','Ependymal','Macrophage','Lymphocyte','Fibroblast','Pericytes','vSMC',
                                                     'Astrocyte','Oligo','OPC','Neurons','Doublets'], inplace=True)


sc.pl.umap(adata,color=['Condition','CellType_Subcell','Gender','batch'],title=['Condition','Cell Types','Gender','Batch'],ncols=2,wspace=0.3)


# # Add subcell info to adata


meta_data = adata.obs
meta_data = pd.merge(meta_data, neuroglia.obs.loc[:,'NeuroGlia_subcell'], left_index=True, right_index=True, how='left')
meta_data = pd.merge(meta_data, immune.obs.loc[:,'Immune_subcell'], left_index=True, right_index=True, how='left')
meta_data = pd.merge(meta_data, epith.obs.loc[:,'Epithelial_Subcell'], left_index=True, right_index=True, how='left')
meta_data = pd.merge(meta_data, mesen.obs.loc[:,'Mesenchymal_subcells'], left_index=True, right_index=True, how='left')
meta_data = pd.merge(meta_data, endoth.obs.loc[:,'Endothelial_Subcells'], left_index=True, right_index=True, how='left')
#update meta data in adata
adata.obs = meta_data


# ## Visualize cell types


sc.pl.umap(adata,color='CellType_Subcell',title='Cell Types')
sc.pl.umap(adata,color='NeuroGlia_subcell',title='NeuroGlia subcell')
sc.pl.umap(adata,color='Immune_subcell',title='Immune subcell')
sc.pl.umap(adata,color='Epithelial_Subcell',title='Epithelial subcell')
sc.pl.umap(adata,color='Mesenchymal_subcells',title='Mesenchymal subcell')
sc.pl.umap(adata,color='Endothelial_Subcells',title='Endothelial subcell')


# ## Visualize Subcells


sc.pl.umap(neuroglia,color='NeuroGlia_subcell')
sc.pl.umap(mesen,color='Mesenchymal_subcells')
sc.pl.umap(immune,color='Immune_subcell')
sc.pl.umap(epith,color='Epithelial_Subcell')
sc.pl.umap(endoth,color='Endothelial_Subcells')


# ## Top markers


sc.tl.rank_genes_groups(adata,'CellType_Subcell',method='wilcoxon',key_added='wilcoxon_CellType_Subcell')
sc.tl.dendrogram(adata,groupby='CellType_Subcell')
sc.tl.filter_rank_genes_groups(adata, groupby='CellType_Subcell', key = "wilcoxon_CellType_Subcell", key_added='wilcoxon_filtered_CellType_Subcell',
                               min_in_group_fraction=0.1,max_out_group_fraction=0.5,min_fold_change=0.5)



sc.pl.rank_genes_groups_dotplot(adata,n_genes=10,key='wilcoxon_filtered_CellType_Subcell', groupby='CellType_Subcell',show=True)


sc.pl.rank_genes_groups_dotplot(adata,n_genes=10,key='wilcoxon_filtered_CellType_Subcell', groupby='CellType_Subcell',show=True,
                               groups = ['Epithelial','Macrophage','Lymphocyte','Endothelial','Fibroblast','Pericytes','vSMC','Oligo','OPC','Neurons','Ependymal','Astrocyte'])


# ## Cell Summary


fig, (ax1,ax2) = plt.subplots(nrows=1,ncols=2,constrained_layout=True,figsize=(12,4))
sns.countplot(y = 'batch',data = adata.obs,ax=ax1)
for container in ax1.containers:
    ax1.bar_label(container)
    
sns.countplot(y = 'Condition',data = adata.obs,ax=ax2)
for container in ax2.containers:
    ax2.bar_label(container)

fig, (ax1) = plt.subplots(nrows=1,ncols=1,constrained_layout=True,figsize=(8,4))
sns.countplot(y = 'CellType_Subcell',data = adata.obs,ax=ax1)
for container in ax1.containers:
    ax1.bar_label(container)

def batch_pct_cellType(adata_file,groupby_key,legend_key,plot_title,pos):
    tmp = pd.crosstab(adata_file.obs[groupby_key],adata_file.obs[legend_key], normalize='index')   
    tmp.plot.bar(stacked=True,figsize=(12,15),title=plot_title, rot=90, grid=False,
                 xlabel='',ax=pos).legend(bbox_to_anchor=[1.2,1])
    
fig, (ax1,ax2,ax3) = plt.subplots(3, 1, constrained_layout=True)
batch_pct_cellType(adata,'CellType_Subcell','batch','',ax1)
batch_pct_cellType(adata,'CellType_Subcell','Condition','',ax2)
batch_pct_cellType(adata,'Condition','CellType_Subcell','',ax3)


# ## Remove doublets

adata_noDbl = adata[adata.obs['CellType_Subcell']!='Doublets']


# ### Top markers

sc.tl.rank_genes_groups(adata_noDbl,'CellType_Subcell',method='wilcoxon',key_added='wilcoxon_CellType_Subcell')
sc.tl.dendrogram(adata_noDbl,groupby='CellType_Subcell')
sc.tl.filter_rank_genes_groups(adata_noDbl, groupby='CellType_Subcell', key = "wilcoxon_CellType_Subcell", key_added='wilcoxon_filtered_CellType_Subcell',
                               min_in_group_fraction=0.1,max_out_group_fraction=0.5,min_fold_change=0.5)
sc.pl.rank_genes_groups_dotplot(adata_noDbl,n_genes=10,key='wilcoxon_filtered_CellType_Subcell', groupby='CellType_Subcell',show=True)


def rank_genes_groups_df(adata_file, wilcoxon_key):
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

df = rank_genes_groups_df(adata_noDbl,'wilcoxon_filtered_CellType_Subcell')
df = df[df['pvals_adj']<0.05]
df = df[~df['names'].isnull()]

#only keep the top 100 genes for each cluster. Genes were sorted by the Score. 
new_df = pd.DataFrame([], columns = df.columns)
for index in np.unique(df.index):
    new_df = pd.concat([new_df, df.loc[index].head(n=200)])
    new_df.to_csv('../../results/Final_CellType_markers_filtered.csv')


# ### Cell Summary


fig, (ax1,ax2) = plt.subplots(nrows=1,ncols=2,constrained_layout=True,figsize=(12,4))
sns.countplot(y = 'batch',data = adata_noDbl.obs,ax=ax1)
for container in ax1.containers:
    ax1.bar_label(container)
    
sns.countplot(y = 'Condition',data = adata_noDbl.obs,ax=ax2)
for container in ax2.containers:
    ax2.bar_label(container)

fig, (ax1) = plt.subplots(nrows=1,ncols=1,constrained_layout=True,figsize=(8,4))
sns.countplot(y = 'CellType_Subcell',data = adata_noDbl.obs,ax=ax1)
for container in ax1.containers:
    ax1.bar_label(container)

def batch_pct_cellType(adata_file,groupby_key,legend_key,plot_title,pos):
    tmp = pd.crosstab(adata_file.obs[groupby_key],adata_file.obs[legend_key], normalize='index')   
    tmp.plot.bar(stacked=True,figsize=(12,15),title=plot_title, rot=90, grid=False,
                 xlabel='',ax=pos).legend(bbox_to_anchor=[1.25,1])
    
fig, (ax1,ax2,ax3) = plt.subplots(3, 1, constrained_layout=True)
batch_pct_cellType(adata_noDbl,'CellType_Subcell','batch','',ax1)
batch_pct_cellType(adata_noDbl,'CellType_Subcell','Condition','',ax2)
batch_pct_cellType(adata_noDbl,'Condition','CellType_Subcell','',ax3)


# # Save the adata

#remove redundant info
wilcoxon_key_list=[]
for key in adata.uns.keys():
    if key.startswith('wilcoxon_filtered_'):
        wilcoxon_key_list.append(key)

for i in wilcoxon_key_list:
    del adata.uns[i]
    
adata.write_h5ad('../results/Human_harmony_CellType_edited.h5ad')


#remove redundant info
wilcoxon_key_list=[]
for key in adata_noDbl.uns.keys():
    if key.startswith('wilcoxon_filtered_'):
        wilcoxon_key_list.append(key)

for i in wilcoxon_key_list:
    del adata_noDbl.uns[i]
    
adata_noDbl.write_h5ad('../results/Human_harmony_CellType_noDoublets_edited.h5ad')


# # Update adata_raw meta info


adata_raw.obs = adata.obs
adata_raw.write_h5ad('../../results/Human_filtered_rawCounts_edited.h5ad')


# # Subcells

# ## Top markers

# Epithelial cells


#del epith.uns['log1p']
sc.tl.rank_genes_groups(epith,'Epithelial_Subcell',method='wilcoxon',key_added='wilcoxon_Epithelial_Subcell')
sc.tl.dendrogram(epith,groupby='Epithelial_Subcell')
sc.tl.filter_rank_genes_groups(epith, groupby='Epithelial_Subcell', key = "wilcoxon_Epithelial_Subcell", key_added='wilcoxon_filtered_Epithelial_Subcell',
                               min_in_group_fraction=0.1,max_out_group_fraction=0.5,min_fold_change=0.5)


# Endothelial cells


#del endoth.uns['log1p']
sc.tl.rank_genes_groups(endoth,'Endothelial_Subcells',method='wilcoxon',key_added='wilcoxon_Endothelial_Subcells')
sc.tl.dendrogram(endoth,groupby='Endothelial_Subcells')
sc.tl.filter_rank_genes_groups(endoth, groupby='Endothelial_Subcells', key = "wilcoxon_Endothelial_Subcells", key_added='wilcoxon_filtered_Endothelial_Subcells',
                               min_in_group_fraction=0.1,max_out_group_fraction=0.5,min_fold_change=0.5)


# Mesenchymal cells


del mesen.uns['log1p']
sc.tl.rank_genes_groups(mesen,'Mesenchymal_subcells',method='wilcoxon',key_added='wilcoxon_Mesenchymal_subcells')
sc.tl.dendrogram(mesen,groupby='Mesenchymal_subcells')
sc.tl.filter_rank_genes_groups(mesen, groupby='Mesenchymal_subcells', key = "wilcoxon_Mesenchymal_subcells", key_added='wilcoxon_filtered_Mesenchymal_subcells',
                               min_in_group_fraction=0.1,max_out_group_fraction=0.5,min_fold_change=0.5)


# Immune cells


del immune.uns['log1p']
sc.tl.rank_genes_groups(immune,'Immune_subcell',method='wilcoxon',key_added='wilcoxon_Immune_subcell')
sc.tl.dendrogram(immune,groupby='Immune_subcell')
sc.tl.filter_rank_genes_groups(immune, groupby='Immune_subcell', key = "wilcoxon_Immune_subcell", key_added='wilcoxon_filtered_Immune_subcell',
                               min_in_group_fraction=0.1,max_out_group_fraction=0.5,min_fold_change=0.5)


# NeuroGlia cells

del neuroglia.uns['log1p']
sc.tl.rank_genes_groups(neuroglia,'NeuroGlia_subcell',method='wilcoxon',key_added='wilcoxon_NeuroGlia_subcell')
sc.tl.dendrogram(neuroglia,groupby='NeuroGlia_subcell')
sc.tl.filter_rank_genes_groups(neuroglia, groupby='NeuroGlia_subcell', key = "wilcoxon_NeuroGlia_subcell", key_added='wilcoxon_filtered_NeuroGlia_subcell',
                               min_in_group_fraction=0.1,max_out_group_fraction=0.5,min_fold_change=0.5)


# ## DotPlot


sc.pl.rank_genes_groups_dotplot(epith,n_genes=15,key='wilcoxon_filtered_Epithelial_Subcell', groupby='Epithelial_Subcell',show=True)
sc.pl.rank_genes_groups_dotplot(endoth,n_genes=15,key='wilcoxon_filtered_Endothelial_Subcells', groupby='Endothelial_Subcells',show=True)
sc.pl.rank_genes_groups_dotplot(mesen,n_genes=15,key='wilcoxon_filtered_Mesenchymal_subcells', groupby='Mesenchymal_subcells',show=True)
sc.pl.rank_genes_groups_dotplot(immune,n_genes=15,key='wilcoxon_filtered_Immune_subcell', groupby='Immune_subcell',show=True)
sc.pl.rank_genes_groups_dotplot(neuroglia,n_genes=15,key='wilcoxon_filtered_NeuroGlia_subcell', groupby='NeuroGlia_subcell',show=True)

