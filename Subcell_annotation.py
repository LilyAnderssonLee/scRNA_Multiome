#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
pd.set_option("max_rows", None)
pd.set_option("max_columns", None)


# #  raw counts

# In[2]:


adata_raw = sc.read_h5ad('../../results/Mouse_filtered_rawcounts.h5ad')


# In[3]:


adata_raw.obs['CellType'].value_counts()


# #  With doubelts

# In[4]:


epith = adata_raw[(adata_raw.obs['CellType']=='Epithelial') | (adata_raw.obs['CellType']=='Ciliogenesis')|
                 (adata_raw.obs['CellType']=='Epith-Endoth-Mesen') | (adata_raw.obs['CellType']=='Epith-Immune')|
                 (adata_raw.obs['CellType']=='Epith-Astro')]
sc.pp.filter_genes(epith, min_cells=5)
sc.pl.violin(epith,keys=['total_counts','n_genes_by_counts','pct_counts_Mito','pct_counts_Ribo'],rotation=90,stripplot=False,
             jitter=0.4,color = "#006837",multi_panel=True,groupby='batch')


# ## Normalization, DR, integration & clustering

# In[5]:


# normalize to depth 10000
sc.pp.normalize_per_cell(epith,counts_per_cell_after=1e4)
#logaritmize
sc.pp.log1p(epith)
epith.raw = epith

sc.pp.highly_variable_genes(epith, min_mean=0.0125, max_mean=3, min_disp=0.5)
var_genes_all = epith.var.highly_variable
print("Highly variable genes: %d"%sum(var_genes_all))

#Detect variable genes in each dataset separately using the batch_key parameter.
sc.pp.highly_variable_genes(epith, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key = 'batch')
print("Highly variable genes intersection: %d"%sum(epith.var.highly_variable_intersection))

print("Number of batches where gene is variable:")
print(epith.var.highly_variable_nbatches.value_counts())

var_genes_batch = epith.var.highly_variable_nbatches > 0

#Compare overlap of the variable genes.
print("Any batch var genes: %d"%sum(var_genes_batch))
print("All data var genes: %d"%sum(var_genes_all))
print("Overlap: %d"%sum(var_genes_batch & var_genes_all))
print("Variable genes in all batches: %d"%sum(epith.var.highly_variable_nbatches == 10))
print("Overlap batch instersection and all: %d"%sum(var_genes_all & epith.var.highly_variable_intersection))

#Select all genes that are variable in at least 4 datasets and use for remaining analysis.
var_select = epith.var.highly_variable_nbatches >4
var_genes = var_select.index[var_select]
len(var_genes)

epith= epith[:,var_genes]
#sc.pp.regress_out(epith,['total_counts','n_genes_by_counts','pct_counts_Mito','pct_counts_Ribo'])
sc.pp.scale(epith,max_value=10)

sc.tl.pca(epith,svd_solver = 'arpack',n_comps = 50)
sc.pl.pca_variance_ratio(epith,n_pcs=50)
#umap
sc.pp.neighbors(epith,n_pcs = 10)
sc.tl.umap(epith)
sc.tl.tsne(epith,n_pcs=20,random_state=0)
epith.write_h5ad('../results/Epithelial_withBE.h5ad')

#harmony
import harmonypy
import scanpy.external as sce
epith_har = epith
sce.pp.harmony_integrate(epith_har, 'batch',max_iter_harmony=50)
# tsne and umap
sc.pp.neighbors(epith_har, n_pcs =50, use_rep = "X_pca_harmony")
sc.tl.umap(epith_har)
sc.tl.tsne(epith_har, n_pcs = 50, use_rep = "X_pca_harmony")
epith_har.write_h5ad('../results/Epithelial_har.h5ad')

#remove redundant info
wilcoxon_key_list=[]
for key in epith_har.uns.keys():
    if key.startswith('wilcoxon_louvain_'):
        wilcoxon_key_lwspace=t.append(key)

for i in wilcoxon_key_list:
    del epith_har.uns[i]

for name in epith_har.obs.columns:
    if name.startswith('louvain'):
        del epith_har.obs[name]
        
res = [0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0]
for snn_res in res:
    sc.tl.louvain(epith_har, resolution = snn_res, key_added = 'louvain_%s' %snn_res) 

#Top markers of louvain snn
res_list = epith_har.obs.columns[epith_har.obs.columns.str.startswith('louvain')]
for snn_res in res_list:
    sc.tl.rank_genes_groups(epith_har,snn_res,method='wilcoxon',key_added='wilcoxon_%s' %snn_res)
    sc.tl.dendrogram(epith_har,groupby=snn_res)
    sc.tl.filter_rank_genes_groups(epith_har, groupby=snn_res, key = "wilcoxon_%s" %snn_res, key_added='wilcoxon_filtered_%s' %snn_res,
                                   min_in_group_fraction=0.25,max_out_group_fraction=0.5,min_fold_change=0.25)
    


# In[7]:


for snn_res in res_list[:4]:
    sc.pl.rank_genes_groups_dotplot(epith_har,n_genes=10,key='wilcoxon_filtered_%s' %snn_res, groupby=snn_res,show=True)
    
#plot
sc.pl.umap(epith_har,color=res_list,ncols=4, wspace=0.3)
sc.pl.umap(epith_har,color=['CellType','Age','batch'],wspace=0.6)


# # Doublets check

# In[8]:


epithelial = ["Krt18","Htr2c","Car12","Clic6", "Folr1","Igfbp2", "Ttr"]
mesenchymal = ["Col3a1", "Col1a2", "Lum","Col1a1",'Rgs5','Des','Pdgfra','Pdgfrb']
endothelial = ["Egfl7","Kdr","Pecam1", "Cd34",  "Plvap"]
immune = ["C1qa","Coro1a", "Aif1", "Tyrobp", "Fcer1g","Cx3cr1",'Cd74','Cd68','Skap1','Ptprc','Txk']
neurons = ["Rtn1","Ntm"]
glia_like = ["Bcan","Plcd4","Slc1a3",'Gfap','Fa2h','Cnp']
ciliogenesis = ["Deup1","Shisa8","Ccno","Mcidas"]
pericyte = ['Acta2','Des']


# In[9]:


sc.pl.violin(epith_har,keys=['total_counts','n_genes_by_counts'],rotation=90,stripplot=False,jitter=0.4,color = "#006837",multi_panel=True,groupby='CellType')
sc.pl.umap(epith_har,color=['Clic6','Folr1','Acta2','Des','Col3a1','Pecam1','Cd34','C1qa','Ptprc','Gfap','Cnp','Deup1'],ncols=4)


# # Remove doublets

# In[10]:


epith_noDoublets = epith_har[(epith_har.obs['CellType']=='Epithelial') |(epith_har.obs['CellType']=='Ciliogenesis') ]

epith1 = adata_raw[epith_noDoublets.obs.index]
sc.pp.filter_genes(epith1, min_cells=5)
sc.pl.violin(epith1,keys=['total_counts','n_genes_by_counts','pct_counts_Mito','pct_counts_Ribo'],rotation=90,stripplot=False,
             jitter=0.4,color = "#006837",multi_panel=True,groupby='batch')


# ## Normalization, DR, integration & clustering

# In[11]:


# normalize to depth 10000
sc.pp.normalize_per_cell(epith1,counts_per_cell_after=1e4)
#logaritmize
sc.pp.log1p(epith1)
epith1.raw = epith1

sc.pp.highly_variable_genes(epith1, min_mean=0.0125, max_mean=3, min_disp=0.5)
var_genes_all = epith1.var.highly_variable
print("Highly variable genes: %d"%sum(var_genes_all))

#Detect variable genes in each dataset separately using the batch_key parameter.
sc.pp.highly_variable_genes(epith1, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key = 'batch')
print("Highly variable genes intersection: %d"%sum(epith1.var.highly_variable_intersection))

print("Number of batches where gene is variable:")
print(epith1.var.highly_variable_nbatches.value_counts())

var_genes_batch = epith1.var.highly_variable_nbatches > 0

#Compare overlap of the variable genes.

print("Any batch var genes: %d"%sum(var_genes_batch))
print("All data var genes: %d"%sum(var_genes_all))
print("Overlap: %d"%sum(var_genes_batch & var_genes_all))
print("Variable genes in all batches: %d"%sum(epith1.var.highly_variable_nbatches == 10))
print("Overlap batch instersection and all: %d"%sum(var_genes_all & epith1.var.highly_variable_intersection))

#Select all genes that are variable in at least 4 datasets and use for remaining analysis.
var_select = epith1.var.highly_variable_nbatches >4
var_genes = var_select.index[var_select]
len(var_genes)

epith1= epith1[:,var_genes]
#sc.pp.regress_out(epith1,['total_counts','n_genes_by_counts','pct_counts_Mito','pct_counts_Ribo'])
sc.pp.scale(epith1,max_value=10)

sc.tl.pca(epith1,svd_solver = 'arpack',n_comps = 50)
sc.pl.pca_variance_ratio(epith1,n_pcs=50)
#umap
sc.pp.neighbors(epith1,n_pcs = 10, n_neighbors=15)
sc.tl.umap(epith1)
sc.tl.tsne(epith1,n_pcs=20,random_state=0)
epith1.write_h5ad('../results/Epithelial_withBE_noDbl.h5ad')

#harmony
import harmonypy
import scanpy.external as sce
epith_har1 = epith1
sce.pp.harmony_integrate(epith_har1, 'batch',max_iter_harmony=20)
# tsne and umap
sc.pp.neighbors(epith_har1, n_pcs =50, use_rep = "X_pca_harmony")
sc.tl.umap(epith_har1)
sc.tl.tsne(epith_har1, n_pcs = 50, use_rep = "X_pca_harmony")
epith_har1.write_h5ad('../results/Epithelial_har_noDbl.h5ad')


# In[12]:


#remove redundant info
wilcoxon_key_list=[]
for key in epith_har1.uns.keys():
    if key.startswith('wilcoxon_louvain_'):
        wilcoxon_key_list.append(key)

for i in wilcoxon_key_list:
    del epith_har1.uns[i]

for name in epith_har1.obs.columns:
    if name.startswith('louvain'):
        del epith_har1.obs[name]
        
res = [0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0]
for snn_res in res:
    sc.tl.louvain(epith_har1, resolution = snn_res, key_added = 'louvain_%s' %snn_res) 

#Top markers of louvain snn
res_list = epith_har1.obs.columns[epith_har1.obs.columns.str.startswith('louvain')]
for snn_res in res_list:
    sc.tl.rank_genes_groups(epith_har1,snn_res,method='wilcoxon',key_added='wilcoxon_%s' %snn_res)
    sc.tl.dendrogram(epith_har1,groupby=snn_res)
    sc.tl.filter_rank_genes_groups(epith_har1, groupby=snn_res, key = "wilcoxon_%s" %snn_res, key_added='wilcoxon_filtered_%s' %snn_res,
                                   min_in_group_fraction=0.25,max_out_group_fraction=0.5,min_fold_change=0.25)


# In[13]:


for snn_res in res_list[:4]:
    sc.pl.rank_genes_groups_dotplot(epith_har1,n_genes=10,key='wilcoxon_filtered_%s' %snn_res, groupby=snn_res,show=True)
    
#plot
sc.pl.umap(epith_har1,color=res_list,ncols=4, wspace=0.3)
sc.pl.umap(epith_har1,color=['CellType','Age','batch'],wspace=0.6)


# In[14]:


#remove redundant info
wilcoxon_key_list=[]
for key in epith_har1.uns.keys():
    if key.startswith('wilcoxon_filtered_'):
        wilcoxon_key_list.append(key)

for i in wilcoxon_key_list:
    del epith_har1.uns[i]
    
epith_har1.write_h5ad('../results/Epithelial_har_noDbl.h5ad')


# In[15]:


sc.pl.umap(epith_har1,color = ['Acta2','Myl9','Cnn1',
                               'louvain_0.5','Age','batch'],ncols=3)#myoepithelial from human organoid


# !!!! Check if the cells of cluster 8 under louvain_0.5 are mesenchymal-epithelial doublets

# In[16]:


sc.pl.umap(epith_har1,color=['total_counts','n_genes_by_counts','pct_counts_Mito','pct_counts_Ribo'])
sc.pl.violin(epith_har1,keys=['total_counts','n_genes_by_counts'],rotation=90,stripplot=False,jitter=0.4,color = "#006837",
             multi_panel=True,groupby='louvain_0.5')


# ## Doublets recheck

# In[17]:


sc.pl.umap(epith_har1,color=['Clic6','Folr1','Acta2','Des','Col3a1','Pecam1','Cd34','Mki67','Top2a','Cdk1','Shisa8','Deup1'],ncols=4)


# !!!! If combine the violin plot and the UMAP plots, it looks like C7 of louvain_0.5 is doublets. it could be a good idea to validate.

# In[20]:


sc.pl.umap(epith_har1,color='Gata3')


# In[39]:


#myoepithelial cell markers from Laura organoids
sc.pl.umap(epith_har1,color = ['Acta2','Myl9'])


# # Subcells_Version1
# Include cells have both Epithelial and Mesenchymal features.

# In[22]:


epith_har1 = sc.read_h5ad('../results/Epithelial_har_noDbl.h5ad')
epith_har1.obs['Epithelial_subcell']  = epith_har1.obs['louvain_0.5']

#epith_har1.obs['Epithelial_subcell'] = epith_har1.obs['Epithelial_subcell'].cat.add_categories(['Myoepithelial'])
#for cell_id in epith_har1.obs.index[epith_har1.obs['louvain_0.5']=='7']:
    #epith_har1.obs.loc[cell_id,'Epithelial_subcell']='Myoepithelial'
#epith_har1.obs['Epithelial_subcell']  = epith_har1.obs['Epithelial_subcell'].astype('category')

epith_har1.obs['Epithelial_subcell'].replace({'0':'EP','1':'EP','2':'EP','3':'EP','4':'EP',
                                              '5':'Ciliogenesis/EC Proliferating','6':'Gata3+Pon1','7':'Myoepithelial'},inplace=True)


# In[23]:


sc.pl.umap(epith_har1,color=['Epithelial_subcell','Age','batch'],wspace=0.7)


# In[24]:


sc.tl.rank_genes_groups(epith_har1,'Epithelial_subcell',method='wilcoxon',key_added='wilcoxon_Epithelial_subcell')
sc.tl.dendrogram(epith_har1,groupby='Epithelial_subcell')
sc.tl.filter_rank_genes_groups(epith_har1, groupby='Epithelial_subcell', key = "wilcoxon_Epithelial_subcell", key_added='wilcoxon_filtered_Epithelial_subcell',
                               min_in_group_fraction=0.1,max_out_group_fraction=0.5,min_fold_change=0.25)
sc.pl.rank_genes_groups_dotplot(epith_har1,n_genes=15,key='wilcoxon_filtered_Epithelial_subcell', groupby='Epithelial_subcell',show=True)


# ## Save the data

# In[25]:


#remove redundant info
wilcoxon_key_list=[]
for key in epith_har1.uns.keys():
    if key.startswith('wilcoxon_filtered_'):
        wilcoxon_key_list.append(key)

for i in wilcoxon_key_list:
    del epith_har1.uns[i]
    
epith_har1.write_h5ad('../results/Epithelial_har_noDbl_1.h5ad')


# # Remove Myoepithelial
# 
# Endothelial, vSMC and prolifetating cell markers also expressed.

# In[27]:


epith_noDoublets2 = epith_har1[(epith_har1.obs['louvain_0.5']!='7')]
epith2 = adata_raw[epith_noDoublets2.obs.index]
sc.pp.filter_genes(epith2, min_cells=5)
sc.pl.violin(epith2,keys=['total_counts','n_genes_by_counts','pct_counts_Mito','pct_counts_Ribo'],rotation=90,stripplot=False,
             jitter=0.4,color = "#006837",multi_panel=True,groupby='batch')


# ## Normalization, DR, integration & clustering

# In[28]:


# normalize to depth 10000
sc.pp.normalize_per_cell(epith2,counts_per_cell_after=1e4)
#logaritmize
sc.pp.log1p(epith2)
epith2.raw = epith2

sc.pp.highly_variable_genes(epith2, min_mean=0.0125, max_mean=3, min_disp=0.5)
var_genes_all = epith2.var.highly_variable
print("Highly variable genes: %d"%sum(var_genes_all))

#Detect variable genes in each dataset separately using the batch_key parameter.
sc.pp.highly_variable_genes(epith2, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key = 'batch')
print("Highly variable genes intersection: %d"%sum(epith2.var.highly_variable_intersection))

print("Number of batches where gene is variable:")
print(epith2.var.highly_variable_nbatches.value_counts())

var_genes_batch = epith2.var.highly_variable_nbatches > 0

#Compare overlap of the variable genes.

print("Any batch var genes: %d"%sum(var_genes_batch))
print("All data var genes: %d"%sum(var_genes_all))
print("Overlap: %d"%sum(var_genes_batch & var_genes_all))
print("Variable genes in all batches: %d"%sum(epith2.var.highly_variable_nbatches == 10))
print("Overlap batch instersection and all: %d"%sum(var_genes_all & epith2.var.highly_variable_intersection))

#Select all genes that are variable in at least 4 datasets and use for remaining analysis.
var_select = epith2.var.highly_variable_nbatches >4
var_genes = var_select.index[var_select]
len(var_genes)

epith2= epith2[:,var_genes]
#sc.pp.regress_out(epith2,['total_counts','n_genes_by_counts','pct_counts_Mito','pct_counts_Ribo'])
sc.pp.scale(epith2,max_value=10)

sc.tl.pca(epith2,svd_solver = 'arpack',n_comps = 50)
sc.pl.pca_variance_ratio(epith2,n_pcs=50)
#umap
sc.pp.neighbors(epith2,n_pcs = 10, n_neighbors=15)
sc.tl.umap(epith2)
sc.tl.tsne(epith2,n_pcs=20,random_state=0)
epith2.write_h5ad('../results/Epithelial_withBE_noDbl_2.h5ad')

#harmony
import harmonypy
import scanpy.external as sce
epith_har2 = epith2
sce.pp.harmony_integrate(epith_har2, 'batch',max_iter_harmony=20)
# tsne and umap
sc.pp.neighbors(epith_har2, n_pcs =50, use_rep = "X_pca_harmony")
sc.tl.umap(epith_har2)
sc.tl.tsne(epith_har2, n_pcs = 50, use_rep = "X_pca_harmony")
epith_har2.write_h5ad('../results/Epithelial_har_noDbl_2.h5ad')


# In[55]:


#remove redundant info
wilcoxon_key_list=[]
for key in epith_har2.uns.keys():
    if key.startswith('wilcoxon_louvain_'):
        wilcoxon_key_list.append(key)

for i in wilcoxon_key_list:
    del epith_har2.uns[i]

for name in epith_har2.obs.columns:
    if name.startswith('louvain'):
        del epith_har2.obs[name]
        
res = [0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.0]
for snn_res in res:
    sc.tl.louvain(epith_har2, resolution = snn_res, key_added = 'louvain_%s' %snn_res) 

#Top markers of louvain snn
res_list = epith_har2.obs.columns[epith_har2.obs.columns.str.startswith('louvain')]
for snn_res in res_list:
    sc.tl.rank_genes_groups(epith_har2,snn_res,method='wilcoxon',key_added='wilcoxon_%s' %snn_res)
    sc.tl.dendrogram(epith_har2,groupby=snn_res)
    sc.tl.filter_rank_genes_groups(epith_har2, groupby=snn_res, key = "wilcoxon_%s" %snn_res, key_added='wilcoxon_filtered_%s' %snn_res,
                                   min_in_group_fraction=0.1,max_out_group_fraction=0.5,min_fold_change=1)


# In[56]:


for snn_res in res_list:
    sc.pl.rank_genes_groups_dotplot(epith_har2,n_genes=10,key='wilcoxon_filtered_%s' %snn_res, groupby=snn_res,show=True)
#plot
sc.pl.umap(epith_har2,color=res_list,ncols=4, wspace=0.3)


# ## Save top markers of louvian 

# In[57]:


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

wilcoxon_key_list=[]
for key in epith_har2.uns.keys():
    if key.startswith('wilcoxon_filtered_louvain_'):
        wilcoxon_key_list.append(key)

for i in wilcoxon_key_list:
    df = rank_genes_groups_df(epith_har2,i)
    df = df[~df['names'].isnull()] 
    df = df[df['pvals_adj']<0.05]
    
    #only keep the top 100 genes for each cluster. Genes were sorted by the Score. 
    new_df = pd.DataFrame([], columns = df.columns)
    for index in np.unique(df.index):
        new_df = pd.concat([new_df, df.loc[index].head(n=100)])
    
    new_df.to_csv('../results/Epithelial_louvain_markers/%s_markers.csv'%i)
    
    del df, new_df


# ## Doublets recheck

# In[34]:


sc.pl.umap(epith_har2,color=['Clic6','Folr1','Acta2','Des','Col3a1','Pecam1','Cd34','C1qa','Ptprc','Bcan','Cnp','Deup1'],ncols=4)


# In[35]:


#remove redundant info
wilcoxon_key_list=[]
for key in epith_har2.uns.keys():
    if key.startswith('wilcoxon_filtered_'):
        wilcoxon_key_list.append(key)

for i in wilcoxon_key_list:
    del epith_har2.uns[i]
    
epith_har2.write_h5ad('../results/Epithelial_har_noDbl_2.h5ad')


# # Subcells_Version2
# ## Reference markers
# ### Lehtinen

# In[37]:


ciliogenesis = ['Deup1','Shisa8','Cdc20'] #ciliogenesis
proliferating = ['Hnrnpa2b1','Tcp1','Rsph1']# Proliferating
IEG = ['Fos','Dusp1','Jun']# IEG expression
sc.pl.umap(epith_har2,color = ciliogenesis+proliferating+IEG,wspace=0.6,ncols=3)


# ## Markers from Laura organoids
# ![image.png](attachment:3c1c7feb-c808-4ce9-8080-430cabd935b5.png)

# In[40]:


sc.pl.umap(epith_har2,color = ['Ppargc1a','Htr2c','Cldn5','Foxj1','Kif3a','Tppp3','Acta2','Myl9','Mki67','Top2a','Cenpe','Ube2c'])


# ## Assign subcells

# In[79]:


epith_har2 = sc.read_h5ad('../results/Epithelial_har_noDbl_2.h5ad')


# In[80]:


epith_har2.obs['Epithelial_subcell']=epith_har2.obs['louvain_0.2']
epith_har2.obs['Epithelial_subcell'] = epith_har2.obs['Epithelial_subcell'].cat.add_categories(['Ciliogenesis','Cidea+'])
epith_har2.obs['Epithelial_subcell'].replace({'0':'EP','1':'EP','2':'EP Proliferating','3':'Gata3+'},inplace=True)
for cell_id in epith_har2.obs.index[(epith_har2.obs['louvain_0.7']=='7') & (epith_har2.obs['Epithelial_subcell']=='EP Proliferating')]:
    epith_har2.obs.loc[cell_id,'Epithelial_subcell']='Ciliogenesis'
for cell_id in epith_har2.obs.index[(epith_har2.obs['louvain_0.8']=='7')]:
    epith_har2.obs.loc[cell_id,'Epithelial_subcell']='Cidea+'

sc.pl.umap(epith_har2,color='Epithelial_subcell')

    


# ### Top markers

# In[81]:


sc.pl.umap(epith_har2,color=['Cidea'])


# ### Top markers

# In[82]:


sc.tl.rank_genes_groups(epith_har2,'Epithelial_subcell',method='wilcoxon',key_added='wilcoxon_Epithelial_subcell')
sc.tl.dendrogram(epith_har2,groupby='Epithelial_subcell')
sc.tl.filter_rank_genes_groups(epith_har2, groupby='Epithelial_subcell', key = "wilcoxon_Epithelial_subcell", key_added='wilcoxon_filtered_Epithelial_subcell',
                               min_in_group_fraction=0.1,max_out_group_fraction=0.5,min_fold_change=0.25)
sc.pl.rank_genes_groups_dotplot(epith_har2,n_genes=15,key='wilcoxon_filtered_Epithelial_subcell', groupby='Epithelial_subcell',show=True)


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

df = rank_genes_groups_df(epith_har2,'wilcoxon_filtered_Epithelial_subcell')
df = df[df['pvals_adj']<0.05]
df = df[~df['names'].isnull()]

#only keep the top 100 genes for each cluster. Genes were sorted by the Score. 
new_df = pd.DataFrame([], columns = df.columns)
for index in np.unique(df.index):
    new_df = pd.concat([new_df, df.loc[index].head(n=200)])
    new_df.to_csv('Epithelial_SubcellType_markers_filtered.csv')


# ## Save the data

# In[83]:


#remove redundant info
wilcoxon_key_list=[]
for key in epith_har2.uns.keys():
    if key.startswith('wilcoxon_filtered_'):
        wilcoxon_key_list.append(key)

for i in wilcoxon_key_list:
    del epith_har2.uns[i]
    
epith_har2.write_h5ad('../results/Epithelial_har_noDbl_2.h5ad')


# # Subcells_Version3
# Regress out IEG genes

# In[60]:


epith_har2 = sc.read_h5ad('../results/Epithelial_har_noDbl_2.h5ad')
IEGs = pd.read_csv('../../../IEGs/IEGs_mm.csv',delimiter=' ')
IEGs_mm = []
for value in IEGs.x:
    if value in epith_har2.var_names:
        IEGs_mm.append(value)
# estimate the IEGs score
epith_har2_raw =  adata_raw[epith_har2.obs.index]    
pct_IEGs = np.sum(epith_har2_raw[:,IEGs_mm].X,axis=1).A1/np.sum(epith_har2_raw.X,axis=1).A1*100
epith_har2.obs['pct_IEGs'] = pct_IEGs
del epith_har2_raw
gc.collect()

epith3 = epith_har2.raw.to_adata()
epith3.raw = epith3
epith3 = epith3[:,epith_har2.var_names]

sc.pp.regress_out(epith3,['pct_IEGs'])
sc.pp.scale(epith3,max_value=10)

sc.tl.pca(epith3,svd_solver = 'arpack',n_comps = 50)
sc.pl.pca_variance_ratio(epith3,n_pcs=50)
#umap
sc.pp.neighbors(epith3,n_pcs = 10, n_neighbors=15)
sc.tl.umap(epith3)
#sc.tl.tsne(epith3,n_pcs=20,random_state=0)
epith3.write_h5ad('../results/Epithelial_withBE_noDbl_3.h5ad')

#harmony
import harmonypy
import scanpy.external as sce
epith_har3 = epith3
sce.pp.harmony_integrate(epith_har3, 'batch',max_iter_harmony=20)
# tsne and umap
sc.pp.neighbors(epith_har3, n_pcs =50, use_rep = "X_pca_harmony")
sc.tl.umap(epith_har3)
sc.tl.tsne(epith_har3, n_pcs = 50, use_rep = "X_pca_harmony")
epith_har3.write_h5ad('../results/Epithelial_har_noDbl_3.h5ad')

#remove redundant info
wilcoxon_key_list=[]
for key in epith_har3.uns.keys():
    if key.startswith('wilcoxon_louvain_'):
        wilcoxon_key_list.append(key)

for i in wilcoxon_key_list:
    del epith_har3.uns[i]

for name in epith_har3.obs.columns:
    if name.startswith('louvain'):
        del epith_har3.obs[name]
        
res = [0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.0]
for snn_res in res:
    sc.tl.louvain(epith_har3, resolution = snn_res, key_added = 'louvain_%s' %snn_res) 

#Top markers of louvain snn
res_list = epith_har3.obs.columns[epith_har3.obs.columns.str.startswith('louvain')]
for snn_res in res_list:
    sc.tl.rank_genes_groups(epith_har3,snn_res,method='wilcoxon',key_added='wilcoxon_%s' %snn_res)
    sc.tl.dendrogram(epith_har3,groupby=snn_res)
    sc.tl.filter_rank_genes_groups(epith_har3, groupby=snn_res, key = "wilcoxon_%s" %snn_res, key_added='wilcoxon_filtered_%s' %snn_res,
                                   min_in_group_fraction=0.1,max_out_group_fraction=0.5,min_fold_change=0.25)
sc.pl.rank_genes_groups_dotplot(epith_har3,n_genes=10,key='wilcoxon_filtered_louvain_0.2' , groupby='louvain_0.2',show=True)
    
#plot
sc.pl.umap(epith_har3,color=res_list,ncols=4, wspace=0.3)


# ## Assign subcells

# In[62]:


for snn_res in res_list:
    sc.pl.rank_genes_groups_dotplot(epith_har3,n_genes=10,key='wilcoxon_filtered_%s' %snn_res, groupby=snn_res,show=True)


# In[64]:


sc.pl.umap(epith_har3,color = ['Ppargc1a','Htr2c','Deup1','Foxj1','Kif3a','Tppp3','Acta2','Myl9','Mki67','Top2a','Cenpe','Ube2c'])


# In[84]:


Epithelial_subcell = epith_har2.obs.loc[epith_har3.obs.index,'Epithelial_subcell']
epith_har3.obs['Epithelial_subcell']  =Epithelial_subcell
epith_har3.obs['Epithelial_subcell_V3']=epith_har3.obs['louvain_0.4']
epith_har3.obs['Epithelial_subcell_V3'].replace({'0':'EP','1':'EP','2':'EP','3':'Gata3+',
                                                 '4':'EP Proliferating','5':'Cidea+','6':'Neuro?'},inplace=True)
for cell_id in epith_har3.obs.index[epith_har3.obs['louvain_0.2']=='3']:
    epith_har3.obs.loc[cell_id,'Epithelial_subcell_V3'] ='Ciliogenesis'
    #if cell_id in epith_har3.obs.index[epith_har3.obs['Epithelial_subcell_V3']=='Ciliogenesis']:
    #    epith_har3.obs.loc[cell_id,'Epithelial_subcell_V3'] ='EP Proliferating'


# In[85]:


sc.pl.umap(epith_har3,color=['Epithelial_subcell','Epithelial_subcell_V3'],wspace=0.7,size=10)


# !!!! Regress out IEGs doesn't cause a very different clusterings results. So I think we don't need to regress out IEGs during subcell clustering.
# 
# SO I will keep IEGs in whole process.

# In[86]:


#remove redundant info
wilcoxon_key_list=[]
for key in epith_har3.uns.keys():
    if key.startswith('wilcoxon_filtered_'):
        wilcoxon_key_list.append(key)

for i in wilcoxon_key_list:
    del epith_har3.uns[i]
    
epith_har3.write_h5ad('../results/Epithelial_har_noDbl_3.h5ad')


# # Final version: Remove neuro-epith doublets in subcell_Version2

# In[87]:


# find cells ID after remove neuros
epith_har4 = epith_har3[epith_har3.obs['Epithelial_subcell_V3']!='Neuro?']
epith_har2_edit = epith_har2[epith_har4.obs.index]


# ## Top markers

# In[88]:


sc.tl.rank_genes_groups(epith_har2_edit,'Epithelial_subcell',method='wilcoxon',key_added='wilcoxon_Epithelial_subcell')
sc.tl.dendrogram(epith_har2_edit,groupby='Epithelial_subcell')
sc.tl.filter_rank_genes_groups(epith_har2_edit, groupby='Epithelial_subcell', key = "wilcoxon_Epithelial_subcell", key_added='wilcoxon_filtered_Epithelial_subcell',
                               min_in_group_fraction=0.1,max_out_group_fraction=0.5,min_fold_change=0.25)
sc.pl.rank_genes_groups_dotplot(epith_har2_edit,n_genes=15,key='wilcoxon_filtered_Epithelial_subcell', groupby='Epithelial_subcell',show=True)


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

df = rank_genes_groups_df(epith_har2_edit,'wilcoxon_filtered_Epithelial_subcell')
df = df[df['pvals_adj']<0.05]
df = df[~df['names'].isnull()]

#only keep the top 100 genes for each cluster. Genes were sorted by the Score. 
new_df = pd.DataFrame([], columns = df.columns)
for index in np.unique(df.index):
    new_df = pd.concat([new_df, df.loc[index].head(n=200)])
    new_df.to_csv('Epithelial_SubcellType_edited_markers_filtered.csv')


# ## Cell summary

# In[91]:


fig, (ax1,ax2) = plt.subplots(nrows=1,ncols=2,constrained_layout=True,figsize=(12,4))
sns.countplot(y = 'batch',data = epith_har2_edit.obs,ax=ax1)
for container in ax1.containers:
    ax1.bar_label(container)
    
sns.countplot(y = 'Age',data = epith_har2_edit.obs,ax=ax2)
for container in ax2.containers:
    ax2.bar_label(container)

fig, (ax1) = plt.subplots(nrows=1,ncols=1,constrained_layout=True,figsize=(8,4))
sns.countplot(y = 'Epithelial_subcell',data = epith_har2_edit.obs,ax=ax1)
for container in ax1.containers:
    ax1.bar_label(container)

def batch_pct_cellType(adata_file,groupby_key,legend_key,plot_title,pos):
    tmp = pd.crosstab(adata_file.obs[groupby_key],adata_file.obs[legend_key], normalize='index')   
    tmp.plot.bar(stacked=True,figsize=(12,10),title=plot_title, rot= 0, grid=False,
                 xlabel='',ax=pos).legend(bbox_to_anchor=[1.15,1])
    
fig, (ax1,ax2,ax3) = plt.subplots(3, 1, constrained_layout=True)
batch_pct_cellType(epith_har2_edit,'Epithelial_subcell','batch','',ax1)
batch_pct_cellType(epith_har2_edit,'Epithelial_subcell','Age','',ax2)
batch_pct_cellType(epith_har2_edit,'Age','Epithelial_subcell','',ax3)


# ## Go analysis

# In[92]:


def Enrichr_Plot(geneset_name,enr_res,clt):
    df = enr_res.results[enr_res.results['Gene_set']==geneset_name]
    df_sig = df[df['Adjusted P-value']<0.05]
    x_values = df_sig['Adjusted P-value'].apply(lambda x: -np.log10(x)).to_list()
    y_values = df_sig['Term']
    
    if df_sig.shape[0]>20:
        ax = sns.barplot(x = x_values[:20], y =y_values[:20],color='#78c679')
        ax.set_title('Enrichment from %s were found in %s' %(geneset_name,clt))
        ax.set_xlabel('-log10(Adjusted P-value)')
    else:
        ax = sns.barplot(x = x_values, y =y_values,color='#78c679')
        ax.set_title('Enrichment from %s were found in %s' %(geneset_name,clt))
        ax.set_xlabel('-log10(Adjusted P-value)')


# In[93]:


gene_set_names = gseapy.get_library_name(database='Mouse')


# ### EP Proliferating

# In[94]:


glist_pro = sc.get.rank_genes_groups_df(epith_har2_edit, group='EP Proliferating', key='wilcoxon_filtered_Epithelial_subcell')[['names','logfoldchanges','pvals_adj','scores']]
glist_pro = glist_pro[~glist_pro['names'].isnull()]
enr_res_pro = gseapy.enrichr(gene_list=glist_pro['names'][glist_pro['pvals_adj']<0.05].tolist(),
                         organism='mouse',
                         gene_sets=['GO_Biological_Process_2021','GO_Molecular_Function_2021'],
                         description='pathway',
                         cutoff = 0.5,
                         outdir ='../GoTerm_subcell/Enrichr/Epithelial/Subcell/EP_Pro')
# check if enrichment were found in all geneset
for geneset_name in ['GO_Biological_Process_2021','GO_Molecular_Function_2021']:
    df = enr_res_pro.results[enr_res_pro.results['Gene_set']==geneset_name]
    df_sig = df[df['Adjusted P-value']<0.05]
    
    if df_sig.shape[0]==0:
        print('No enrichments from %s!' %geneset_name)
        
    else:
        print('Oh YEAH there are enrichments in %s!' %geneset_name)
    
enr_res_pro_sig = enr_res_pro.results[enr_res_pro.results['Adjusted P-value']<0.05]
print(enr_res_pro_sig.shape)
enr_res_pro_sig.to_csv('../GoTerm_subcell/Enrichr/Epithelial/Subcell/enr_res_EP_pro_sig.csv')


# In[95]:


# plot geneset with enrichments
for geneset_name in ['GO_Biological_Process_2021']:
    fig, ax = plt.subplots(1,1, figsize=(10, 8))
    Enrichr_Plot(geneset_name,enr_res_pro,'EP Proliferating')


# ### Ciliogenesis

# In[96]:


glist_cilio = sc.get.rank_genes_groups_df(epith_har2_edit, group='Ciliogenesis', key='wilcoxon_filtered_Epithelial_subcell')[['names','logfoldchanges','pvals_adj','scores']]
glist_cilio = glist_cilio[~glist_cilio['names'].isnull()]
enr_res_cilio = gseapy.enrichr(gene_list=glist_cilio['names'][glist_cilio['pvals_adj']<0.05].tolist(),
                         organism='mouse',
                         gene_sets=['GO_Biological_Process_2021','GO_Molecular_Function_2021'],
                         description='pathway',
                         cutoff = 0.5,
                         outdir ='../GoTerm_subcell/Enrichr/Epithelial/Subcell/Ciliogenesis')
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
enr_res_cilio_sig.to_csv('../GoTerm_subcell/Enrichr/Epithelial/Subcell/enr_res_Ciliogenesis_sig.csv')


# In[98]:


#plot geneset with enrichments
for geneset_name in ['GO_Biological_Process_2021']:
    fig, ax = plt.subplots(1,1, figsize=(10, 8))
    Enrichr_Plot(geneset_name,enr_res_cilio,'Ciliogenesis')


# ### Gata3

# In[99]:


glist_gata3 = sc.get.rank_genes_groups_df(epith_har2_edit, group='Gata3+', key='wilcoxon_filtered_Epithelial_subcell')[['names','logfoldchanges','pvals_adj','scores']]
glist_gata3 = glist_gata3[~glist_gata3['names'].isnull()]
enr_res_gata3 = gseapy.enrichr(gene_list=glist_gata3['names'][glist_gata3['pvals_adj']<0.05].tolist(),
                         organism='mouse',
                         gene_sets=['GO_Biological_Process_2021','GO_Molecular_Function_2021'],
                         description='pathway',
                         cutoff = 0.5,
                         outdir ='../GoTerm_subcell/Enrichr/Epithelial/Subcell/Gata3+')
# check if enrichment were found in all geneset
for geneset_name in ['GO_Biological_Process_2021','GO_Molecular_Function_2021']:
    df = enr_res_gata3.results[enr_res_gata3.results['Gene_set']==geneset_name]
    df_sig = df[df['Adjusted P-value']<0.05]
    
    if df_sig.shape[0]==0:
        print('No enrichments from %s!' %geneset_name)
        
    else:
        print('Oh YEAH there are enrichments in %s!' %geneset_name)
    
enr_res_gata3_sig = enr_res_gata3.results[enr_res_gata3.results['Adjusted P-value']<0.05]
print(enr_res_gata3_sig.shape)
enr_res_gata3_sig.to_csv('../GoTerm_subcell/Enrichr/Epithelial/Subcell/enr_res_Gata3_sig.csv')


# In[100]:


# plot geneset with enrichments
for geneset_name in ['GO_Biological_Process_2021','GO_Molecular_Function_2021']:
    fig, ax = plt.subplots(1,1, figsize=(10, 8))
    Enrichr_Plot(geneset_name,enr_res_gata3,'Gata3+')


# ### Cidea

# In[102]:


glist_cidea = sc.get.rank_genes_groups_df(epith_har2_edit, group='Cidea+', key='wilcoxon_filtered_Epithelial_subcell')[['names','logfoldchanges','pvals_adj','scores']]
glist_cidea = glist_cidea[~glist_cidea['names'].isnull()]
enr_res_cidea = gseapy.enrichr(gene_list=glist_cidea['names'][glist_cidea['pvals_adj']<0.05].tolist(),
                         organism='mouse',
                         gene_sets=['GO_Biological_Process_2021','GO_Molecular_Function_2021'],
                         description='pathway',
                         cutoff = 0.5,
                         outdir ='../GoTerm_subcell/Enrichr/Epithelial/Subcell/cidea')
# check if enrichment were found in all geneset
for geneset_name in ['GO_Biological_Process_2021','GO_Molecular_Function_2021']:
    df = enr_res_cidea.results[enr_res_cidea.results['Gene_set']==geneset_name]
    df_sig = df[df['Adjusted P-value']<0.05]
    
    if df_sig.shape[0]==0:
        print('No enrichments from %s!' %geneset_name)
        
    else:
        print('Oh YEAH there are enrichments in %s!' %geneset_name)
    
enr_res_cidea_sig = enr_res_cidea.results[enr_res_cidea.results['Adjusted P-value']<0.05]
print(enr_res_cidea_sig.shape)
enr_res_cidea_sig.to_csv('../GoTerm_subcell/Enrichr/Epithelial/Subcell/enr_res_cidea_sig.csv')


# ### EP

# In[103]:


glist_EP = sc.get.rank_genes_groups_df(epith_har2_edit, group='EP', key='wilcoxon_filtered_Epithelial_subcell')[['names','logfoldchanges','pvals_adj','scores']]
glist_EP = glist_EP[~glist_EP['names'].isnull()]
enr_res_EP = gseapy.enrichr(gene_list=glist_EP['names'][glist_EP['pvals_adj']<0.05].tolist(),
                         organism='mouse',
                         gene_sets=['GO_Biological_Process_2021','GO_Molecular_Function_2021'],
                         description='pathway',
                         cutoff = 0.5,
                         outdir ='../GoTerm_subcell/Enrichr/Epithelial/Subcell/EP2')
# check if enrichment were found in all geneset
for geneset_name in ['GO_Biological_Process_2021','GO_Molecular_Function_2021']:
    df = enr_res_EP.results[enr_res_EP.results['Gene_set']==geneset_name]
    df_sig = df[df['Adjusted P-value']<0.05]
    
    if df_sig.shape[0]==0:
        print('No enrichments from %s!' %geneset_name)
        
    else:
        print('Oh YEAH there are enrichments in %s!' %geneset_name)
    
enr_res_EP_sig = enr_res_EP.results[enr_res_EP.results['Adjusted P-value']<0.05]
print(enr_res_EP_sig.shape)
enr_res_EP_sig.to_csv('../GoTerm_subcell/Enrichr/Epithelial/Subcell/enr_res_EP_sig.csv')


# In[104]:


# plot geneset with enrichments
for geneset_name in ['GO_Biological_Process_2021']:
    fig, ax = plt.subplots(1,1, figsize=(10, 8))
    Enrichr_Plot(geneset_name,enr_res_EP,'EP')


# ## Save the data

# In[111]:


#remove redundant info
wilcoxon_key_list=[]
for key in epith_har2_edit.uns.keys():
    if key.startswith('wilcoxon_filtered_'):
        wilcoxon_key_list.append(key)

for i in wilcoxon_key_list:
    del epith_har2_edit.uns[i]
    
epith_har2_edit.write_h5ad('../results/Epithelial_har_noDbl_2_edited.h5ad')


# # Correct CellType for file including all cells

# In[105]:


epith = sc.read_h5ad('../results/Epithelial_har.h5ad')
adata_har = sc.read_h5ad('../../results/Mouse_harmony_rmBatchbetweenSamples.h5ad')


# In[106]:


new_doublets = []
for cell_id in epith.obs.index:
    if cell_id not in epith_har2_edit.obs.index:
        new_doublets.append(cell_id)


# In[107]:


adata_har.obs['CellType_edited4']=adata_har.obs['CellType_edited3']
for cell_id in new_doublets:
    adata_har.obs.loc[cell_id,'CellType_edited4']='Doublets'


# In[108]:


sc.pl.umap(adata_har,color='CellType_edited4')


# In[109]:


adata_har.write_h5ad('../../results/Mouse_harmony_rmBatchbetweenSamples.h5ad')


# In[110]:


adata_har_noDbl = adata_har[adata_har.obs['CellType_edited4']!='Doublets']
adata_har_noDbl.write_h5ad('../../results/Mouse_harmony_rmBatchbetweenSamples_rmDoublets.h5ad')


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




