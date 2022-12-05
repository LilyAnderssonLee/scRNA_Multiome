#!/usr/bin/env python
# coding: utf-8

# scDRS analysis for human sampels was ran under conda env: _/proj/snic2021-23-715/private/Lili/analysis/Multiome_scRNA_ATAC/multiome_env_
# 
# Single cell Disease Related Score tutorial: 
# https://github.com/martinjzhang/scDRS
# 

# In[1]:


import scdrs
import scanpy as sc
from anndata import AnnData
from scipy import stats
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import warnings


sc.settings.verbosity = 3             
sc.settings.set_figure_params(dpi=150,figsize=(4,4))


# # Human

# In[2]:


human = sc.read_h5ad('../../results/old/Human_harmony_CellType_edit_noDBL.h5ad')


# **Remove doublets**

# In[3]:


human_noDBL = human[(human.obs['Immune_subcell_Ana']!='Mac-Epithelial') & (human.obs['Immune_subcell_Ana']!='Mac-Mesenchymal')]


# In[4]:


sc.pl.umap(human_noDBL,color=['Immune_subcell_Ana','Epithelial_subcell_Ana','Mesenchymal_subcell_Ana','Endothelial_Subcells',
                              'NeuroGlia_subcell_Ana','CellType_intergrate_edited5'],wspace=0.7,ncols=3)


# ## PP h5ad file
# 
# **scDRS use the raw counts**

# **_Combine all subcell info_**

# In[5]:


human_noDBL.obs['CellType_subcell'] = human_noDBL.obs.CellType_intergrate_edited5


# In[6]:


human_noDBL.obs['CellType_subcell'] = human_noDBL.obs['CellType_subcell'].cat.add_categories(['Pericytes','Fibroblast','Myoepithelial','Neurons','OPC',
                                                                                              'Oligodendrocytes','Astrocyte','T cells','NK cells','B cells','Plasmocytes','Dendritic',
                                                                                             'Macrophage'])


# In[7]:


human_noDBL.obs['CellType_subcell'].replace({'Mesenchymal':'Fibroblast'},inplace=True)


# In[8]:


for cell_id in human_noDBL.obs.index[(human_noDBL.obs.Mesenchymal_subcell_Ana=='vSMC') | (human_noDBL.obs.Mesenchymal_subcell_Ana=='Pericytes')]:
    human_noDBL.obs.loc[cell_id,'CellType_subcell'] = 'Pericytes'
for cell_id in human_noDBL.obs.index[human_noDBL.obs.Epithelial_subcell_Ana=='Myoepithelial']:
    human_noDBL.obs.loc[cell_id,'CellType_subcell'] = 'Myoepithelial'


# **For Ependymal cells, I will only treat ones seperated from Astrocytes cells.**

# In[9]:


#for cell_id in human_noDBL.obs.index[human_noDBL.obs.NeuroGlia_subcell_Ana=='Astrocyte']:
#    human_noDBL.obs.loc[cell_id,'CellType_subcell'] = 'Astrocyte'

for cell_id in human_noDBL.obs.index[human_noDBL.obs.NeuroGlia_subcell_Ana=='OPC']:
    human_noDBL.obs.loc[cell_id,'CellType_subcell'] = 'OPC'
for cell_id in human_noDBL.obs.index[human_noDBL.obs.NeuroGlia_subcell_Ana=='Oligodendrocytes']:
    human_noDBL.obs.loc[cell_id,'CellType_subcell'] = 'Oligodendrocytes'

for cell_id in human_noDBL.obs.index[(human_noDBL.obs.NeuroGlia_subcell_Ana=='Neurons 1') | (human_noDBL.obs.NeuroGlia_subcell_Ana=='Neurons 2')]:
    human_noDBL.obs.loc[cell_id,'CellType_subcell'] = 'Neurons'


# In[10]:


human_noDBL.obs['CellType_subcell'].replace({'NeuronGlial':'Astrocytes'},inplace=True)


# In[11]:


for cell_id in human_noDBL.obs.index[human_noDBL.obs.Immune_subcell_Ana=='NK cells']:
    human_noDBL.obs.loc[cell_id,'CellType_subcell'] = 'NK cells'
for cell_id in human_noDBL.obs.index[human_noDBL.obs.Immune_subcell_Ana=='T cells']:
    human_noDBL.obs.loc[cell_id,'CellType_subcell'] = 'T cells'
for cell_id in human_noDBL.obs.index[human_noDBL.obs.Immune_subcell_Ana=='B cells']:
    human_noDBL.obs.loc[cell_id,'CellType_subcell'] = 'B cells'
for cell_id in human_noDBL.obs.index[human_noDBL.obs.Immune_subcell_Ana=='T cells']:
    human_noDBL.obs.loc[cell_id,'CellType_subcell'] = 'T cells'
for cell_id in human_noDBL.obs.index[human_noDBL.obs.Immune_subcell_Ana=='Plasmocytes']:
    human_noDBL.obs.loc[cell_id,'CellType_subcell'] = 'Plasmocytes'
for cell_id in human_noDBL.obs.index[(human_noDBL.obs.Immune_subcell_Ana=='Dendritic') | (human_noDBL.obs.Immune_subcell_Ana=='Dendritic cells ITGAX')]:
    human_noDBL.obs.loc[cell_id,'CellType_subcell'] = 'Dendritic'
for cell_id in human_noDBL.obs.index[(human_noDBL.obs.Immune_subcell_Ana=='Mac Proliferation') | (human_noDBL.obs.Immune_subcell_Ana=='Mac') | 
                                    (human_noDBL.obs.Immune_subcell_Ana=='Phagocyting Mac') | (human_noDBL.obs.Immune_subcell_Ana=='Migratory Mac') |
                                    (human_noDBL.obs.Immune_subcell_Ana=='Activated Mac: chemotaxis high')]:
    human_noDBL.obs.loc[cell_id,'CellType_subcell'] = 'Macrophage'


# In[12]:


human_noDBL.obs.CellType_subcell = human_noDBL.obs.CellType_subcell.cat.remove_categories(['Astrocyte','Immune'])


# In[13]:


sc.pl.umap(human_noDBL,color=['CellType_intergrate_edited5','CellType_subcell'],wspace=0.7)


# **Load raw counts**

# In[15]:


human_raw = sc.read_h5ad('../../results/old/Human_filtered_rawCounts.h5ad')


# In[16]:


human_noDBL_raw = human_raw[human_noDBL.obs.index]
human_noDBL_raw.obs = human_noDBL.obs
human_noDBL_raw.obsm = human_noDBL.obsm
human_noDBL_raw.obsp = human_noDBL.obsp


# In[18]:


human_noDBL_raw.write_h5ad('../scDRS/Ana_human_ms/human_noDBL_scDRS.h5ad')


# In[20]:


human_noDBL_raw = sc.read_h5ad('../scDRS/Ana_human_ms/human_noDBL_scDRS.h5ad')


# ## PP cov
# 
#     scDRS covariate file for the .h5ad single-cell data. .tsv file.
# 
#     First column: cell names, consistent with adata.obs_names.
# 
#     Other comlumns: covariates with numerical values.
# 

# In[19]:


human_cov_dict = {'n_genes':human_noDBL_raw.obs.n_genes}


# In[20]:


human_cov=pd.DataFrame.from_dict(human_cov_dict,orient='columns')


# In[105]:


#human_cov.rename(columns={'CellType_intergrate_edited5':'level1Class','CellType_subcell':'level2Class'},inplace=True)
#identify all categorical variables
#cat_columns = human_cov.select_dtypes(['category']).columns
#convert all categorical variables to numeric
#human_cov[cat_columns] = human_cov[cat_columns].apply(lambda x: pd.factorize(x)[0])


# In[21]:


human_cov.to_csv('../scDRS/Ana_human_ms/human_noDBL_cov.tsv',sep='\t',header=True)


# ## PP gene set file
# Include the height related geneset from the tutorial, as a control

# In[6]:


ref_gs = pd.read_csv('../scDRS/data/geneset.gs',index_col=0,delimiter='\t')


# **Convert mouse genes into human genes**

# In[7]:


matched_genes = pd.read_csv('../scDRS/scdrs/data/mouse_human_homologs.txt',sep='\t',index_col=0)


# ### Gene set of height

# In[8]:


ref_gs_height= ref_gs[(ref_gs.index=='UKB_460K.body_HEIGHTz')]
tmp_splt = ref_gs_height.GENESET[0].split(',')

gene_name = []
gene_value = []
for item in tmp_splt:
    gene_name.append(item.split(":")[0])
    gene_value.append(item.split(':')[1])
gene_set ={}
for key in gene_name:
    for value in gene_value:
        gene_set[key] = value
        gene_value.remove(value)
        break

gene_set_df_height = pd.DataFrame(gene_set.items(),columns=['Gene','Score'])
gene_set_df_height['human_gene']= gene_set_df_height.Gene
for gene in gene_set_df_height.Gene:
    index_id = gene_set_df_height.index[gene_set_df_height.Gene==gene]
    if gene not in matched_genes.index:
        gene_set_df_height.drop(index_id, axis=0)
    else: 
        gene_set_df_height.loc[index_id,'human_gene'] = matched_genes.HUMAN_GENE_SYM[matched_genes.index==gene].values


# In[9]:


gene_set_df_height.Gene = gene_set_df_height.human_gene
gene_set_df_height.drop('human_gene',axis=1,inplace=True)


# Intersect genes between gene set and human_noDBL_raw.var_names

# In[10]:


for item in gene_set_df_height.Gene:
    if item not in human_noDBL_raw.var_names:
        gene_set_df_height.drop(gene_set_df_height.index[gene_set_df_height.Gene==item],axis=0,inplace=True)


# In[12]:


gene_set_df_height = gene_set_df_height.astype({'Score':'float'})


# ### Gene set of Alzheimers

# In[13]:


ref_gs_Alzheimers= ref_gs[(ref_gs.index=='PASS_Alzheimers_Jansen2019')]
tmp_splt = ref_gs_Alzheimers.GENESET[0].split(',')

gene_name = []
gene_value = []
for item in tmp_splt:
    gene_name.append(item.split(":")[0])
    gene_value.append(item.split(':')[1])
gene_set ={}
for key in gene_name:
    for value in gene_value:
        gene_set[key] = value
        gene_value.remove(value)
        break

gene_set_df_Alzheimers = pd.DataFrame(gene_set.items(),columns=['Gene','Score'])
gene_set_df_Alzheimers['human_gene']= gene_set_df_Alzheimers.Gene
for gene in gene_set_df_Alzheimers.Gene:
    index_id = gene_set_df_Alzheimers.index[gene_set_df_Alzheimers.Gene==gene]
    if gene not in matched_genes.index:
        gene_set_df_Alzheimers.drop(index_id, axis=0)
    else: 
        gene_set_df_Alzheimers.loc[index_id,'human_gene'] = matched_genes.HUMAN_GENE_SYM[matched_genes.index==gene].values


# In[14]:


gene_set_df_Alzheimers.Gene = gene_set_df_Alzheimers.human_gene
gene_set_df_Alzheimers.drop('human_gene',axis=1,inplace=True)


# Intersect genes between gene set and human_noDBL_raw.var_names

# In[ ]:


for item in gene_set_df_Alzheimers.Gene:
    if item not in human_noDBL_raw.var_names:
        gene_set_df_Alzheimers.drop(gene_set_df_Alzheimers.index[gene_set_df_Alzheimers.Gene==item],axis=0,inplace=True)


# In[16]:


gene_set_df_Alzheimers = gene_set_df_Alzheimers.astype({'Score':'float'})


# ## Combine gene sets from height, Alzheimer and MS

# In[4]:


#MS
gene_set= scdrs.util.load_gs('../scDRS/Ana_human_ms/ms_gwas_1_scDRS.gs',to_intersect=human_noDBL_raw.var_names)


# In[17]:


gene_set['Height'] = (gene_set_df_height.Gene.tolist(),gene_set_df_height.Score.tolist())
gene_set['Alzheimers'] = (gene_set_df_Alzheimers.Gene.tolist(),gene_set_df_Alzheimers.Score.tolist())


# ## scDRS analysis of disease enrichment for individual cells

# In[2]:


adata = scdrs.util.load_h5ad('../scDRS/Ana_human_ms/human_noDBL_scDRS.h5ad',flag_filter_data=True,flag_raw_count=True)


# In[3]:


human_cov = pd.read_table('../scDRS/Ana_human_ms/human_noDBL_cov.tsv',index_col=0)


# In[18]:


# preprocess data to
# (1) regress out from covariates
# (2) group genes into bins by mean and variance
scdrs.preprocess(adata, cov=human_cov, n_mean_bin=20, n_var_bin=20, copy=False)


# In[22]:


human_noDBL_raw.uns['SCDRS_PARAM'] = adata.uns['SCDRS_PARAM']


# ## Computing scDRS scores

# In[131]:


dict_df_score = dict()
for trait in gene_set:
    gene_list, gene_weights = gene_set[trait]
    dict_df_score[trait] = scdrs.score_cell(
        data=adata,
        gene_list=gene_list,
        gene_weight=gene_weights,
        ctrl_match_key="mean_var",
        n_ctrl=1000,
        weight_opt="vs",
        return_ctrl_raw_score=False,
        return_ctrl_norm_score=True,
        verbose=False,
    )


# In[146]:


dict_df_score['MS'].to_csv('../scDRS/Ana_human_ms/human_MS_scDRS_score.tsv',header=True,sep='\t')


# In[148]:


dict_df_score['Height'].to_csv('../scDRS/Ana_human_ms/human_Height_scDRS_score.tsv',header=True,sep='\t')


# In[149]:


dict_df_score['Alzheimers'].to_csv('../scDRS/Ana_human_ms/human_Alzheimers_scDRS_score.tsv',header=True,sep='\t')


# In[155]:


human_noDBL_raw.obs['scDRS_MS_norm_score'] = dict_df_score['MS'].norm_score
human_noDBL_raw.obs['scDRS_height_norm_score'] = dict_df_score['Height'].norm_score
human_noDBL_raw.obs['scDRS_Alzheimers_norm_score'] = dict_df_score['Alzheimers'].norm_score


# In[161]:


sc.pl.umap(human_noDBL_raw,color=['scDRS_MS_norm_score','scDRS_height_norm_score','scDRS_Alzheimers_norm_score',
                                 'CellType_subcell'],
           color_map="RdBu_r",vmin=-5,vmax=5,wspace=0.8,ncols=2,size=3)


# **Conclusion**
# 
# The immune cells have enriched MS related gene signals. So I will zoom in Immune cells to check the geneset enrichment in each subcells.
# 
# No cell types have strong height signals which is reasonable since the height is controled by low allele frequency shift. 
# 
# Macropage and dendritic cells have enriched signals of Alzheimer.

# In[162]:


human_noDBL_raw.write_h5ad('../scDRS/Ana_human_ms/human_noDBL_scDRS_normScore.h5ad')


# ## Downstream analyses
# 
# ### MS

# In[25]:


dict_df_score_MS = pd.read_csv('../scDRS/Ana_human_ms/human_MS_scDRS_score.tsv',index_col=0,delimiter="\t")


# In[27]:


df_stats_ms = scdrs.method.downstream_group_analysis(
    adata=adata,
    df_full_score=dict_df_score_MS,
    group_cols=["CellType_subcell"],
)["CellType_subcell"]


# Column names in the table df_stats_ms:
# 
# **_assoc_mcp:_**
# 
#     significance of cell type-disease association
#     
# **_hetero_mcp:_** 
#     
#     significance heterogeneity in association with disease across individual cells within a given cell type
#     
# **_n_fdr_0.05:_** 
# 
#     number of significantly associated cells (with FDR<0.05 across all cells for a given disease)

# In[28]:


display(df_stats_ms.style.set_caption("Group-level statistics for MS"))


# In[49]:


df_stats_ms.to_csv('../scDRS/Ana_human_ms/human_MS_scDRS_stats.tsv',header=True,sep='\t')


# ### Alzheimers 

# In[52]:


dict_df_score_alzheimers = pd.read_csv('../scDRS/Ana_human_ms/human_Alzheimers_scDRS_score.tsv',index_col=0,delimiter="\t")
df_stats_alzheimers = scdrs.method.downstream_group_analysis(
    adata=adata,
    df_full_score=dict_df_score_alzheimers,
    group_cols=["CellType_subcell"],
)["CellType_subcell"]


# In[53]:


display(df_stats_alzheimers.style.set_caption("Group-level statistics for Alzheimers"))


# In[54]:


df_stats_alzheimers.to_csv('../scDRS/Ana_human_ms/human_Alzheimers_scDRS_stats.tsv',header=True,sep='\t')


# ### Height

# In[55]:


dict_df_score_height = pd.read_csv('../scDRS/Ana_human_ms/human_Height_scDRS_score.tsv',index_col=0,delimiter="\t")
df_stats_height = scdrs.method.downstream_group_analysis(
    adata=adata,
    df_full_score=dict_df_score_height,
    group_cols=["CellType_subcell"],
)["CellType_subcell"]


# In[56]:


display(df_stats_height.style.set_caption("Group-level statistics for Alzheimers"))


# In[57]:


df_stats_height.to_csv('../scDRS/Ana_human_ms/human_Height_scDRS_stats.tsv',header=True,sep='\t')


# Combine stats data from three differnt traits.

# In[59]:


dict_df_stats = {'MS':df_stats_ms,
                'Alzheimers': df_stats_alzheimers,
                'Height': df_stats_height}


# In[61]:


scdrs.util.plot_group_stats(dict_df_stats,df_fdr_prop=0.05)


# **Conclusion**
# 
# we observed that Macrophage and Dendritic cells show the strongest cell type association as well as significant heterogeneity.
# 
# Now we focus on this subset of cells and further understand sources of heterogeneity.

# In[70]:


adata = sc.read_h5ad('../scDRS/Ana_human_ms/human_noDBL_scDRS_normScore.h5ad')


# In[73]:


# extract immune cells and perform a re-clustering
adata_immune = adata[adata.obs.CellType_intergrate_edited5.isin(["Immune"])].copy()
sc.pp.filter_cells(adata_immune, min_genes=0)
sc.pp.filter_genes(adata_immune, min_cells=1)
sc.pp.normalize_total(adata_immune, target_sum=1e4)
sc.pp.log1p(adata_immune)

sc.pp.highly_variable_genes(adata_immune, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata_immune = adata_immune[:, adata_immune.var.highly_variable]
sc.pp.scale(adata_immune, max_value=10)
sc.tl.pca(adata_immune, svd_solver="arpack")

sc.pp.neighbors(adata_immune, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_immune, n_components=2)


# In[88]:


sc.pl.umap(adata_immune,color=['scDRS_MS_norm_score','scDRS_height_norm_score','scDRS_Alzheimers_norm_score',
                               'Condition','CellType_subcell','Immune_subcell_Ana',],
           color_map="RdBu_r",vmin=-5,vmax=5,s=20,ncols=3,wspace=0.5)


# In[89]:


adata_immune.write_h5ad('../scDRS/Ana_human_ms/Human_scDRS_Immune.h5ad')


# **_Attention_**
# 
# You will find the batch effect between samples here since I didn't remove it.

# In[85]:


adata_immune.obs.CellType_subcell.value_counts()


# **Plot the average MS disease score**

# In[ ]:


df_plot = adata.obs[['scDRS_MS_norm_score','scDRS_height_norm_score','scDRS_Alzheimers_norm_score']].copy()
df_plot["MS quintile"] = pd.qcut(df_plot["scDRS_MS_norm_score"], 5, labels=np.arange(5))

fig, ax = plt.subplots(figsize=(3.5, 3.5))
for trait in ["MS", "Height"]:
    sns.lineplot(
        data=df_plot,
        x="Dorsal quintile",
        y=trait,
        label=trait,
        err_style="bars",
        marker="o",
        ax=ax,
    )
ax.set_xticks(np.arange(5))
ax.set_xlabel("Dorsal quintile")
ax.set_ylabel("Mean scDRS disease score")
fig.show()


# In[67]:


adata.obs['scDRS_MS_norm'] = dict_df_score_MS.norm_score
adata.obs['scDRS_Alzheimer_norm'] = dict_df_score_alzheimers.norm_score
adata.obs['scDRS_Height_norm'] = dict_df_score_height.norm_score


# In[69]:


human_noDBL_raw.obs['scDRS_MS_norm'] = dict_df_score_MS.norm_score
human_noDBL_raw.obs['scDRS_Alzheimer_norm'] = dict_df_score_alzheimers.norm_score
human_noDBL_raw.obs['scDRS_Height_norm'] = dict_df_score_height.norm_score


# In[72]:


adata


# In[166]:


import csv
w = csv.writer(open("../scDRS/Ana_human_ms/SCDRS_PARAM.csv", "w"))
# loop over dictionary keys and values
for key, val in adata.uns['SCDRS_PARAM'].items():
    # write every key and value to file
    w.writerow([key, val])


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




