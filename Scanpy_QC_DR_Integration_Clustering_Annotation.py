#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#load packages
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import gc
import os
import platform
sc.settings.verbosity = 3             
sc.settings.set_figure_params(dpi=100,figsize=(5,5))


# In[ ]:


print(platform.platform())


# In[ ]:


pip list


# # Load raw counts

# In[ ]:


adata_dict = {}

adata_dict['Spain_18_006_Ctrl_M_Age51'] = sc.read_10x_h5('../../../../raw_counts_scRNA/xxx/outs/filtered_feature_bc_matrix.h5')
adata_dict['R100_63_Ctrl_F_Age56'] = sc.read_10x_h5('../../../../raw_counts_scRNA/xxx/outs/filtered_feature_bc_matrix.h5')
adata_dict['Spain_19_005_Ctrl_M_Age60'] = sc.read_10x_h5('../../../../raw_counts_scRNA/xxx/outs/filtered_feature_bc_matrix.h5')
adata_dict['R100_64_Ctrl_M_Age64'] = sc.read_10x_h5('../../../../raw_counts_scRNA/xxx/outs/filtered_feature_bc_matrix.h5')
adata_dict['R77_61_Ctrl_F_Age79'] = sc.read_10x_h5('../../../../raw_counts_scRNA/xxx/outs/filtered_feature_bc_matrix.h5')
adata_dict['R100_72_Ctrl_F_Age82'] = sc.read_10x_h5('../../../../raw_counts_scRNA/xxx/outs/filtered_feature_bc_matrix.h5')

adata_dict['R77_60_MS_M_Age46'] = sc.read_10x_h5('../../../../raw_counts_scRNA/xxx/outs/filtered_feature_bc_matrix.h5')
adata_dict['R100_71_MS_F_Age52'] = sc.read_10x_h5('../../../../raw_counts_scRNA/xxx/outs/filtered_feature_bc_matrix.h5')
adata_dict['R100_62_MS_F_Age62'] = sc.read_10x_h5('../../../../raw_counts_scRNA/xxx/outs/filtered_feature_bc_matrix.h5')
adata_dict['R100_61_MS_F_Age75'] = sc.read_10x_h5('../../../../raw_counts_scRNA/xxx/outs/filtered_feature_bc_matrix.h5')

adata_dict['R108_30_Ctrl_M_Age88'] = sc.read_10x_h5('../../../../../raw_data/xxx/outs/filtered_feature_bc_matrix.h5')
adata_dict['R108_31_Ctrl_F_Age82'] = sc.read_10x_h5('../../../../../raw_data/xxx/outs/filtered_feature_bc_matrix.h5')


# In[ ]:


for sample in adata_dict.keys():
    adata_dict[sample].var_names_make_unique()


# # Add metadata

# In[ ]:


#Spain_18_006_Ctrl_M
adata_dict['Spain_18_006_Ctrl_M_Age51'].obs['Sample'] = "18TN00006"
adata_dict['Spain_18_006_Ctrl_M_Age51'].obs['Condition'] = 'Control'
adata_dict['Spain_18_006_Ctrl_M_Age51'].obs ['Ventricle'] = 'LV'
adata_dict['Spain_18_006_Ctrl_M_Age51'].obs['Gender'] ='Male'
adata_dict['Spain_18_006_Ctrl_M_Age51'].obs['Age'] =51
adata_dict['Spain_18_006_Ctrl_M_Age51'].obs['project'] = 'Spain'

#Spain_19_005_Ctrl_M
adata_dict['Spain_19_005_Ctrl_M_Age60'].obs['Sample'] = "19TN00005"
adata_dict['Spain_19_005_Ctrl_M_Age60'].obs['Condition'] = 'Control'
adata_dict['Spain_19_005_Ctrl_M_Age60'].obs['Gender'] ='Male'
adata_dict['Spain_19_005_Ctrl_M_Age60'].obs['Age'] =60
adata_dict['Spain_19_005_Ctrl_M_Age60'].obs ['Ventricle'] = 'LV'
adata_dict['Spain_19_005_Ctrl_M_Age60'].obs['project'] = 'Spain'


#R100_64_Ctrl_M
adata_dict['R100_64_Ctrl_M_Age64'].obs['Sample'] = "1468"
adata_dict['R100_64_Ctrl_M_Age64'].obs['Condition'] = 'Control'
adata_dict['R100_64_Ctrl_M_Age64'].obs ['Ventricle'] = 'LV'
adata_dict['R100_64_Ctrl_M_Age64'].obs['Gender'] ='Male'
adata_dict['R100_64_Ctrl_M_Age64'].obs['Age'] =64
adata_dict['R100_64_Ctrl_M_Age64'].obs['project'] = 'R100'


#R77_61_Ctrl_F
adata_dict['R77_61_Ctrl_F_Age79'].obs['Sample'] = "CS2255"
adata_dict['R77_61_Ctrl_F_Age79'].obs['Condition'] = 'Control'
adata_dict['R77_61_Ctrl_F_Age79'].obs ['Ventricle'] = 'LV'
adata_dict['R77_61_Ctrl_F_Age79'].obs['Gender'] ='Female'
adata_dict['R77_61_Ctrl_F_Age79'].obs['Age'] =79
adata_dict['R77_61_Ctrl_F_Age79'].obs['project'] = 'R77'

#R100_72_Ctrl_F
adata_dict['R100_72_Ctrl_F_Age82'].obs['Sample'] = "18TN00005"
adata_dict['R100_72_Ctrl_F_Age82'].obs['Condition'] = 'Control'
adata_dict['R100_72_Ctrl_F_Age82'].obs ['Ventricle'] = 'LV'
adata_dict['R100_72_Ctrl_F_Age82'].obs['Gender'] ='Female'
adata_dict['R100_72_Ctrl_F_Age82'].obs['Age'] =82
adata_dict['R100_72_Ctrl_F_Age82'].obs['project'] = 'R100'

#R100_63_Ctrl_F
adata_dict['R100_63_Ctrl_F_Age56'].obs['Sample'] = "1541"
adata_dict['R100_63_Ctrl_F_Age56'].obs['Condition'] = 'Control'
adata_dict['R100_63_Ctrl_F_Age56'].obs ['Ventricle'] = 'LV'
adata_dict['R100_63_Ctrl_F_Age56'].obs['Gender'] ='Female'
adata_dict['R100_63_Ctrl_F_Age56'].obs['Age'] =56
adata_dict['R100_63_Ctrl_F_Age56'].obs['project'] = 'R100'

#R77_60_MS_M
adata_dict['R77_60_MS_M_Age46'].obs['Sample'] = "CS639"
adata_dict['R77_60_MS_M_Age46'].obs['Condition'] = 'MS'
adata_dict['R77_60_MS_M_Age46'].obs ['Ventricle'] = 'LV'
adata_dict['R77_60_MS_M_Age46'].obs['Gender'] ='Male'
adata_dict['R77_60_MS_M_Age46'].obs['Age'] =46
adata_dict['R77_60_MS_M_Age46'].obs['project'] = 'R77'

#R100_71_MS_F
adata_dict['R100_71_MS_F_Age52'].obs['Sample'] = "MS921"
adata_dict['R100_71_MS_F_Age52'].obs['Condition'] = 'MS'
adata_dict['R100_71_MS_F_Age52'].obs ['Ventricle'] = 'LV'
adata_dict['R100_71_MS_F_Age52'].obs['Gender'] ='Female'
adata_dict['R100_71_MS_F_Age52'].obs['Age'] =52
adata_dict['R100_71_MS_F_Age52'].obs['project'] = 'R100'

#R100_62_MS_F
adata_dict['R100_62_MS_F_Age62'].obs['Sample'] = "MS748"
adata_dict['R100_62_MS_F_Age62'].obs['Condition'] = 'MS'
adata_dict['R100_62_MS_F_Age62'].obs ['Ventricle'] = 'LV'
adata_dict['R100_62_MS_F_Age62'].obs['Gender'] ='Female'
adata_dict['R100_62_MS_F_Age62'].obs['Age'] =62
adata_dict['R100_62_MS_F_Age62'].obs['project'] = 'R100'

#R100_61_MS_F
adata_dict['R100_61_MS_F_Age75'].obs['Sample'] = "MS745"
adata_dict['R100_61_MS_F_Age75'].obs['Condition'] = 'MS'
adata_dict['R100_61_MS_F_Age75'].obs ['Ventricle'] = 'LV'
adata_dict['R100_61_MS_F_Age75'].obs['Gender'] ='Female'
adata_dict['R100_61_MS_F_Age75'].obs['Age'] =75
adata_dict['R100_61_MS_F_Age75'].obs['project'] = 'R100'

#R108_30_Ctrl_M_Age88
adata_dict['R108_30_Ctrl_M_Age88'].obs['Sample'] = "PDC164"
adata_dict['R108_30_Ctrl_M_Age88'].obs['Condition'] = 'Control'
adata_dict['R108_30_Ctrl_M_Age88'].obs ['Ventricle'] = 'LV'
adata_dict['R108_30_Ctrl_M_Age88'].obs['Gender'] ='Male'
adata_dict['R108_30_Ctrl_M_Age88'].obs['Age'] = 88
adata_dict['R108_30_Ctrl_M_Age88'].obs['project'] = 'R108'

#R108_31_Ctrl_NA_AgeNA
adata_dict['R108_31_Ctrl_F_Age82'].obs['Sample'] = "18TN00005"
adata_dict['R108_31_Ctrl_F_Age82'].obs['Condition'] = 'Control'
adata_dict['R108_31_Ctrl_F_Age82'].obs ['Ventricle'] = 'LV'
adata_dict['R108_31_Ctrl_F_Age82'].obs['Gender'] ='Female'
adata_dict['R108_31_Ctrl_F_Age82'].obs['Age'] = 82
adata_dict['R108_31_Ctrl_F_Age82'].obs['project'] = 'R108'


# In[ ]:


for sample in adata_dict.keys():
    adata_dict[sample].obs['Sample'] = adata_dict[sample].obs['Sample'].astype('category')
    adata_dict[sample].obs['Condition'] = adata_dict[sample].obs['Condition'].astype('category')
    adata_dict[sample].obs['Gender'] = adata_dict[sample].obs['Gender'].astype('category')
    adata_dict[sample].obs['project'] = adata_dict[sample].obs['project'].astype('category')


# # Raw counts

# In[ ]:


adata_merge = adata_dict['Spain_18_006_Ctrl_M_Age51'].concatenate(
    adata_dict['R100_63_Ctrl_F_Age56'],adata_dict['Spain_19_005_Ctrl_M_Age60'],adata_dict['R100_64_Ctrl_M_Age64'],adata_dict['R77_61_Ctrl_F_Age79'],
    adata_dict['R100_72_Ctrl_F_Age82'],adata_dict['R108_30_Ctrl_M_Age88'],adata_dict['R108_31_Ctrl_F_Age82'],adata_dict['R77_60_MS_M_Age46'],
    adata_dict['R100_71_MS_F_Age52'],adata_dict['R100_62_MS_F_Age62'],adata_dict['R100_61_MS_F_Age75'],
    batch_categories=['Spain_18_006_Ctrl_M_Age51','R100_63_Ctrl_F_Age56','Spain_19_005_Ctrl_M_Age60','R100_64_Ctrl_M_Age64','R77_61_Ctrl_F_Age79',
                      'R100_72_Ctrl_F_Age82','R108_30_Ctrl_M_Age88','R108_31_Ctrl_F_Age82','R77_60_MS_M_Age46','R100_71_MS_F_Age52','R100_62_MS_F_Age62', 'R100_61_MS_F_Age75'])


# In[ ]:


#add mito, ribo and hemo gene info
adata_merge.var['Mito'] = adata_merge.var_names.str.startswith('MT-') 
adata_merge.var['Ribo'] = adata_merge.var_names.str.startswith(("RPS","RPL"))
adata_merge.var['Hb'] = adata_merge.var_names.str.contains(("^HB[^(P)]"))
adata_merge.obs['n_counts'] = adata_merge.X.sum(axis = 1).A1
adata_merge.obs['n_genes'] = np.sum(adata_merge.X > 0, 1)
sc.pp.calculate_qc_metrics(adata_merge, qc_vars=['Mito','Ribo','Hb'], percent_top=None, log1p=False, inplace=True)

#Sample
annot = sc.queries.biomart_annotations(
        "hsapiens",
        ["ensembl_gene_id", "external_gene_name", "start_position", "end_position", "chromosome_name"],
    ).set_index("external_gene_name")
    
#sample sex
chrY_genes = adata_merge.var_names.intersection(annot.index[annot.chromosome_name == "Y"])
adata_merge.obs['percent_chrY'] = np.sum(adata_merge[:, chrY_genes].X, axis=1).A1 / np.sum(adata_merge.X, axis=1).A1 * 100
adata_merge.obs["XIST-counts"] = adata_merge.X[:,adata_merge.var_names.str.match('XIST')].toarray()



# In[ ]:


adata_merge.write_h5ad('../results/Human_unfiltered.h5ad')


# ## Plot QC before filtering

# In[ ]:


sc.pl.violin(adata_merge,keys=['n_counts','n_genes','pct_counts_Mito','pct_counts_Ribo','pct_counts_Hb'],
                 rotation=90,stripplot=False,jitter=0.4,color = "#006837",multi_panel=True,groupby='batch')

sc.pl.violin(adata_merge, ["XIST-counts", "percent_chrY"], jitter=0, groupby = 'batch', rotation= 90)

#sns.set_context("paper", rc={"font.size":8,"axes.titlesize":12,"axes.labelsize":8})
ax = sns.countplot(y = 'batch',data = adata_merge.obs)
ax.set_ylabel('# The number of cells')

for container in ax.containers:
    ax.bar_label(container)
    
sc.pl.highest_expr_genes(adata_merge,n_top=20)



# In[ ]:


del adata_merge
gc.collect()


# # Run the pipeline individually

# ## Define filtering, DR and clustering function

# In[ ]:


def filtering_DR_clustering(adata,sample):
    #filtering
    print('Filtering for sample ',sample, ':')
    #calculate QC
    adata.var['Mito'] = adata.var_names.str.startswith("MT-")
    adata.var['Hb'] = adata.var_names.str.contains(("^HB[^(P)]"))
    adata.var ['Ribo'] = adata.var_names.str.startswith(("RPS","RPL"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['Mito','Ribo','Hb'], percent_top=None, log1p=False, inplace=True)
    ##add the total counts per cell as observations-annotation to adata
    adata.obs['n_counts'] = adata.X.sum(axis = 1).A1
    adata.obs['n_genes'] = np.sum(adata.X > 0, 1)
    
    # plot unfiltered QC
    #sc.pl.violin(adata,keys=['n_genes','n_counts','pct_counts_Mito','pct_counts_Ribo','pct_counts_Hb'],
    #             rotation=45,stripplot=False,jitter=0.4,color = "#006837",multi_panel=True)
    #plt.scatter(adata.obs['n_counts'],adata.obs['n_genes'],alpha = 0.5, color = "#006837")
    #plt.savefig('../plots/%s_correlation_CountGenes.pdf'%sample)
    sc.pl.highest_expr_genes(adata,n_top=20)
    
    #Filtering
    ##Low quality cells and genes
    print(adata.shape)
    sc.pp.filter_cells(adata,min_genes = 500)
    #sc.pp.filter_genes(adata,min_cells = 5)
    print(adata.n_obs,adata.n_vars)
    
    ## Mito/Hb filtering
    adata = adata[adata.obs['pct_counts_Mito']<2,:]
    adata = adata[adata.obs['pct_counts_Hb'] < 0.5,:]
    adata = adata[adata.obs['pct_counts_Ribo'] < 10,:]
    print("Remaining cells %d" % adata.n_obs)
    
    ## Filter MT genes and MALAT1
    malat1 = adata.var_names.str.startswith('MALAT1')
    mito_genes = adata.var_names.str.startswith('MT-')
    remove = np.add(mito_genes, malat1)
    keep = np.invert(remove)
    adata = adata[:,keep]
    print(adata.shape)
    
    #Calculate cell-cycle scores
    get_ipython().system('if [ ! -f data/regev_lab_cell_cycle_genes.txt ]; then curl -o regev_lab_cell_cycle_genes.txt https://raw.githubusercontent.com/theislab/scanpy_usage/master/180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt; fi')
    
    cell_cycle_genes = [x.strip() for x in open('regev_lab_cell_cycle_genes.txt')]
    print(len(cell_cycle_genes))
    ## Split into 2 lists
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
    print(len(cell_cycle_genes))

    adata.raw = adata
    ## normalize to depth 10 000
    sc.pp.normalize_per_cell(adata,counts_per_cell_after=1e4)
    ## logaritmize
    sc.pp.log1p(adata)
    ## scale
    sc.pp.scale(adata)
    ## calculating cell cycle phase
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    
    
    ##  Predict doublets
    import scrublet as scr
    scrub = scr.Scrublet(adata.raw.X)
    adata.obs['doublet_score'] ,adata.obs['predicted_doublets'] = scrub.scrub_doublets()
    scrub.plot_histogram()
    sum(adata.obs['predicted_doublets'])
    
    ## add in column with singlet/doublet instead of True/False
    adata.obs['doublet_info'] = adata.obs["predicted_doublets"].astype(str)
    sc.pl.violin(adata, 'n_genes_by_counts',
                 jitter=0.4, groupby = 'doublet_info', rotation=45)
    
    adata = adata.raw.to_adata() 
    adata = adata[adata.obs['doublet_info'] == 'False',:]
    print(adata.shape)
    
    
    ##Remove highly expressed genes
    sc.pp.filter_cells(adata,max_genes = 6000)
    sc.pp.filter_cells(adata,max_counts = 20000)
    
    #Plot filtered QC
    sc.pl.violin(adata,['n_genes','n_counts','pct_counts_Mito','pct_counts_Ribo', 'pct_counts_Hb'],
                 multi_panel=True,stripplot=False,color='#006837')
    print('The number of cells features after filtering is: \n',adata.shape)
    
    #return adata
    
    ##DR 
    # normalize to depth 10000
    sc.pp.normalize_per_cell(adata,counts_per_cell_after=1e4)
    #logaritmize
    sc.pp.log1p(adata)
    #store normalized counts in the raw slot, we will subset adata.X for variable genes, but want to keep all genes matrix as well
    adata.raw = adata
    
    #HVGs
    sc.pp.highly_variable_genes(adata,min_mean=0.0125, max_mean=3, min_disp=0.5)
    print('Highly variable genes: %d' %sum(adata.var.highly_variable))
    #plot variable genes
    sc.pl.highly_variable_genes(adata)
    
    #subset for variable genes in the dataset
    adata = adata[:,adata.var['highly_variable']]

    #Z-score transformation
    # regress out unwanted variables
    #sc.pp.regress_out(adata,['pct_counts_Mito'])
    # scale data, clip values exceeding standard deviation 10.
    sc.pp.scale(adata,max_value=10)
    
    #PCA
    sc.tl.pca(adata,svd_solver = 'arpack',n_comps = 50)
    sc.pl.pca_variance_ratio(adata,log = True,n_pcs = 50)
    
    #UMAP
    sc.pp.neighbors(adata,n_pcs = 20) #Calculate neighborhood graph
    sc.tl.umap(adata)
    
    # clustering
    sc.tl.louvain(adata, resolution = 0.1, key_added = "louvain_0.1") 
    sc.tl.louvain(adata, resolution = 0.2, key_added = "louvain_0.2")
    sc.tl.louvain(adata, resolution = 0.3, key_added = "louvain_0.3")
    sc.tl.louvain(adata, resolution = 0.4, key_added = "louvain_0.4")
    sc.tl.louvain(adata, resolution = 0.5, key_added = "louvain_0.5")
    sc.tl.louvain(adata, resolution = 0.6, key_added = "louvain_0.6")
    sc.tl.louvain(adata, resolution = 0.7, key_added = "louvain_0.7")
    sc.tl.louvain(adata, resolution = 0.8, key_added = "louvain_0.8")
    sc.tl.louvain(adata, key_added = "louvain_1.0") # default resolution in 1.0
    
    #plot
    print('The clustering of sample: ',sample)
    sc.pl.umap(adata,color=['louvain_0.1','louvain_0.2','louvain_0.3','louvain_0.4','louvain_0.5','louvain_0.6'],
               title=['louvain_0.1','louvain_0.2','louvain_0.3','louvain_0.4','louvain_0.5','louvain_0.6'],
               ncols=3, wspace=0.5)
    
    return adata
 


# ## Canonical markers

# In[ ]:


Epithelial_marker = ["OTX2","CLIC6","FOLR1","IGFBP2","PDGFA"]
Endothelial_marker = ["PLVAP","CD34","ICAM1","CLDN5"]
Immune_marker = ["AIF1","TYROBP","FCER1G","CD74",'CD68']
Mesenchyma_marker = ["LUM","COL1A1","COL1A2","COL3A1","CD44"]
T_cells_NK = ["SKAP1","PTPRC","MBNL1","B2M","CDC42SE2"]
Astrocyte = ['GFAP']
Oligos = ["FA2H","CNP","UGT8"]
Neuros = ["KCNIP4","NRXN1","KCND2","RIMS1","RIMS2","DSCAM","TNR","VCAN"]
OPCs = ['PDGFRA', 'SOX6','SOX10']


# ### Spain_18_006_Ctrl_M_Age51

# In[ ]:


Spain_18_006_Ctrl_M_Age51_adata = filtering_DR_clustering(adata_dict['Spain_18_006_Ctrl_M_Age51'],'Spain_18_006_Ctrl_M_Age51')


# In[ ]:


sc.pl.dotplot(Spain_18_006_Ctrl_M_Age51_adata, Epithelial_marker, groupby='louvain_0.4', dendrogram=True,title='Epithelial')
sc.pl.dotplot(Spain_18_006_Ctrl_M_Age51_adata, Endothelial_marker, groupby='louvain_0.4', dendrogram=True,title='Endothelial')
sc.pl.dotplot(Spain_18_006_Ctrl_M_Age51_adata, Mesenchyma_marker, groupby='louvain_0.4', dendrogram=True,title='Mesenchymal')
sc.pl.dotplot(Spain_18_006_Ctrl_M_Age51_adata, Immune_marker, groupby='louvain_0.4', dendrogram=True,title='Immune')
sc.pl.dotplot(Spain_18_006_Ctrl_M_Age51_adata, T_cells_NK, groupby='louvain_0.4', dendrogram=True,title='T_cells')
sc.pl.dotplot(Spain_18_006_Ctrl_M_Age51_adata, Neuros+Astrocyte+Oligos, groupby='louvain_0.4', dendrogram=True,title='Neurons+Astro+Oligos')

sc.pl.dotplot(Spain_18_006_Ctrl_M_Age51_adata, Epithelial_marker+Endothelial_marker+Mesenchyma_marker+Immune_marker+T_cells_NK+Neuros+Astrocyte+Oligos, groupby='louvain_0.4', dendrogram=True,
              title='Epithelial_marker+Endothelial_marker+Mesenchyma_marker+Immune_marker+T_cells_NK+Neuros+Astrocyte+Oligos')


# In[ ]:


Spain_18_006_Ctrl_M_Age51_adata.obs['CellType'] = Spain_18_006_Ctrl_M_Age51_adata.obs['louvain_0.4'].astype('category')
Spain_18_006_Ctrl_M_Age51_adata.obs['CellType'].replace({'0': 'Mesenchymal',
                                                      '1': 'Mesenchymal',
                                                      '2': 'Mesenchymal',
                                                      '3': 'Immune',
                                                      '4': 'Epithelial',
                                                      '5': 'Immune',
                                                      '6': 'Mesenchymal',
                                                      '7': 'Mesenchymal',
                                                         '8':'Endothelial'},
                                                     inplace=True)

sc.pl.umap(Spain_18_006_Ctrl_M_Age51_adata,color='CellType',title='Spain_18_006_Ctrl_M_Age51')

Spain_18_006_Ctrl_M_Age51_adata.write_h5ad('../results/individual_samples/Spain_18_006_Ctrl_M_Age51_CellType.h5ad')


# In[ ]:


sc.pl.umap(Spain_18_006_Ctrl_M_Age51_adata,color=['CLIC6','OTX2','PECAM1','CLDN5','COL3A1','ACTA2','PTPRC','VCAN'],ncols=4) 


# ### Spain_19_005_Ctrl_M_Age60

# In[ ]:


Spain_19_005_Ctrl_M_Age60_adata = filtering_DR_clustering(adata_dict['Spain_19_005_Ctrl_M_Age60'],'Spain_19_005_Ctrl_M_Age60')


# In[ ]:


sc.pl.dotplot(Spain_19_005_Ctrl_M_Age60_adata, Epithelial_marker, groupby='louvain_0.3', dendrogram=True,title='Epithelial')
sc.pl.dotplot(Spain_19_005_Ctrl_M_Age60_adata, Endothelial_marker, groupby='louvain_0.3', dendrogram=True,title='Endothelial')
sc.pl.dotplot(Spain_19_005_Ctrl_M_Age60_adata, Mesenchyma_marker, groupby='louvain_0.3', dendrogram=True,title='Mesenchymal')
sc.pl.dotplot(Spain_19_005_Ctrl_M_Age60_adata, Immune_marker, groupby='louvain_0.3', dendrogram=True,title='Immune')
sc.pl.dotplot(Spain_19_005_Ctrl_M_Age60_adata, T_cells_NK, groupby='louvain_0.3', dendrogram=True,title='T_cells')
sc.pl.dotplot(Spain_19_005_Ctrl_M_Age60_adata, Neuros+Astrocyte+Oligos, groupby='louvain_0.3', dendrogram=True,title='Neurons+Astro+Oligos')

sc.pl.dotplot(Spain_19_005_Ctrl_M_Age60_adata, Mesenchyma_marker+Neuros+Astrocyte+Oligos, groupby='louvain_0.3', dendrogram=True,title='Mesenchymal+Neurons+Astro+Oligos')


# In[ ]:


Spain_19_005_Ctrl_M_Age60_adata.obs['CellType'] = Spain_19_005_Ctrl_M_Age60_adata.obs['louvain_0.3'].astype('category')
Spain_19_005_Ctrl_M_Age60_adata.obs['CellType'].replace({'0': 'Mesenchymal','1': 'Mesenchymal','2': 'Mesenchymal','3': 'Mesenchymal',
                                                         '4': 'Immune','5': 'Epithelial','6': 'Mesenchymal',
                                                         '7':'Immune','8':'Endothelial'},
                                                     inplace=True)

sc.pl.umap(Spain_19_005_Ctrl_M_Age60_adata,color='CellType',title='Spain_19_005_Ctrl_M_Age60')

Spain_19_005_Ctrl_M_Age60_adata.write_h5ad('../results/individual_samples/Spain_19_005_Ctrl_M_Age60_CellType.h5ad')


# In[ ]:


sc.pl.umap(Spain_19_005_Ctrl_M_Age60_adata,color=['CLIC6','OTX2','PECAM1','CLDN5','COL3A1','ACTA2','PTPRC','VCAN'],ncols=4) 


# ### R77_61_Ctrl_F_Age79

# In[ ]:


R77_61_Ctrl_F_Age79_adata = filtering_DR_clustering(adata_dict['R77_61_Ctrl_F_Age79'],'R77_61_Ctrl_F_Age79')


# In[ ]:


sc.pl.dotplot(R77_61_Ctrl_F_Age79_adata, Epithelial_marker, groupby='louvain_0.2', dendrogram=True,title='Epithelial')
sc.pl.dotplot(R77_61_Ctrl_F_Age79_adata, Endothelial_marker, groupby='louvain_0.2', dendrogram=True,title='Endothelial')
sc.pl.dotplot(R77_61_Ctrl_F_Age79_adata, Mesenchyma_marker, groupby='louvain_0.2', dendrogram=True,title='Mesenchymal')
sc.pl.dotplot(R77_61_Ctrl_F_Age79_adata, Immune_marker, groupby='louvain_0.2', dendrogram=True,title='Immune')
sc.pl.dotplot(R77_61_Ctrl_F_Age79_adata, T_cells_NK, groupby='louvain_0.2', dendrogram=True,title='T_cells')
sc.pl.dotplot(R77_61_Ctrl_F_Age79_adata, Neuros+Astrocyte+Oligos, groupby='louvain_0.2', dendrogram=True,title='Neurons+Astro+Oligos')

sc.pl.dotplot(R77_61_Ctrl_F_Age79_adata, Epithelial_marker+Endothelial_marker+Mesenchyma_marker+Immune_marker+T_cells_NK+Neuros+Astrocyte+Oligos, 
              groupby='louvain_0.2', dendrogram=True,title='Epithelial_marker+Endothelial_marker+Mesenchyma_marker+Immune_marker+T_cells_NK+Neurons+Astro+Oligos')

#sc.pl.dotplot(Spain_19_005_Ctrl_M_adata, Mesenchyma_marker+Neuros+Astrocyte+Oligos, groupby='louvain_0.2', dendrogram=True,title='Mesenchymal+Neurons+Astro+Oligos')


# In[ ]:


R77_61_Ctrl_F_Age79_adata.obs['CellType'] = R77_61_Ctrl_F_Age79_adata.obs['louvain_0.2'].astype('category')
R77_61_Ctrl_F_Age79_adata.obs['CellType'].replace({'0': 'Mesenchymal','1': 'Immune','2': 'Epithelial','3': 'Immune',
                                                   '4': 'Mesenchymal','5': 'Immune','6': 'NeuronGlial',
                                                   '7': 'NeuronGlial','8': 'NeuronGlial','9':'Mesenchymal'},
                                               inplace=True)
#cluster9 is the doublets of pericyte and endothelial
sc.pl.umap(R77_61_Ctrl_F_Age79_adata,color='CellType',title='R77_61_Ctrl_F_Age79')

R77_61_Ctrl_F_Age79_adata.write_h5ad('../results/individual_samples/R77_61_Ctrl_F_Age79_CellType.h5ad')


# In[ ]:


sc.pl.umap(R77_61_Ctrl_F_Age79_adata,color=['CLIC6','OTX2','PECAM1','CLDN5','COL3A1','ACTA2','PTPRC','VCAN'],ncols=4) 


# ### R100_63_Ctrl_F_Age56

# In[ ]:


R100_63_Ctrl_F_Age56_adata = filtering_DR_clustering(adata_dict['R100_63_Ctrl_F_Age56'],'R100_63_Ctrl_F_Age56')


# In[ ]:


sc.pl.dotplot(R100_63_Ctrl_F_Age56_adata, Epithelial_marker, groupby='louvain_0.6', dendrogram=True,title='Epithelial')
sc.pl.dotplot(R100_63_Ctrl_F_Age56_adata, Endothelial_marker, groupby='louvain_0.6', dendrogram=True,title='Endothelial')
sc.pl.dotplot(R100_63_Ctrl_F_Age56_adata, Mesenchyma_marker, groupby='louvain_0.6', dendrogram=True,title='Mesenchymal')
sc.pl.dotplot(R100_63_Ctrl_F_Age56_adata, Immune_marker, groupby='louvain_0.6', dendrogram=True,title='Immune')
sc.pl.dotplot(R100_63_Ctrl_F_Age56_adata, T_cells_NK, groupby='louvain_0.6', dendrogram=True,title='T_cells')
sc.pl.dotplot(R100_63_Ctrl_F_Age56_adata, Neuros+Astrocyte+Oligos, groupby='louvain_0.6', dendrogram=True,title='Neurons+Astro+Oligos')

sc.pl.dotplot(R100_63_Ctrl_F_Age56_adata, Epithelial_marker+Endothelial_marker+Mesenchyma_marker+Immune_marker+T_cells_NK+Neuros+Astrocyte+Oligos, groupby='louvain_0.6', dendrogram=True,
              title='Epithelial_marker+Endothelial_marker+Mesenchyma_marker+Immune_marker+T_cells_NK+Neuros+Astrocyte+Oligos')


# In[ ]:


R100_63_Ctrl_F_Age56_adata.obs['CellType'] = R100_63_Ctrl_F_Age56_adata.obs['louvain_0.6'].astype('category')
R100_63_Ctrl_F_Age56_adata.obs['CellType'].replace({'0': 'Mesenchymal',
                                                '1': 'Epithelial',
                                                '2': 'Mesenchymal',
                                                '3': 'Mesenchymal',
                                                '4': 'Immune',
                                                '5': 'Mesenchymal',
                                                '6': 'Immune',
                                                '7': 'Endothelial'},
                                               inplace=True)

sc.pl.umap(R100_63_Ctrl_F_Age56_adata,color='CellType',title='R100_63_Ctrl_F_Age56')

R100_63_Ctrl_F_Age56_adata.write_h5ad('../results/individual_samples/R100_63_Ctrl_F_Age56_CellType.h5ad')


# In[ ]:


sc.pl.umap(R100_63_Ctrl_F_Age56_adata,color=['CLIC6','OTX2','PECAM1','CLDN5','COL3A1','ACTA2','PTPRC','VCAN'],ncols=4) 


# ### R100_64_Ctrl_M_Age64

# In[ ]:


R100_64_Ctrl_M_Age64_adata = filtering_DR_clustering(adata_dict['R100_64_Ctrl_M_Age64'],'R100_64_Ctrl_M_Age64')


# In[ ]:


sc.pl.dotplot(R100_64_Ctrl_M_Age64_adata, Epithelial_marker, groupby='louvain_0.2', dendrogram=True,title='Epithelial')
sc.pl.dotplot(R100_64_Ctrl_M_Age64_adata, Endothelial_marker, groupby='louvain_0.2', dendrogram=True,title='Endothelial')
sc.pl.dotplot(R100_64_Ctrl_M_Age64_adata, Mesenchyma_marker, groupby='louvain_0.2', dendrogram=True,title='Mesenchymal')
sc.pl.dotplot(R100_64_Ctrl_M_Age64_adata, Immune_marker, groupby='louvain_0.2', dendrogram=True,title='Immune')
sc.pl.dotplot(R100_64_Ctrl_M_Age64_adata, T_cells_NK, groupby='louvain_0.2', dendrogram=True,title='T_cells')
sc.pl.dotplot(R100_64_Ctrl_M_Age64_adata, Neuros+Astrocyte+Oligos, groupby='louvain_0.2', dendrogram=True,title='Neurons+Astro+Oligos')

#sc.pl.dotplot(Spain_19_005_Ctrl_M_adata, Mesenchyma_marker+Neuros+Astrocyte+Oligos, groupby='louvain_0.2', dendrogram=True,title='Mesenchymal+Neurons+Astro+Oligos')


# In[ ]:


R100_64_Ctrl_M_Age64_adata.obs['CellType'] = R100_64_Ctrl_M_Age64_adata.obs['louvain_0.2'].astype('category')
R100_64_Ctrl_M_Age64_adata.obs['CellType'].replace({'0': 'Epithelial',
                                                '1': 'NeuronGlial',
                                                '2': 'Immune',
                                                '3': 'Mesenchymal',
                                                '4': 'NeuronGlial',
                                                '5': 'Immune',
                                                '6': 'Epithelial',
                                                '7': 'Endothelial',
                                                '8': 'Mesenchymal',
                                                '9': 'NeuronGlial',
                                                '10': 'NeuronGlial'},
                                               inplace=True)

sc.pl.umap(R100_64_Ctrl_M_Age64_adata,color='CellType',title='R100_64_Ctrl_M_Age64')

R100_64_Ctrl_M_Age64_adata.write_h5ad('../results/individual_samples/R100_64_Ctrl_M_Age64_CellType.h5ad')


# In[ ]:


sc.pl.umap(R100_64_Ctrl_M_Age64_adata,color=['CLIC6','OTX2','PECAM1','CLDN5','COL3A1','ACTA2','PTPRC','VCAN'],ncols=4) 


# ### R100_72_Ctrl_F_Age82

# In[ ]:


R100_72_Ctrl_F_Age82_adata = filtering_DR_clustering(adata_dict['R100_72_Ctrl_F_Age82'],'R100_72_Ctrl_F_Age82')


# In[ ]:


sc.pl.dotplot(R100_72_Ctrl_F_Age82_adata, Epithelial_marker, groupby='louvain_0.3', dendrogram=True,title='Epithelial')
sc.pl.dotplot(R100_72_Ctrl_F_Age82_adata, ["PLVAP","ICAM1","CLDN5"], groupby='louvain_0.3', dendrogram=True,title='Endothelial')
sc.pl.dotplot(R100_72_Ctrl_F_Age82_adata, Mesenchyma_marker, groupby='louvain_0.3', dendrogram=True,title='Mesenchymal')
sc.pl.dotplot(R100_72_Ctrl_F_Age82_adata, Immune_marker, groupby='louvain_0.3', dendrogram=True,title='Immune')
sc.pl.dotplot(R100_72_Ctrl_F_Age82_adata, T_cells_NK, groupby='louvain_0.3', dendrogram=True,title='T_cells')
sc.pl.dotplot(R100_72_Ctrl_F_Age82_adata, Neuros+Astrocyte+Oligos, groupby='louvain_0.3', dendrogram=True,title='Neurons+Astro+Oligos')

#sc.pl.dotplot(Spain_19_005_Ctrl_M_adata, Mesenchyma_marker+Neuros+Astrocyte+Oligos, groupby='louvain_0.2', dendrogram=True,title='Mesenchymal+Neurons+Astro+Oligos')


# In[ ]:


R100_72_Ctrl_F_Age82_adata.obs['CellType'] = R100_72_Ctrl_F_Age82_adata.obs['louvain_0.3'].astype('category')
R100_72_Ctrl_F_Age82_adata.obs['CellType'].replace({'0': 'Epithelial',
                                                    '1': 'Immune',
                                                    '2': 'Epithelial',
                                                    '3': 'Mesenchymal',
                                                    '4': 'Endothelial',
                                                    '5': 'Immune',
                                                    '6': 'NeuronGlial',
                                                    '7': 'Mesenchymal'},
                                               inplace=True)

sc.pl.umap(R100_72_Ctrl_F_Age82_adata,color='CellType',title='R100_72_Ctrl_F_Age82')

R100_72_Ctrl_F_Age82_adata.write_h5ad('../results/individual_samples/R100_72_Ctrl_F_Age82_CellType.h5ad')


# In[ ]:


sc.pl.umap(R100_72_Ctrl_F_Age82_adata,color=['CLIC6','OTX2','PECAM1','CLDN5','COL3A1','ACTA2','PTPRC','VCAN'],ncols=4) 


# ### R108_30_Ctrl_M_Age88

# In[ ]:


R108_30_Ctrl_M_Age88_adata = filtering_DR_clustering(adata_dict['R108_30_Ctrl_M_Age88'],'R108_30_Ctrl_M_Age88')


# In[ ]:


sc.pl.dotplot(R108_30_Ctrl_M_Age88_adata, Epithelial_marker, groupby='louvain_0.3', dendrogram=True,title='Epithelial')
sc.pl.dotplot(R108_30_Ctrl_M_Age88_adata, ["PLVAP","ICAM1","CLDN5"], groupby='louvain_0.3', dendrogram=True,title='Endothelial')
sc.pl.dotplot(R108_30_Ctrl_M_Age88_adata, Mesenchyma_marker, groupby='louvain_0.3', dendrogram=True,title='Mesenchymal')
sc.pl.dotplot(R108_30_Ctrl_M_Age88_adata, Immune_marker, groupby='louvain_0.3', dendrogram=True,title='Immune')
sc.pl.dotplot(R108_30_Ctrl_M_Age88_adata, T_cells_NK, groupby='louvain_0.3', dendrogram=True,title='T_cells')
sc.pl.dotplot(R108_30_Ctrl_M_Age88_adata, Neuros+Astrocyte+Oligos, groupby='louvain_0.3', dendrogram=True,title='Neurons+Astro+Oligos')

sc.pl.dotplot(R108_30_Ctrl_M_Age88_adata, Epithelial_marker+["PLVAP","ICAM1","CLDN5"]+Mesenchyma_marker+Immune_marker+T_cells_NK+Neuros+Astrocyte+Oligos, groupby='louvain_0.3', dendrogram=True,title='All cell types')


# In[ ]:


R108_30_Ctrl_M_Age88_adata.obs['CellType'] = R108_30_Ctrl_M_Age88_adata.obs['louvain_0.3'].astype('category')
R108_30_Ctrl_M_Age88_adata.obs['CellType'].replace({'0': 'Epithelial','1': 'Epithelial','2': 'Immune','3': 'Epithelial','4': 'Mesenchymal','5': 'NeuronGlial','6': 'NeuronGlial','7': 'NeuronGlial','8':'Endothelial','9':'Immune'},
                                               inplace=True)

sc.pl.umap(R108_30_Ctrl_M_Age88_adata,color='CellType',title='R108_30_Ctrl_M_Age88')

R108_30_Ctrl_M_Age88_adata.write_h5ad('../results/individual_samples/R108_30_Ctrl_M_Age88_CellType.h5ad')


# In[ ]:


# canonical markers expression
sc.pl.umap(R108_30_Ctrl_M_Age88_adata,color=['CLIC6','OTX2','PECAM1','CLDN5','COL3A1','ACTA2','PTPRC','VCAN','CNP','GFAP','RIMS2','SOX10'],ncols=4) 


# ### R108_31_Ctrl_F_Age82

# In[ ]:


R108_31_Ctrl_F_Age82_adata = filtering_DR_clustering(adata_dict['R108_31_Ctrl_F_Age82'],'R108_31_Ctrl_F_Age82')


# In[ ]:


sc.pl.dotplot(R108_31_Ctrl_F_Age82_adata, Epithelial_marker, groupby='louvain_0.4', dendrogram=True,title='Epithelial')
sc.pl.dotplot(R108_31_Ctrl_F_Age82_adata, ["PLVAP","ICAM1","CLDN5"], groupby='louvain_0.4', dendrogram=True,title='Endothelial')
sc.pl.dotplot(R108_31_Ctrl_F_Age82_adata, Mesenchyma_marker, groupby='louvain_0.4', dendrogram=True,title='Mesenchymal')
sc.pl.dotplot(R108_31_Ctrl_F_Age82_adata, Immune_marker, groupby='louvain_0.4', dendrogram=True,title='Immune')
sc.pl.dotplot(R108_31_Ctrl_F_Age82_adata, T_cells_NK, groupby='louvain_0.4', dendrogram=True,title='T_cells')
sc.pl.dotplot(R108_31_Ctrl_F_Age82_adata, Neuros+Astrocyte+Oligos, groupby='louvain_0.4', dendrogram=True,title='Neurons+Astro+Oligos')

sc.pl.dotplot(R108_31_Ctrl_F_Age82_adata, Epithelial_marker+["PLVAP","ICAM1","CLDN5"]+Mesenchyma_marker+Immune_marker+T_cells_NK+Neuros+Astrocyte+Oligos, groupby='louvain_0.4', dendrogram=True,title='All cell types')


# In[ ]:


R108_31_Ctrl_F_Age82_adata.obs['CellType'] = R108_31_Ctrl_F_Age82_adata.obs['louvain_0.4'].astype('category')
R108_31_Ctrl_F_Age82_adata.obs['CellType'].replace({'0': 'Epithelial','1': 'Epithelial','2': 'Epithelial','3': 'Immune','4': 'Mesenchymal','5': 'Immune',
                                                    '6': 'Endothelial','7': 'NeuronGlial','8':'Mesenchymal','9':'NeuronGlial','10':'Epithelial'},
                                               inplace=True)

sc.pl.umap(R108_31_Ctrl_F_Age82_adata,color='CellType',title='R108_31_Ctrl_F_Age82')

R108_31_Ctrl_F_Age82_adata.write_h5ad('../results/individual_samples/R108_31_Ctrl_F_Age82_CellType.h5ad')


# In[ ]:


# canonical markers expression
sc.pl.umap(R108_31_Ctrl_F_Age82_adata,color=['CLIC6','OTX2','PECAM1','CLDN5','COL3A1','ACTA2','PTPRC','VCAN','CNP','GFAP','RIMS2','SOX10'],ncols=4) 


# ### R77_60_MS_M_Age46

# scrub_doublets doesn't perform well, use call_doublets(threshold=new_threshold) instead

# In[ ]:


def filtering_DR_clustering_v2(adata,sample):
    #filtering
    print('Filtering for sample ',sample, ':')
    #calculate QC
    adata.var['Mito'] = adata.var_names.str.startswith("MT-")
    adata.var['Hb'] = adata.var_names.str.contains(("^HB[^(P)]"))
    adata.var ['Ribo'] = adata.var_names.str.startswith(("RPS","RPL"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['Mito','Ribo','Hb'], percent_top=None, log1p=False, inplace=True)
    ##add the total counts per cell as observations-annotation to adata
    adata.obs['n_counts'] = adata.X.sum(axis = 1).A1
    adata.obs['n_genes'] = np.sum(adata.X > 0, 1)
    
    # plot unfiltered QC
    #sc.pl.violin(adata,keys=['n_genes','n_counts','pct_counts_Mito','pct_counts_Ribo','pct_counts_Hb'],
    #             rotation=45,stripplot=False,jitter=0.4,color = "#006837",multi_panel=True)
    #plt.scatter(adata.obs['n_counts'],adata.obs['n_genes'],alpha = 0.5, color = "#006837")
    #plt.savefig('../plots/%s_correlation_CountGenes.pdf'%sample)
    sc.pl.highest_expr_genes(adata,n_top=20)
    
    #Filtering
    ##Low quality cells and genes
    print(adata.shape)
    sc.pp.filter_cells(adata,min_genes = 500)
    #sc.pp.filter_genes(adata,min_cells = 5)
    print(adata.n_obs,adata.n_vars)
    
    ## Mito/Hb filtering
    adata = adata[adata.obs['pct_counts_Mito']<2,:]
    adata = adata[adata.obs['pct_counts_Hb'] < 0.5,:]
    adata = adata[adata.obs['pct_counts_Ribo'] < 10,:]
    print("Remaining cells %d" % adata.n_obs)
    
    ## Filter MT genes and MALAT1
    malat1 = adata.var_names.str.startswith('MALAT1')
    mito_genes = adata.var_names.str.startswith('MT-')
    remove = np.add(mito_genes, malat1)
    keep = np.invert(remove)
    adata = adata[:,keep]
    print(adata.shape)
    
    #Calculate cell-cycle scores
    get_ipython().system('if [ ! -f data/regev_lab_cell_cycle_genes.txt ]; then curl -o regev_lab_cell_cycle_genes.txt https://raw.githubusercontent.com/theislab/scanpy_usage/master/180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt; fi')
    
    cell_cycle_genes = [x.strip() for x in open('regev_lab_cell_cycle_genes.txt')]
    print(len(cell_cycle_genes))
    ## Split into 2 lists
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
    print(len(cell_cycle_genes))

    adata.raw = adata
    ## normalize to depth 10 000
    sc.pp.normalize_per_cell(adata,counts_per_cell_after=1e4)
    ## logaritmize
    sc.pp.log1p(adata)
    ## scale
    sc.pp.scale(adata)
    ## calculating cell cycle phase
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    
    
    ##  Predict doublets
    import scrublet as scr
    scrub = scr.Scrublet(adata.raw.X)
    adata.obs['doublet_score'] ,doublets = scrub.scrub_doublets()
    adata.obs['predicted_doublets'] = scrub.call_doublets(threshold=0.32)
    scrub.plot_histogram()
    sum(adata.obs['predicted_doublets'])
    
    ## add in column with singlet/doublet instead of True/False
    adata.obs['doublet_info'] = adata.obs["predicted_doublets"].astype(str)
    sc.pl.violin(adata, 'n_genes_by_counts',
                 jitter=0.4, groupby = 'doublet_info', rotation=45)
    
    adata = adata.raw.to_adata() 
    adata = adata[adata.obs['doublet_info'] == 'False',:]
    print(adata.shape)
    
    
    ##Remove highly expressed genes
    sc.pp.filter_cells(adata,max_genes = 6000)
    sc.pp.filter_cells(adata,max_counts = 20000)
    
    #Plot filtered QC
    sc.pl.violin(adata,['n_genes','n_counts','pct_counts_Mito','pct_counts_Ribo', 'pct_counts_Hb'],
                 multi_panel=True,stripplot=False,color='#006837')
    print('The number of cells features after filtering is: \n',adata.shape)
    
    #return adata
    
    ##DR 
    # normalize to depth 10000
    sc.pp.normalize_per_cell(adata,counts_per_cell_after=1e4)
    #logaritmize
    sc.pp.log1p(adata)
    #store normalized counts in the raw slot, we will subset adata.X for variable genes, but want to keep all genes matrix as well
    adata.raw = adata
    
    #HVGs
    sc.pp.highly_variable_genes(adata,min_mean=0.0125, max_mean=3, min_disp=0.5)
    print('Highly variable genes: %d' %sum(adata.var.highly_variable))
    #plot variable genes
    sc.pl.highly_variable_genes(adata)
    
    #subset for variable genes in the dataset
    adata = adata[:,adata.var['highly_variable']]

    #Z-score transformation
    # regress out unwanted variables
    #sc.pp.regress_out(adata,['pct_counts_Mito'])
    # scale data, clip values exceeding standard deviation 10.
    sc.pp.scale(adata,max_value=10)
    
    #PCA
    sc.tl.pca(adata,svd_solver = 'arpack',n_comps = 50)
    sc.pl.pca_variance_ratio(adata,log = True,n_pcs = 50)
    
    #UMAP
    sc.pp.neighbors(adata,n_pcs = 20) #Calculate neighborhood graph
    sc.tl.umap(adata)
    
    # clustering
    sc.tl.louvain(adata, resolution = 0.1, key_added = "louvain_0.1") 
    sc.tl.louvain(adata, resolution = 0.2, key_added = "louvain_0.2")
    sc.tl.louvain(adata, resolution = 0.3, key_added = "louvain_0.3")
    sc.tl.louvain(adata, resolution = 0.4, key_added = "louvain_0.4")
    sc.tl.louvain(adata, resolution = 0.5, key_added = "louvain_0.5")
    sc.tl.louvain(adata, resolution = 0.6, key_added = "louvain_0.6")
    sc.tl.louvain(adata, resolution = 0.7, key_added = "louvain_0.7")
    sc.tl.louvain(adata, resolution = 0.8, key_added = "louvain_0.8")
    sc.tl.louvain(adata, key_added = "louvain_1.0") # default resolution in 1.0
    
    #plot
    print('The clustering of sample: ',sample)
    sc.pl.umap(adata,color=['louvain_0.1','louvain_0.2','louvain_0.3','louvain_0.4',
                            'louvain_0.5','louvain_0.6','louvain_0.7','louvain_0.8','louvain_1.0'],
               title=['louvain_0.1','louvain_0.2','louvain_0.3','louvain_0.4','louvain_0.5',
                      'louvain_0.6','louvain_0.7','louvain_0.8','louvain_1.0'],
               ncols=3, wspace=0.5)
    
    return adata
 
    


# In[ ]:


R77_60_MS_M_Age46_adata = filtering_DR_clustering_v2(adata_dict['R77_60_MS_M_Age46'],'R77_60_MS_M_Age46')


# In[ ]:


sc.pl.dotplot(R77_60_MS_M_Age46_adata, Epithelial_marker, groupby='louvain_0.4', dendrogram=True,title='Epithelial')
sc.pl.dotplot(R77_60_MS_M_Age46_adata, Endothelial_marker, groupby='louvain_0.4', dendrogram=True,title='Endothelial')
sc.pl.dotplot(R77_60_MS_M_Age46_adata, Mesenchyma_marker, groupby='louvain_0.4', dendrogram=True,title='Mesenchymal')
sc.pl.dotplot(R77_60_MS_M_Age46_adata, Immune_marker, groupby='louvain_0.4', dendrogram=True,title='Immune')
sc.pl.dotplot(R77_60_MS_M_Age46_adata, T_cells_NK, groupby='louvain_0.4', dendrogram=True,title='T_cells')
sc.pl.dotplot(R77_60_MS_M_Age46_adata, Neuros+Astrocyte+Oligos, groupby='louvain_0.4', dendrogram=True,title='Neurons+Astro+Oligos')
sc.pl.dotplot(R77_60_MS_M_Age46_adata, Epithelial_marker+Neuros+Astrocyte+Oligos, groupby='louvain_0.4', dendrogram=True,title='Epithelial+Neurons+Astro+Oligos')


# In[ ]:


R77_60_MS_M_Age46_adata.obs['CellType'] = R77_60_MS_M_Age46_adata.obs['louvain_0.4'].astype('category')
R77_60_MS_M_Age46_adata.obs['CellType'].replace({'0': 'Epithelial','1': 'Epithelial','2': 'Immune','3': 'NeuronGlial','4': 'NeuronGlial','5': 'Mesenchymal',
                                                 '6': 'Endothelial','7': 'Immune','8':'NeuronGlial','9':'Mesenchymal','10':'NeuronGlial','11':'Mesenchymal'},
                                               inplace=True)

sc.pl.umap(R77_60_MS_M_Age46_adata,color='CellType',title='R77_60_MS_M_Age46')

R77_60_MS_M_Age46_adata.write_h5ad('../results/individual_samples/R77_60_MS_M_Age46_CellType.h5ad')


# In[ ]:


# canonical markers expression
sc.pl.umap(R77_60_MS_M_Age46_adata,color=['CLIC6','OTX2','PECAM1','CLDN5','COL3A1','ACTA2','PTPRC','VCAN','CNP','GFAP','RIMS2','SOX10'],ncols=4) 


# ### R100_61_MS_F_Age75

# In[ ]:


R100_61_MS_F_Age75_adata = filtering_DR_clustering(adata_dict['R100_61_MS_F_Age75'],'R100_61_MS_F_Age75')


# In[ ]:


sc.pl.dotplot(R100_61_MS_F_Age75_adata, Epithelial_marker, groupby='louvain_0.3', dendrogram=True,title='Epithelial')
sc.pl.dotplot(R100_61_MS_F_Age75_adata, Endothelial_marker, groupby='louvain_0.3', dendrogram=True,title='Endothelial')
sc.pl.dotplot(R100_61_MS_F_Age75_adata, Mesenchyma_marker, groupby='louvain_0.3', dendrogram=True,title='Mesenchymal')
sc.pl.dotplot(R100_61_MS_F_Age75_adata, Immune_marker, groupby='louvain_0.3', dendrogram=True,title='Immune')
sc.pl.dotplot(R100_61_MS_F_Age75_adata, T_cells_NK, groupby='louvain_0.3', dendrogram=True,title='T_cells')
sc.pl.dotplot(R100_61_MS_F_Age75_adata, Neuros+Astrocyte+Oligos, groupby='louvain_0.3', dendrogram=True,title='Neurons+Astro+Oligos')


# In[ ]:


R100_61_MS_F_Age75_adata.obs['CellType'] = R100_61_MS_F_Age75_adata.obs['louvain_0.3'].astype('category')
R100_61_MS_F_Age75_adata.obs['CellType'].replace({'0': 'Epithelial','1': 'Epithelial','2': 'Mesenchymal','3': 'Mesenchymal','4': 'Immune','5': 'NeuronGlial',
                                                  '6': 'NeuronGlial','7':'NeuronGlial','8':'Immune','9':'Endothelial'},
                                               inplace=True)

sc.pl.umap(R100_61_MS_F_Age75_adata,color='CellType',title='R100_61_MS_F_Age75')

R100_61_MS_F_Age75_adata.write_h5ad('../results/individual_samples/R100_61_MS_F_Age75_CellType.h5ad')


# In[ ]:


sc.pl.umap(R100_61_MS_F_Age75_adata,color=['CLIC6','OTX2','PECAM1','CLDN5','COL3A1','ACTA2','PTPRC','VCAN'],ncols=4) 


# ### R100_62_MS_F_Age62

# In[ ]:


R100_62_MS_F_Age62_adata = filtering_DR_clustering(adata_dict['R100_62_MS_F_Age62'],'R100_62_MS_F_Age62')


# In[ ]:


sc.pl.dotplot(R100_62_MS_F_Age62_adata, Epithelial_marker, groupby='louvain_0.4', dendrogram=True,title='Epithelial')
sc.pl.dotplot(R100_62_MS_F_Age62_adata, Endothelial_marker, groupby='louvain_0.4', dendrogram=True,title='Endothelial')
sc.pl.dotplot(R100_62_MS_F_Age62_adata, Mesenchyma_marker, groupby='louvain_0.4', dendrogram=True,title='Mesenchymal')
sc.pl.dotplot(R100_62_MS_F_Age62_adata, Immune_marker, groupby='louvain_0.4', dendrogram=True,title='Immune')
sc.pl.dotplot(R100_62_MS_F_Age62_adata, T_cells_NK, groupby='louvain_0.4', dendrogram=True,title='T_cells')
sc.pl.dotplot(R100_62_MS_F_Age62_adata, Neuros+Astrocyte+Oligos, groupby='louvain_0.4', dendrogram=True,title='Neurons+Astro+Oligos')

#sc.pl.dotplot(Spain_19_005_Ctrl_M_adata, Mesenchyma_marker+Neuros+Astrocyte+Oligos, groupby='louvain_0.2', dendrogram=True,title='Mesenchymal+Neurons+Astro+Oligos')


# In[ ]:


R100_62_MS_F_Age62_adata.obs['CellType'] = R100_62_MS_F_Age62_adata.obs['louvain_0.4'].astype('category')
R100_62_MS_F_Age62_adata.obs['CellType'].replace({'0': 'Epithelial','1': 'Epithelial','2': 'Epithelial','3': 'Immune','4': 'NeuronGlial','5': 'NeuronGlial',
                                                  '6': 'Mesenchymal','7':'Mesenchymal','8':'Mesenchymal','9':'Endothelial','10':'Immune'},
                                               inplace=True)

sc.pl.umap(R100_62_MS_F_Age62_adata,color='CellType',title='R100_62_MS_F_Age62')

R100_62_MS_F_Age62_adata.write_h5ad('../results/individual_samples/R100_62_MS_F_Age62_CellType.h5ad')


# In[ ]:


sc.pl.umap(R100_62_MS_F_Age62_adata,color=['CLIC6','OTX2','PECAM1','CLDN5','COL3A1','ACTA2','PTPRC','VCAN'],ncols=4) 


# ### R100_71_MS_F_Age52

# In[ ]:


R100_71_MS_F_Age52_adata = filtering_DR_clustering(adata_dict['R100_71_MS_F_Age52'],'R100_71_MS_F_Age52')


# In[ ]:


sc.pl.dotplot(R100_71_MS_F_Age52_adata, Epithelial_marker, groupby='louvain_0.2', dendrogram=True,title='Epithelial')
sc.pl.dotplot(R100_71_MS_F_Age52_adata, Endothelial_marker, groupby='louvain_0.2', dendrogram=True,title='Endothelial')
sc.pl.dotplot(R100_71_MS_F_Age52_adata, Mesenchyma_marker, groupby='louvain_0.2', dendrogram=True,title='Mesenchymal')
sc.pl.dotplot(R100_71_MS_F_Age52_adata, Immune_marker, groupby='louvain_0.2', dendrogram=True,title='Immune')
sc.pl.dotplot(R100_71_MS_F_Age52_adata, T_cells_NK, groupby='louvain_0.2', dendrogram=True,title='T_cells')
sc.pl.dotplot(R100_71_MS_F_Age52_adata, Neuros+Astrocyte+Oligos, groupby='louvain_0.2', dendrogram=True,title='Neurons+Astro+Oligos')


# In[ ]:


R100_71_MS_F_Age52_adata.obs['CellType'] = R100_71_MS_F_Age52_adata.obs['louvain_0.2'].astype('category')
R100_71_MS_F_Age52_adata.obs['CellType'].replace({'0': 'Mesenchymal',
                                               '1': 'Epithelial',
                                               '2': 'Mesenchymal',
                                               '3': 'Immune',
                                               '4': 'Immune',
                                               '5': 'NeuronGlial',
                                               '6': 'Endothelial',
                                               '7': 'Mesenchymal',
                                               '8':'NeuronGlial',
                                               '9':'NeuronGlial'},
                                               inplace=True)

sc.pl.umap(R100_71_MS_F_Age52_adata,color='CellType',title='R100_71_MS_F_Age52')

R100_71_MS_F_Age52_adata.write_h5ad('../results/individual_samples/R100_71_MS_F_Age52_CellType.h5ad')


# In[ ]:


sc.pl.umap(R100_71_MS_F_Age52_adata,color=['CLIC6','OTX2','PECAM1','CLDN5','COL3A1','ACTA2','PTPRC','VCAN'],ncols=4) 


# ## QC: after filtering

# In[ ]:


Spain_18_006_Ctrl_M_Age51_adata = sc.read_h5ad('../results/individual_samples/Spain_18_006_Ctrl_M_Age51_CellType.h5ad')
R100_63_Ctrl_F_Age56_adata = sc.read_h5ad('../results/individual_samples/R100_63_Ctrl_F_Age56_CellType.h5ad')
Spain_19_005_Ctrl_M_Age60_adata = sc.read_h5ad('../results/individual_samples/Spain_19_005_Ctrl_M_Age60_CellType.h5ad')
R100_64_Ctrl_M_Age64_adata = sc.read_h5ad('../results/individual_samples/R100_64_Ctrl_M_Age64_CellType.h5ad')
R77_61_Ctrl_F_Age79_adata = sc.read_h5ad('../results/individual_samples/R77_61_Ctrl_F_Age79_CellType.h5ad')
R100_72_Ctrl_F_Age82_adata = sc.read_h5ad('../results/individual_samples/R100_72_Ctrl_F_Age82_CellType.h5ad')
R77_60_MS_M_Age46_adata = sc.read_h5ad('../results/individual_samples/R77_60_MS_M_Age46_CellType.h5ad')
R100_71_MS_F_Age52_adata = sc.read_h5ad('../results/individual_samples/R100_71_MS_F_Age52_CellType.h5ad')
R100_62_MS_F_Age62_adata = sc.read_h5ad('../results/individual_samples/R100_62_MS_F_Age62_CellType.h5ad')
R100_61_MS_F_Age75_adata = sc.read_h5ad('../results/individual_samples/R100_61_MS_F_Age75_CellType.h5ad')

R108_30_Ctrl_M_Age88_adata = sc.read_h5ad('../results/individual_samples/R108_30_Ctrl_M_Age88_CellType.h5ad')
R108_31_Ctrl_F_Age82_adata = sc.read_h5ad('../results/individual_samples/R108_31_Ctrl_F_Age82_CellType.h5ad')


# In[ ]:


#merge filtered sample
adata_merge_filt = Spain_18_006_Ctrl_M_Age51_adata.concatenate(
    R100_63_Ctrl_F_Age56_adata,Spain_19_005_Ctrl_M_Age60_adata,R100_64_Ctrl_M_Age64_adata,R77_61_Ctrl_F_Age79_adata,R100_72_Ctrl_F_Age82_adata,R108_30_Ctrl_M_Age88_adata,R108_31_Ctrl_F_Age82_adata,
    R77_60_MS_M_Age46_adata,R100_71_MS_F_Age52_adata,R100_62_MS_F_Age62_adata,R100_61_MS_F_Age75_adata,
    
    batch_categories=['Spain_18_006_Ctrl_M_Age51','R100_63_Ctrl_F_Age56','Spain_19_005_Ctrl_M_Age60','R100_64_Ctrl_M_Age64','R77_61_Ctrl_F_Age79',
                      'R100_72_Ctrl_F_Age82','R108_30_Ctrl_M_Age88','R108_31_Ctrl_F_Age82','R77_60_MS_M_Age46','R100_71_MS_F_Age52','R100_62_MS_F_Age62','R100_61_MS_F_Age75'])

sc.pl.violin(adata_merge_filt,keys=['n_counts','n_genes','pct_counts_Mito','pct_counts_Ribo','pct_counts_Hb'],
                 rotation=90,stripplot=False,jitter=0.4,color = "#006837",multi_panel=True,groupby='batch')
sc.pl.violin(adata_merge_filt,keys=['S_score','G2M_score'],
             rotation=90,stripplot=False,jitter=0.4,color = "#006837",multi_panel=True,groupby='batch')


# ### The number of cells after filtering

# In[ ]:


ax = sns.countplot(y = 'batch',data = adata_merge_filt.obs)
#ax.set_ylabel('# Cells in each sample')
for container in ax.containers:
    ax.bar_label(container)


# ### Cell type across samples

# In[ ]:


def batch_pct_cellType(adata_file,groupby_key,legend_key,plot_title,pos):
    tmp = pd.crosstab(adata_file.obs[groupby_key],adata_file.obs[legend_key], normalize='index')   
    tmp.plot.bar(stacked=True,figsize=(10,10),title=plot_title, rot= 0, grid=False,
                 xlabel='',ax=pos).legend(bbox_to_anchor=[1.05,1])
    
fig, (ax1,ax2,ax3) = plt.subplots(3, 1, constrained_layout=True,figsize=(10,10))
batch_pct_cellType(adata_merge_filt,'CellType','batch','Cell summary before integration',ax1)
batch_pct_cellType(adata_merge_filt,'CellType','Condition','Cell summary before integration',ax2)
batch_pct_cellType(adata_merge_filt,'CellType','Gender','Cell summary before integration',ax3)


# ### Cell types

# In[ ]:


def cellType_plot(adata_file,count_key, pos,ylabel_title):
    import seaborn as sns
    df1 = adata_file.obs[['batch','Condition','Gender','CellType']]
    sns.countplot(y = count_key,data = df1,ax=pos).set_ylabel(ylabel_title)

fig, ax = plt.subplots(nrows=1,ncols=1,constrained_layout=True,figsize=(10,3))
cellType_plot(adata_merge_filt,'CellType',ax,'Cell types')
for container in ax.containers:
    ax.bar_label(container)


# In[ ]:


def batch_pct_cellType(adata_file,groupby_key,legend_key,plot_title,pos):
    tmp = pd.crosstab(adata_file.obs[groupby_key],adata_file.obs[legend_key], normalize='index')   
    tmp.plot.bar(stacked=True,figsize=(10,6),title=plot_title, rot= 90, grid=False,
                 xlabel='',ax=pos).legend(bbox_to_anchor=[1.05,1])

fig, ax = plt.subplots(1, 1, constrained_layout=True,figsize=(10,6))
batch_pct_cellType(adata_merge_filt,'batch','CellType','Cell summary before integration',ax)


# In[ ]:


import math
log10_count = adata_merge_filt.obs['n_counts'].apply(lambda x: math.log10(x))

fig, ax = plt.subplots(1, 1, constrained_layout=True,figsize=(24,6))
ax = sns.violinplot(x="CellType", y=log10_count, hue="batch",data=adata_merge_filt.obs)
ax.set_ylabel('log10(UMI)')

plt.legend(bbox_to_anchor=[1.01,1])
locs, labels = plt.xticks()
#plt.setp(labels, rotation=45)



# ### Remove genes expressed in less than 5 cells

# In[ ]:


sc.pp.filter_genes(adata_merge_filt,min_cells = 5)


# In[ ]:


adata_merge_filt.write_h5ad('../results/Human_filtered.h5ad')


# # Integrate samples
# Merge all samples and remove the batch effect.

# In[ ]:


#remove clustering infomation 
for name in adata_merge_filt.obs.columns:
    if name.startswith('louvain'):
        del adata_merge_filt.obs[name]


# Variable genes can be detected across the full dataset, but then we run the risk of getting many batch-specific genes that will drive a lot of the variation. Or we can select variable genes from each batch separately to get only celltype variation

# In[ ]:


human_raw = sc.read_h5ad('../results/Human_unfiltered.h5ad')


# In[ ]:


adata_merge_filt = sc.read_h5ad('../results/Human_filtered.h5ad')


# In[ ]:


keep_id = []
for i in human_raw.var_names:
    if i in adata_merge_filt.raw.var_names:
        keep_id.append(True)
    else:
        keep_id.append(False)
        


# In[ ]:


adata = human_raw[adata_merge_filt.obs.index]
#adata = human_raw[human_raw.obs.index.isin(adata_merge_filt.obs.index.tolist())]
adata= adata[:,keep_id]


# In[ ]:


adata.obs['S_score'] = adata_merge_filt.obs['S_score']
adata.obs['G2M_score'] = adata_merge_filt.obs['G2M_score']
adata.obs['CellType']=adata_merge_filt.obs['CellType']


# In[ ]:


adata.write_h5ad('../results/Human_filtered_rawCounts.h5ad')


# In[ ]:


del adata_merge_filt, human_raw


# In[ ]:


adata=sc.read_h5ad('../results/Human_filtered_rawCounts.h5ad')


# ## Normalization

# In[ ]:


# normalize to depth 10000
sc.pp.normalize_per_cell(adata,counts_per_cell_after=1e4)
#logaritmize
sc.pp.log1p(adata)
#store normalized counts in the raw slot, we will subset adata.X for variable genes, but want to keep all genes matrix as well
adata.raw = adata


# ## Select HVGs

# In[ ]:


sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
var_genes_all = adata.var.highly_variable
print("Highly variable genes: %d"%sum(var_genes_all))

#Detect variable genes in each dataset separately using the batch_key parameter.
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key = 'batch')
print("Highly variable genes intersection: %d"%sum(adata.var.highly_variable_intersection))

print("Number of batches where gene is variable:")
print(adata.var.highly_variable_nbatches.value_counts())

var_genes_batch = adata.var.highly_variable_nbatches > 0

#Compare overlap of the variable genes.

print("Any batch var genes: %d"%sum(var_genes_batch))
print("All data var genes: %d"%sum(var_genes_all))
print("Overlap: %d"%sum(var_genes_batch & var_genes_all))
print("Variable genes in all batches: %d"%sum(adata.var.highly_variable_nbatches == 12))
print("Overlap batch instersection and all: %d"%sum(var_genes_all & adata.var.highly_variable_intersection))

#Select all genes that are variable in at least 4 datasets and use for remaining analysis.
var_select = adata.var.highly_variable_nbatches > 4
var_genes = var_select.index[var_select]
len(var_genes)


# ## DR

# In[ ]:


adata2 = adata 
#store lognormalized matrix in adata.raw slot
adata2.raw=adata2
#subset for variable genes in the dataset
adata2 = adata2[:,var_genes]

sc.pp.regress_out(adata2,['percent_chrY'])
# scale data, clip values exceeding standard deviation 10.
sc.pp.scale(adata2,max_value=10)
sc.tl.pca(adata2,svd_solver = 'arpack',n_comps = 50)
sc.pl.pca_variance_ratio(adata2,log = True,n_pcs = 50)

#umap
sc.pp.neighbors(adata2,n_pcs = 30)
sc.tl.umap(adata2)
sc.tl.tsne(adata2,n_pcs=30,random_state=0)


# In[ ]:


adata2.write_h5ad('../results/Human_withBatchEffect.h5ad')


# In[ ]:


sc.pl.umap(adata2,color=['batch','CellType','Gender','Condition'],ncols=2,wspace=0.7)


# ## Harmony Integration

# In[ ]:


import harmonypy
import scanpy.external as sce
adata_har = adata2
sce.pp.harmony_integrate(adata_har, 'batch')
# tsne and umap
sc.pp.neighbors(adata_har, n_pcs =50, use_rep = "X_pca_harmony")
sc.tl.umap(adata_har)
sc.tl.tsne(adata_har, n_pcs = 50, use_rep = "X_pca_harmony")


# In[ ]:


adata_har.write_h5ad('../results/Human_harmony.h5ad')


# ## Integration comparison

# In[ ]:


adata2 = sc.read_h5ad('../results/Human_withBatchEffect.h5ad')
fig, (ax1,ax2) = plt.subplots(2, 1, constrained_layout=True, figsize=(8,10),gridspec_kw={'wspace':0.1,'hspace': 0.1})
sc.pl.umap(adata2,color='batch',title='Without integration',ax=ax1,show=False)
sc.pl.umap(adata_har,color='batch',title='Harmony integration',ax=ax2,show=False)


# In[ ]:


del adata2,adata_har
gc.collect()


# # Clustering

# In[ ]:


adata = sc.read_h5ad('../results/Human_harmony.h5ad')


# In[ ]:


sc.pl.umap(adata,color=['Condition','CellType','Gender','batch'],wspace=0.5,ncols=2)


# In[ ]:


res = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.0]
for snn_res in res:
    sc.tl.louvain(adata, resolution = snn_res, key_added = 'louvain_%s' %snn_res) 


# In[ ]:


#plot
sc.pl.umap(adata,color=['louvain_0.1','louvain_0.2','louvain_0.3','louvain_0.4',
                        'louvain_0.5','louvain_0.6','louvain_0.7','louvain_0.8','louvain_1.0'],
           title=['louvain_0.1','louvain_0.2','louvain_0.3','louvain_0.4','louvain_0.5',
                  'louvain_0.6','louvain_0.7','louvain_0.8','louvain_1.0'],
           ncols=3, wspace=0.5)


# ## Top markers

# In[ ]:


del adata.uns['log1p']


# In[ ]:


res_list = adata.obs.columns[adata.obs.columns.str.startswith('louvain')]
for snn_res in res_list:
    sc.tl.rank_genes_groups(adata,snn_res,method='wilcoxon',key_added='wilcoxon_%s' %snn_res,pts=True) 


# In[ ]:


adata.write_h5ad('../results/Human_harmony_clustering.h5ad')


# In[ ]:


for snn_res in res_list:
    sc.tl.dendrogram(adata, use_rep='X_pca_harmony',groupby=snn_res)
    sc.pl.rank_genes_groups_dotplot(adata,n_genes=10,key='wilcoxon_%s' %snn_res, groupby=snn_res,show=True)


# ### Save top markers

# In[ ]:


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


# In[ ]:


wilcoxon_key_list=[]
for key in adata.uns.keys():
    if key.startswith('wilcoxon_louvain_'):
        wilcoxon_key_list.append(key)

for i in wilcoxon_key_list:
    df = rank_genes_groups_df(adata,i)
    df = df[df['pvals_adj']<0.05]
    
    #only keep the top 100 genes for each cluster. Genes were sorted by the Score. 
    new_df = pd.DataFrame([], columns = df.columns)
    for index in np.unique(df.index):
        new_df = pd.concat([new_df, df.loc[index].head(n=100)])
    
    new_df.to_csv('../results/harmony_markers/%s_harmony_markers.csv'%i)
    
    del df, new_df


# In[ ]:


adata.write_h5ad('../results/Human_harmony_CellType.h5ad')


# ## Filtered top markers

# In[ ]:


adata = sc.read_h5ad('../results/Human_harmony_CellType.h5ad')


# In[ ]:


res_list = adata.obs.columns[adata.obs.columns.str.startswith('louvain')]
for snn_res in res_list:
    sc.tl.filter_rank_genes_groups(adata, groupby=snn_res, key = "wilcoxon_%s" %snn_res, key_added='wilcoxon_filtered_%s' %snn_res,
                                   min_in_group_fraction=0.1,max_out_group_fraction=0.5,min_fold_change=0.25)


# ### Visualization

# In[ ]:


for snn_res in res_list:
    sc.tl.dendrogram(adata, use_rep='X_pca_harmony',groupby=snn_res)
    sc.pl.rank_genes_groups_dotplot(adata,n_genes=10,key='wilcoxon_filtered_%s' %snn_res, groupby=snn_res,show=True)


# In[ ]:


wilcoxon_key_list=[]
for key in adata.uns.keys():
    if key.startswith('wilcoxon_filtered_'):
        wilcoxon_key_list.append(key)

for key in wilcoxon_key_list:
    df = rank_genes_groups_df(adata,i)
    df = df[(df['pvals_adj']<0.05)]
    df = df[~df['names'].isnull()]
    
    #only keep the top 100 genes for each cluster. Genes were sorted by the Score. 
    new_df = pd.DataFrame([], columns = df.columns)
    for index in np.unique(df.index):
        new_df = pd.concat([new_df, df.loc[index].head(n=100)])
    
    new_df.to_csv('../results/harmony_markers/%s_harmony_markers.csv' %key)
    
    del df, new_df


# ## Features plot

# In[ ]:


sc.pl.umap(adata,color=['louvain_0.6','CellType'],ncols=2,wspace=0.5)


# In[ ]:


# canonical markers expression
sc.pl.umap(adata,color=['CLIC6','OTX2','PECAM1','CLDN5','COL3A1','ACTA2','PTPRC','VCAN','CNP','GFAP','RIMS2','SOX10'],ncols=4) 


# # Cell Types

# In[ ]:


sc.pl.dotplot(adata, Epithelial_marker, groupby='louvain_0.4', dendrogram=True,title='Epithelial')
sc.pl.dotplot(adata, Endothelial_marker, groupby='louvain_0.4', dendrogram=True,title='Endothelial')
sc.pl.dotplot(adata, Mesenchyma_marker, groupby='louvain_0.4', dendrogram=True,title='Mesenchymal')
sc.pl.dotplot(adata, Immune_marker, groupby='louvain_0.4', dendrogram=True,title='Immune')
sc.pl.dotplot(adata, T_cells_NK, groupby='louvain_0.4', dendrogram=True,title='T_cells')
sc.pl.dotplot(adata, Neuros+Astrocyte+Oligos, groupby='louvain_0.4', dendrogram=True,title='Neurons+Astro+Oligos')

sc.pl.dotplot(adata, Epithelial_marker+Endothelial_marker+Mesenchyma_marker+Immune_marker+T_cells_NK+Neuros+Astrocyte+Oligos, 
              groupby='louvain_0.4', dendrogram=True)


# ## NeuronGlia

# In[ ]:


#Astrocyte
sc.pl.umap(adata,color=['GFAP','ALDH1L1','louvain_0.4'],ncols=3) # Cluster5 is astrocytes


# In[ ]:


#OPC
sc.pl.umap(adata,color=['PDGFRA','CSPG4','SOX10'])# Cluster4 is OPC


# In[ ]:


#Oligos
sc.pl.umap(adata,color=['FA2H','CNP','UGT8'])


# In[ ]:


# Neurons
sc.pl.umap(adata,color=["KCNIP4","NRXN1","KCND2","RIMS1","RIMS2","DSCAM","TNR","VCAN"],ncols=4)


# In[ ]:


# Ependymal
sc.pl.umap(adata,color=['DNAH9','DCDC1','ADGB'],ncols=3)


# ## Check the DGE between two neuron clusters

# In[ ]:


sc.tl.rank_genes_groups(adata,'louvain_0.4',groups=['10','13'],method='wilcoxon',key_added='wilcoxon_neurons')


# In[ ]:


sc.tl.filter_rank_genes_groups(adata, groupby='louvain_0.4', key = "wilcoxon_neurons", key_added='wilcoxon_filtered_neurons',
                               min_in_group_fraction=0.1,max_out_group_fraction=0.5,min_fold_change=1)
sc.pl.rank_genes_groups_dotplot(adata,n_genes=25,key='wilcoxon_filtered_neurons', groupby='louvain_0.4',show=True)


# In[ ]:


adata.obs['CellType_intergrate'] = adata.obs['louvain_0.4'].astype('category')
adata.obs['CellType_intergrate'].replace({'0': 'Epithelial','1': 'Mesenchymal','2': 'Immune','3': 'Mesenchymal','4': 'OPC+OLG','5': 'Astrocyte',
                               '6': 'Immune','7': 'Mesenchymal','8':'Endothelial','9':'Mesenchymal','10':'Neuron1','11':'Epithelial','12':'Epithelial','13':'Neuron2'},
                              inplace=True)

sc.pl.umap(adata,color=['CellType','CellType_intergrate','louvain_0.4'],title='CellType',wspace=0.5)


# In[ ]:


#adata.obs['CellType_intergrate'] = adata.obs['CellType_intergrate'].cat.add_categories(['Ependymal'])
for cell in adata.obs.index[(adata.obs['louvain_0.6']=='14')]:
    adata.obs.loc[cell,'CellType_intergrate']='Ependymal'


# In[ ]:


sc.pl.umap(adata,color='CellType_intergrate')


# # Save data

# In[ ]:


#remove redundant info
wilcoxon_key_list=[]
for key in adata.uns.keys():
    if key.startswith('wilcoxon_filtered_'):
        wilcoxon_key_list.append(key)

for i in wilcoxon_key_list:
    del adata.uns[i]
    
adata.write_h5ad('../results/Human_harmony_CellType.h5ad')


# ## Summarize cells

# In[ ]:


ax = sns.countplot(y = 'CellType_intergrate',data = adata.obs)
#ax.set_ylabel('# Cells in each sample')
for container in ax.containers:
    ax.bar_label(container)


# In[ ]:


def batch_pct_cellType(adata_file,groupby_key,legend_key,pos):
    tmp = pd.crosstab(adata_file.obs[groupby_key],adata_file.obs[legend_key], normalize='index')   
    tmp.plot.bar(stacked=True,rot= 90, grid=False,
                 xlabel='',ax=pos).legend(bbox_to_anchor=[1.05,1])


# In[ ]:


fig, (ax1,ax2,ax3,ax4) = plt.subplots(4, 1, constrained_layout=True,figsize=(8,16))
batch_pct_cellType(adata,'batch','CellType_intergrate',ax1)
batch_pct_cellType(adata,'Condition','CellType_intergrate',ax2)
batch_pct_cellType(adata,'Gender','CellType_intergrate',ax3)
batch_pct_cellType(adata,'CellType_intergrate','batch',ax4)


# ## QC, mito, cell cycle and gender check

# In[ ]:


sc.pl.umap(adata,color=['CellType_intergrate','pct_counts_Mito','pct_counts_Ribo','n_counts','S_score','G2M_score','percent_chrY'],ncols=4,wspace=0.5)


# In[ ]:




