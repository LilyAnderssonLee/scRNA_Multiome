"""
  save highly variable, pca, harmony, tsne and umap matrix into csv file in python
  
  lv is h5ad file
  lv = sc.read_h5ad('../results/lv_harmony_rmBatchbetweenSamples_CellType_noDBL_edit1_noDBL.h5ad')
  hvg_dict = {'highly_variable':lv.var_names[lv.var['highly_variable']].tolist()}
  df = pd.DataFrame(hvg_dict)
  df.to_csv('../cell_cell_communication/lv_HVGs.csv')
  pd.DataFrame(lv.obsm['X_pca']).to_csv('../cell_cell_communication/lv_X_pca.csv')
  pd.DataFrame(lv.obsm['X_pca_harmony']).to_csv('../cell_cell_communication/lv_X_pca_harmony.csv')
  pd.DataFrame(lv.obsm['X_tsne']).to_csv('../cell_cell_communication/lv_X_tsne.csv')
  pd.DataFrame(lv.obsm['X_umap']).to_csv('../cell_cell_communication/lv_X_umap.csv')
"""

library(Seurat)
setwd('OneDrive - Karolinska Institutet/Mac/Documents/snRNA/MS/mouse/ms_model_LV_SLV/cellChat/input/')
seu<-readRDS('lv_cellchat_input.rds')
hvgs <- read.csv('lv_HVGs.csv')
pca <- read.csv('lv_X_pca.csv')
harmony <- read.csv('lv_X_pca_harmony.csv')  
tsne <- read.csv('lv_X_tsne.csv')
umap <- read.csv('lv_X_umap.csv')

#add DR matrix
##pca
pca <- pca[,2:51]
colnames(pca) <- paste0("PC_", 1:50)
rownames(pca) <- colnames(seu)
seu[['X_pca']]<-CreateDimReducObject(embeddings = as.matrix(pca),key = 'PC_',assay='RNA')
##harmony
harmony <- harmony[,2:51]
colnames(harmony) <- paste0("harmony_", 1:50)
rownames(harmony) <- colnames(seu)
seu[['X_pca_harmony']]<-CreateDimReducObject(embeddings = as.matrix(harmony),key = 'harmony_',assay='RNA')
##tsne
tsne <- tsne[,2:3]
colnames(tsne) <- paste0("tsne_", 1:2)
rownames(tsne) <- colnames(seu)
seu[['X_tsne']]<-CreateDimReducObject(embeddings = as.matrix(tsne),key = 'tsne_',assay='RNA')
##umap
umap <- umap[,2:3]
colnames(umap) <- paste0("UMAP_", 1:2)
rownames(umap) <- colnames(seu)
seu[['X_umap']]<-CreateDimReducObject(embeddings = as.matrix(umap),key = 'UMAP_',assay='RNA')
##HVGs
seu@assays$RNA@var.features <- hvgs$highly_variable

seu$CellType<-NULL

seu$Condition <- gsub('MS','DTA',seu$Condition)
saveRDS(seu,'lv_h5ad_info.rds')

#CellType_correct_edit1 only contains four cell types: Endothelial, Epithelial, Immune and Mesenchymal
#CellType_correct contains six cell types: Endothelial, Epithelial, Lymphocyte, Macrophage, Mesenchymal and Pericyte
#Subcell includes all subcell infomation.




