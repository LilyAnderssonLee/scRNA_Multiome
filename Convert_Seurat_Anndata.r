library(SeuratDisk)
library(Seurat)
library(loomR)
library(biomaRt)

Convert("../input/mouse_epith_raw.h5ad", dest = "h5seurat",overwrite = TRUE)

seu<-LoadH5Seurat('../input/mouse_epith_raw.h5seurat')
DefaultAssay(seu)<-'RNA'
saveRDS(seu,'../input/mouse_epith_raw.rds')
