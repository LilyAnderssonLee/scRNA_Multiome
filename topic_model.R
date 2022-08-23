#! /usr/bin/Rscript --vanilla --default-packages=utils
args = commandArgs(trailingOnly=TRUE)

library(Seurat)
library(CountClust)
library(Matrix)
library(dplyr)

setwd('/proj/xxx/private/Lili/analysis/ms_scRNA/mouse/pipeline_withoutMT/subcell_analysis/subcell_withDoublets/topicmodel/')
#Top model use raw counts
seu<-readRDS('input/Epithelial.rds')

#Remove genes expressed in more than 98% or less than 2% cells
high_pct<-dim(seu)[1]*0.98
low_pct<-dim(seu)[1]*0.02
features_keep<-rownames(seu)[Matrix::rowSums(seu)>low_pct & Matrix::rowSums(seu) <high_pct]
seu.filt<-subset(seu,features = features_keep)
saveRDS(seu.filt,"input/Epithelial_filt_TM.rds")

#Top model
topicModeling <- function(mat, n = 10, tol = 0.1, save.path) {
  n.topics <- as.numeric(n)
  tolerance <- as.numeric(tol)

  # create directories to save the results in
  sub.dir <- paste0("Epithelial/CountClust/",n.topics, "topics_tol", tolerance)
  dir.create(file.path(save.path, sub.dir), recursive = T)

  FitGoM(t(as.matrix(GetAssayData(mat))), K = n.topics, tol = tolerance,
         path_rda = file.path(save.path, sub.dir, paste0('FitGoM_k', n.topics, '_tol', tolerance, '.rda')))
}

out.dir <- file.path("TM_output")
if (!dir.exists(out.dir)){
  dir.create(out.dir,recursive = T)
} else {
  print("Dir already exists!")
}

tol <- 0.1

k=args[1]

topicModeling(seu.filt, k, tol, out.dir)


