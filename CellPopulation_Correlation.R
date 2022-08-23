#########################################cell populations correlations#################
#load data
setwd("input/")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(pheatmap)
  library(Matrix)
  library(cowplot)
  library(venn)
}
)

target<-readRDS("human_hvgs5000.rds")
ref<-readRDS("ref_hcgs5000.rds") # reference

#different metods to estimate the average of each cell group
triMean <- function(x, na.rm = TRUE) {
  mean(stats::quantile(x, probs = c(0.25, 0.50, 0.50, 0.75), na.rm = na.rm))
}

truncatedMean = function(x) mean(x, trim = trim, na.rm = TRUE)

#ommit genes with less than 10% expressed cells with each group
thresholdedMean <- function(x, trim = 0.1, na.rm = TRUE) {
  percent <- Matrix::nnzero(x)/length(x)
  if (percent < trim) {
    return(0)
  } else {
    return(mean(x, na.rm = na.rm))
  }
}

#the following estimation will based on thresholdedMean 

ref_group <- ref@active.ident
#define the top features
ref_hvgs<-ref@assays$RNA@var.features
target_hvgs<-target@assays$RNA@var.features
comm_feature<-intersect(ref_hvgs,target_hvgs)

ref_hvgs<-ref@assays$RNA@data[comm_feature,]
#average expression of each cell group
ref_hvgs_cutoffMean<-aggregate(t(ref_hvgs), list(ref_group), FUN = thresholdedMean)
ref_hvgs_cutoffMean_Matrix<-data.matrix(ref_hvgs_cutoffMean)
ref_hvgs_cutoffMean_Matrix<-t(ref_hvgs_cutoffMean_Matrix)
ref_hvgs_cutoffMean_Matrix<-ref_hvgs_cutoffMean_Matrix[-1,]
colnames(ref_hvgs_cutoffMean_Matrix)<-ref_hvgs_cutoffMean$Group.1

######target data
hvgs<-target@assays$RNA@data[comm_feature,]
#average expression of each cell group
#option1
hvgs_cutoffMean<-aggregate(t(hvgs), list(target@active.ident), FUN = thresholdedMean)
#clusters based on PCs15_CCA_snn_res.0.5 
#option2
hvgs_cutoffMean<-aggregate(t(hvgs), list(target$PCs15_CCA_snn_res.0.5), FUN = thresholdedMean)

hvgs_cutoffMean_Matrix<-data.matrix(hvgs_cutoffMean)
hvgs_cutoffMean_Matrix<-t(hvgs_cutoffMean_Matrix)
hvgs_cutoffMean_Matrix<-hvgs_cutoffMean_Matrix[-1,]
colnames(hvgs_cutoffMean_Matrix)<-hvgs_cutoffMean$Group.1


x<-cor.test(x=hvgs_cutoffMean_Matrix[,1],y=ref_hvgs_cutoffMean_Matrix[,1])

for (i in 1:dim(hvgs_cutoffMean_Matrix)[2]) {
  p_value<-list()
  cor_coeff<-list()
  apply(ref_hvgs_cutoffMean_Matrix,2,function(x){
    cor_coeff[[i]]<-cor(hvgs_cutoffMean_Matrix[,i],ref_hvgs_cutoffMean_Matrix[,1])
    #p_value[[i]]<-cor.test(hvgs_cutoffMean_Matrix[,i],x)$p.value
    #cor_coeff[[i]]<-cor.test(hvgs_cutoffMean_Matrix[,i],x)$estimate
  })
  return(cor_coeff)
}

#estimate pearson correlations 
cor_coeff<-list()
p_value<-list()
for (i in 1:dim(hvgs_cutoffMean_Matrix)[2]) {
  query_cellType<-hvgs_cutoffMean_Matrix[,i]
  
  cor_coeff[[i]]<-NaN*seq(dim(ref_hvgs_cutoffMean_Matrix)[2])
  p_value[[i]]<-NaN*seq(dim(ref_hvgs_cutoffMean_Matrix)[2])
  
  for (j in 1:dim(ref_hvgs_cutoffMean_Matrix)[2]) {
    est_cor<-cor.test(query_cellType,ref_hvgs_cutoffMean_Matrix[,j],method = "pearson")
    cor_coeff[[i]][j]<-est_cor$estimate
    p_value[[i]][j]<-est_cor$p.value
  }
}

names(cor_coeff)<-colnames(hvgs_cutoffMean_Matrix)
names(p_value)<-names(cor_coeff)

cor_df<-as.matrix(do.call(cbind,cor_coeff))
p_value_df<-as.matrix(do.call(cbind,p_value))
rownames(cor_df)<-colnames(ref_hvgs_cutoffMean_Matrix)
rownames(p_value_df)<-rownames(cor_df)

#heatmap

# heatmap.2 function with manhattan method and ward.D2 clustering
library(gplots)
colnames(cor_df)<-paste0("Cluster_",colnames(cor_df))
heatmap.2(x=cor_df,margins = c(10,10),
          hclustfun = function(x) hclust(x, method="ward.D2"),
          distfun = function(x) dist(x,method = "manhattan")
          )
#heatmap(cor_df,margins = c(10,10))
#ggplot
#formate
melt_cor_df<-melt(cor_df)
melt_p_value_df<-melt(p_value_df)
melt_cor_df$p.value<-melt_p_value_df$value
#significant correlations
melt_cor_df$r_sig<-ifelse(melt_cor_df$p.value<0.05,melt_cor_df$value,NA)

#order cell groups of single nucleir RNA seq from ref 
level_order<-c("Oligo6","Oligo4","Oligo3","Oligo2","Oligo1","Oligo5","Astrocytes2","COPs","ImOlGs",
               "Neuron4","Neuron3","Neuron1","Neuron2","OPCs","Neuron5","Astrocytes",
               "Endothelial_cells1","Endothelial_cells2","Pericytes","Vasc_smooth_muscle",
               "Microglia_Macrophages","Macrophages","Immune_cells")
melt_cor_df$Var2<-paste0("Cluster_",melt_cor_df$Var2)
ggplot(melt_cor_df,aes(x=Var2,y=factor(Var1,levels = level_order),fill=value,label=round(r_sig,2)))+
  geom_tile()+geom_text()+theme_classic()+
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  #scale_fill_gradient2(high = "red",low = "blue",mid="white")+
  labs(x = NULL, y = NULL, fill = "Pearson's\nCorrelation", title="Cell groups Correlations", subtitle="Only significant Pearson's correlation coefficients shown") + 
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-0.1,1))




