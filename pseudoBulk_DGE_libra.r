library(DESeq2)
library(edgeR)
library(limma)
library(devtools)
library(Libra)
library(Seurat)
library(ggplot2)
library(ggrepel)

Spain_18_006_Ctrl_M_Age51 = Read10X_h5('../../../../raw_counts_scRNA/spain_human/AU7003-18-006-SI-GA-C12/outs/filtered_feature_bc_matrix.h5',use.names = T)
R100_63_Ctrl_F_Age56 = Read10X_h5('../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_063_new/outs/filtered_feature_bc_matrix.h5',use.names = T)
Spain_19_005_Ctrl_M_Age60 = Read10X_h5('../../../../raw_counts_scRNA/spain_human/AU7006-19-005-SI-GA-D3/outs/filtered_feature_bc_matrix.h5',use.names = T)
R100_64_Ctrl_M_Age64 = Read10X_h5('../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_064_new/outs/filtered_feature_bc_matrix.h5',use.names = T)
R77_61_Ctrl_F_Age79 = Read10X_h5('../../../../raw_counts_scRNA/10X_20_AF_09_R77/10X_20_061_premrna/outs/filtered_feature_bc_matrix.h5',use.names = T)
R100_72_Ctrl_F_Age82 = Read10X_h5('../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_072_new/outs/filtered_feature_bc_matrix.h5',use.names = T)

R77_60_MS_M_Age46 = Read10X_h5('../../../../raw_counts_scRNA/10X_20_AF_09_R77/10X_20_060_premrna/outs/filtered_feature_bc_matrix.h5',use.names = T)
R100_71_MS_F_Age52 = Read10X_h5('../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_071_new/outs/filtered_feature_bc_matrix.h5',use.names = T)
R100_62_MS_F_Age62 = Read10X_h5('../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_062_new/outs/filtered_feature_bc_matrix.h5',use.names = T)
R100_61_MS_F_Age75 = Read10X_h5('../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_061_new/outs/filtered_feature_bc_matrix.h5',use.names = T)

R108_30_Ctrl_M_Age88 = Read10X_h5('../../../../../raw_data/10X_22_AF_05_R108/10X_22_030/outs/filtered_feature_bc_matrix.h5',use.names = T)
R108_31_Ctrl_F_Age82= Read10X_h5('../../../../../raw_data/10X_22_AF_05_R108/10X_22_031/outs/filtered_feature_bc_matrix.h5',use.names = T)

sdata.Spain_18_006_Ctrl_M_Age51 <- CreateSeuratObject(Spain_18_006_Ctrl_M_Age51, project = "Spain_18_006_Ctrl_M_Age51")
sdata.R100_63_Ctrl_F_Age56 <- CreateSeuratObject(R100_63_Ctrl_F_Age56, project = "R100_63_Ctrl_F_Age56")
sdata.Spain_19_005_Ctrl_M_Age60 <- CreateSeuratObject(Spain_19_005_Ctrl_M_Age60, project = "Spain_19_005_Ctrl_M_Age60")
sdata.R100_64_Ctrl_M_Age64 <- CreateSeuratObject(R100_64_Ctrl_M_Age64, project = "R100_64_Ctrl_M_Age64")
sdata.R77_61_Ctrl_F_Age79 <- CreateSeuratObject(R77_61_Ctrl_F_Age79, project = "R77_61_Ctrl_F_Age79")
sdata.R100_72_Ctrl_F_Age82 <- CreateSeuratObject(R100_72_Ctrl_F_Age82, project = "R100_72_Ctrl_F_Age82")
sdata.R77_60_MS_M_Age46 <- CreateSeuratObject(R77_60_MS_M_Age46, project = "R77_60_MS_M_Age46")
sdata.R100_71_MS_F_Age52 <- CreateSeuratObject(R100_71_MS_F_Age52, project = "R100_71_MS_F_Age52")
sdata.R100_62_MS_F_Age62 <- CreateSeuratObject(R100_62_MS_F_Age62, project = "R100_62_MS_F_Age62")
sdata.R100_61_MS_F_Age75 <- CreateSeuratObject(R100_61_MS_F_Age75, project = "R100_61_MS_F_Age75")
sdata.R108_30_Ctrl_M_Age88 <- CreateSeuratObject(R108_30_Ctrl_M_Age88, project = "R108_30_Ctrl_M_Age88")
sdata.R108_31_Ctrl_F_Age82 <- CreateSeuratObject(R108_31_Ctrl_F_Age82, project = "R108_31_Ctrl_F_Age82")

# Merge datasets into one single seurat object
alldata <- merge(sdata.Spain_18_006_Ctrl_M_Age51, c(sdata.R100_63_Ctrl_F_Age56,sdata.Spain_19_005_Ctrl_M_Age60 ,sdata.R100_64_Ctrl_M_Age64 ,sdata.R77_61_Ctrl_F_Age79,
                                                   sdata.R100_72_Ctrl_F_Age82,sdata.R108_30_Ctrl_M_Age88,sdata.R108_31_Ctrl_F_Age82,
                                                    sdata.R77_60_MS_M_Age46,sdata.R100_71_MS_F_Age52,sdata.R100_62_MS_F_Age62,sdata.R100_61_MS_F_Age75),
                 add.cell.ids = c("Spain_18_006_Ctrl_M_Age51", "R100_63_Ctrl_F_Age56", "Spain_19_005_Ctrl_M_Age60", "R100_64_Ctrl_M_Age64",
                                  "R77_61_Ctrl_F_Age79", "R100_72_Ctrl_F_Age82",'R108_30_Ctrl_M_Age88','R108_31_Ctrl_F_Age82',
                                  'R77_60_MS_M_Age46','R100_71_MS_F_Age52','R100_62_MS_F_Age62','R100_61_MS_F_Age75'))

meta = read.csv('old/Human_harmony_CellType_edit_noDBL.csv')

#alldata cell_id format is different from scanpy cell_id,so we need to change the cell_id format in alldata into scanpy cell_id
library(stringr)
print(str_length('AAACCCAAGACTCTTG-1'))
umi = substring(meta$X,first = 1,last = 18)
rownames(meta)<-paste(meta$batch,umi,sep = '_')

meta$X<-NULL

gene_name=read.csv('old/gene_name.csv')

alldata.filt <- subset(alldata, cells = rownames(meta),features = gene_name$gene_name)

# Filter MALAT1
#alldata.filt <- alldata.filt[!grepl("MALAT1", rownames(alldata.filt)), ]

# Filter Mitocondrial
#alldata.filt <- alldata.filt[!grepl("^MT-", rownames(alldata.filt)), ]

merged_meta = merge(alldata.filt@meta.data,meta,by='row.names',all=TRUE)
rownames(merged_meta)<-merged_meta$Row.names
merged_meta$Row.names<-NULL


#check if the cell_id order of meta (scanpy) is same as alldata.filt
sum(rownames(meta)==rownames(alldata.filt@meta.data))

alldata.filt$cell_type = meta$CellType
alldata.filt$NeuroGlia_subcell_Ana = meta$NeuroGlia_subcell_Ana
alldata.filt$Immune_subcell_Ana = meta$Immune_subcell_Ana
alldata.filt$Epithelial_subcell_Ana = meta$Epithelial_subcell_Ana
alldata.filt$Mesenchymal_subcell_Ana = meta$Mesenchymal_subcell_Ana
alldata.filt$Endothelial_Subcells = meta$Endothelial_Subcells
alldata.filt$label = meta$Condition
alldata.filt$Age = meta$Age
alldata.filt$Gender = meta$Gender
alldata.filt$percent_chrY = meta$percent_chrY
alldata.filt$S_score = meta$S_score
alldata.filt$G2M_score = meta$G2M_score
alldata.filt$pct_counts_Mito = meta$pct_counts_Mito
alldata.filt$pct_counts_Ribo = meta$pct_counts_Ribo
alldata.filt$pct_counts_Hb = meta$pct_counts_Hb
alldata.filt$replicate = meta$batch

alldata.filt$label = as.factor(alldata.filt$label)
levels(alldata.filt$label) <- c('Control','MS')

# save data

saveRDS(alldata.filt,'old/Libra_alldata_filt.rds')

alldata.filt = readRDS('old/Libra/Libra_alldata_filt.rds')
levels(alldata.filt$label) <- c('Control','MS')

alldata.filt$project='Human'

matrices = to_pseudobulk(alldata.filt,min_cells = 5,min_reps = 0,cell_type_col = "project")
DE = run_de(alldata.filt,cell_type_col = "project")

write.csv(DE,'old/Libra/Libra_edgeR_LRT_Human_allcells_DE.csv')
highlight_DE = DE[(abs(DE$avg_logFC)>2 & DE$p_val_adj<0.05),]
write.csv(highlight_DE,paste0('old/Libra/Libra_edgeR_LRT_Human_allcells_sig_DGEs.csv'))
print(highlight_DE,n=2000)
#plot
options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 300,ggrepel.max.overlaps = Inf)
print('Control vs MS: negative avg logFC means down-regulated in Control and positive means up-regulated in Control.')

ggplot(data = DE, aes(x=avg_logFC, y=-log10(p_val_adj)),colour=gene[abs(avg_logFC>2)])+ geom_point(pch=20,size=0.8,color='#cccccc')+geom_point(data=highlight_DE,aes(x=avg_logFC,y=-log10(p_val_adj)),color='#0571b0',size=0.8,pch=20)+
    facet_wrap(~cell_type,ncol = 3)+theme_bw() + theme(legend.position="none")+
    geom_vline(xintercept = 2,size=0.3,color='black') + geom_vline(xintercept = -2,size=0.3,color='black') + geom_hline(yintercept = -log10(0.05),size=0.3,color='black') + xlab("avg logFC") + ylab("-log10(padj")+
    geom_text_repel(aes(label=ifelse(p_val_adj<0.05 & abs(avg_logFC)>2,gene,'')),size=1,color='#ca0020',fontface = 'bold', box.padding = unit(0.35, "lines"),point.padding = unit(0.5, "lines"),segment.color = 'grey50')

matrices = to_pseudobulk(alldata.filt,min_cells = 5,min_reps = 0)

#saveRDS(matrices,'old/libra_pseudoBulk_matrix.rds')

DE = run_de(alldata.filt)

write.csv(DE,'old/Libra/Libra_edgeR_LRT_allcells_DE.csv')

highlight_DE = DE[(abs(DE$avg_logFC)>2 & DE$p_val_adj<0.05),]
write.csv(highlight_DE,paste0('old/Libra/Libra_edgeR_LRT_allcells_sig_DGEs.csv'))

print(highlight_DE,n=2000)

options(rePercentageFeatureSet.width = 14, repr.plot.height = 6, repr.plot.res = 300,ggrepel.max.overlaps = Inf)
print('Control vs MS: negative avg logFC means down-regulated in Control and positive means up-regulated in Control.')

ggplot(data = DE, aes(x=avg_logFC, y=-log10(p_val_adj)),colour=gene[abs(avg_logFC>2)])+ geom_point(pch=20,size=0.8,color='#cccccc')+geom_point(data=highlight_DE,aes(x=avg_logFC,y=-log10(p_val_adj)),color='#0571b0',size=0.8,pch=20)+
    facet_wrap(~cell_type,ncol = 3)+theme_bw() + theme(legend.position="none")+
    geom_vline(xintercept = 2,size=0.3,color='black') + geom_vline(xintercept = -2,size=0.3,color='black') + geom_hline(yintercept = -log10(0.05),size=0.3,color='black') + xlab("avg logFC") + ylab("-log10(padj")+
    geom_text_repel(aes(label=ifelse(p_val_adj<0.05 & abs(avg_logFC)>2,gene,'')),size=1,color='#ca0020',fontface = 'bold', box.padding = unit(0.35, "lines"),point.padding = unit(0.5, "lines"),segment.color = 'grey50')

pseudobulk_fun <- function(seu_object,cellType_name,DE_method,DE_type,num_cols,figwidth,figheight){
    print(paste0('The total number of cells from ',cellType_name,' is: '))
    print(dim(seu_object)[2])
    seu_object$cell_type=rep(cellType_name,dim(seu_object)[2])

    matrices = to_pseudobulk(seu_object,min_cells = 0,min_reps = 0)
    #run DE analysis based on specific method 
    DE = run_de(seu_object,de_method=DE_method,de_type=DE_type)
    write.csv(DE,paste0('old/Libra/Libra_',DE_method,'_',DE_type,'_',cellType_name,'_allcells_DE.csv'))
    
    highlight_DE = DE[(abs(DE$avg_logFC)>2 & DE$p_val_adj<0.05),]
    highlight_DE = na.omit(highlight_DE)
    write.csv(highlight_DE,paste0('old/Libra/Libra_',DE_method,'_',DE_type,'_',cellType_name,'_allcells_sig_DGEs.csv'))
    print('The significant DEGs:')
    print(highlight_DE,n=2000)
    
    options(rePercentageFeatureSet.width = figwidth, repr.plot.height = figheight, repr.plot.res = 300,ggrepel.max.overlaps = Inf)
    #plot
    print('Control vs MS: negative avg logFC means down-regulated in Control and positive means up-regulated in Control.')
    
    ggplot(data = DE, aes(x=avg_logFC, y=-log10(p_val_adj)),colour=gene[abs(avg_logFC>2)])+ 
        geom_point(pch=20,size=0.8,color='#cccccc')+
        geom_point(data=highlight_DE,aes(x=avg_logFC,y=-log10(p_val_adj)),color='#0571b0',size=0.8,pch=20)+
        facet_wrap(~cell_type,ncol = num_cols)+
        theme_bw() + theme(legend.position="none")+
        geom_vline(xintercept = 2,size=0.3,color='black') + 
        geom_vline(xintercept = -2,size=0.3,color='black') + 
        geom_hline(yintercept = -log10(0.05),size=0.3,color='black') + 
        xlab("avg logFC") + ylab("-log10(padj")+
        #ggtitle(paste0(cellType_name,' excluding ',outlier))+
        #theme(plot.title = element_text(size=8,hjust = 0.5)) + 
        geom_text_repel(aes(label=ifelse(p_val_adj<0.05 & abs(avg_logFC)>2,gene,'')),
                        size=1.5,color='#ca0020',fontface = 'bold', box.padding = unit(0.35, "lines"),
                        point.padding = unit(0.5, "lines"),segment.color = 'grey50')
}

pseudobulk_fun_subcells <- function(seu_object,cellType_name,DE_method,DE_type,subcell_type,num_cols,figwidth,figheight){
    print(paste0('The total number of cells from ',cellType_name,' is: '))
    print(dim(seu_object)[2])
    
    print(table(seu_object@meta.data[subcell_type]))

    matrices = to_pseudobulk(seu_object,min_cells = 0,min_reps = 0, cell_type_col =subcell_type,replicate_col = 'replicate',label_col = 'label')
    #run DE analysis based on specific method 
    DE = run_de(seu_object,de_method=DE_method,de_type=DE_type,cell_type_col =subcell_type)
    write.csv(DE,paste0('old/Libra/Libra_',DE_method,'_',DE_type,'_',cellType_name,'_DE.csv'))
    
    highlight_DE = DE[(abs(DE$avg_logFC)>2 & DE$p_val_adj<0.05),]
    highlight_DE = na.omit(highlight_DE)
    write.csv(highlight_DE,paste0('old/Libra/Libra_',DE_method,'_',DE_type,'_',cellType_name,'_sig_DGEs.csv'))
    print('The significant DEGs:')
    print(highlight_DE,n=2000)
    
    options(rePercentageFeatureSet.width = figwidth, repr.plot.height = figheight, repr.plot.res = 300,ggrepel.max.overlaps = Inf)
    #plot
    print('Control vs MS: negative avg logFC means down-regulated in Control and positive means up-regulated in Control.')
    
    ggplot(data = DE, aes(x=avg_logFC, y=-log10(p_val_adj)),colour=gene[abs(avg_logFC>2)])+ 
        geom_point(pch=20,size=0.8,color='#cccccc')+
        geom_point(data=highlight_DE,aes(x=avg_logFC,y=-log10(p_val_adj)),color='#0571b0',size=0.8,pch=20)+
        facet_wrap(~cell_type,ncol = num_cols)+
        theme_bw() + theme(legend.position="none")+
        geom_vline(xintercept = 2,size=0.3,color='black') + 
        geom_vline(xintercept = -2,size=0.3,color='black') + 
        geom_hline(yintercept = -log10(0.05),size=0.3,color='black') + 
        xlab("avg logFC") + ylab("-log10(padj")+
        #ggtitle(paste0(cellType_name,' excluding ',outlier))+
        #theme(plot.title = element_text(size=8,hjust = 0.5)) + 
        geom_text_repel(aes(label=ifelse(p_val_adj<0.05 & abs(avg_logFC)>2,gene,'')),
                        size=1.5,color='#ca0020',fontface = 'bold', box.padding = unit(0.35, "lines"),
                        point.padding = unit(0.5, "lines"),segment.color = 'grey50')
}

colnames(alldata.filt@meta.data)

epith_gene_name = read.csv('old/Epithelial_gene_name.csv')

colnames(alldata.filt@meta.data)

epith <- subset(alldata.filt,features = epith_gene_name$Epithelia_gene,cells = which(alldata.filt$Epithelial_subcell_Ana!=alldata.filt$NeuroGlia_subcell_Ana[1]))

pseudobulk_fun(seu_object = epith,cellType_name  = 'Epithelial',DE_method = 'edgeR',DE_type = 'LRT',
                        num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun_subcells(seu_object = epith,cellType_name  = 'Epithelial',DE_method = 'edgeR',DE_type = 'LRT',
                        num_cols = 2,figwidth = 10,figheight = 6,subcell_type='Epithelial_subcell_Ana')

immune_gene_name = read.csv('old/Immune_gene_name.csv')

immune <- subset(alldata.filt,features = immune_gene_name$Immune_gene,cells = which(alldata.filt$Immune_subcell_Ana!=alldata.filt$NeuroGlia_subcell_Ana[1]))

pseudobulk_fun(seu_object = immune,cellType_name  = 'Immune',DE_method = 'edgeR',DE_type = 'LRT',num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun_subcells(seu_object = immune,cellType_name  = 'Immune',DE_method = 'edgeR',DE_type = 'LRT',
                        num_cols = 4,figwidth = 16,figheight = 10,subcell_type='Immune_subcell_Ana')

pseudobulk_fun_subcells_2 <- function(seu_object,cellType_name,DE_method,DE_type,num_cols,figwidth,figheight){
    print(paste0('The total number of cells from ',cellType_name,' is: '))
    print(dim(seu_object)[2])

    seu_object$cell_type=rep(cellType_name,dim(seu_object)[2])
    
    matrices = to_pseudobulk(seu_object,min_cells = 5,min_reps = 0)
    #run DE analysis based on specific method 
    DE = run_de(seu_object,de_method=DE_method,de_type=DE_type)
    write.csv(DE,paste0('old/Libra/Libra_',DE_method,'_',DE_type,'_',cellType_name,'_DE.csv'))
    
    highlight_DE = DE[(abs(DE$avg_logFC)>2 & DE$p_val_adj<0.05),]
    highlight_DE = na.omit(highlight_DE)
    write.csv(highlight_DE,paste0('old/Libra/Libra_',DE_method,'_',DE_type,'_',cellType_name,'_sig_DGEs.csv'))
    print('The significant DEGs:')
    print(highlight_DE,n=2000)
    
    options(rePercentageFeatureSet.width = figwidth, repr.plot.height = figheight, repr.plot.res = 300,ggrepel.max.overlaps = Inf)
    #plot
    print('Control vs MS: negative avg logFC means down-regulated in Control and positive means up-regulated in Control.')
    
    ggplot(data = DE, aes(x=avg_logFC, y=-log10(p_val_adj)),colour=gene[abs(avg_logFC>2)])+ 
        geom_point(pch=20,size=0.8,color='#cccccc')+
        geom_point(data=highlight_DE,aes(x=avg_logFC,y=-log10(p_val_adj)),color='#0571b0',size=0.8,pch=20)+
        facet_wrap(~cell_type,ncol = num_cols)+
        theme_bw() + theme(legend.position="none")+
        geom_vline(xintercept = 2,size=0.3,color='black') + 
        geom_vline(xintercept = -2,size=0.3,color='black') + 
        geom_hline(yintercept = -log10(0.05),size=0.3,color='black') + 
        xlab("avg logFC") + ylab("-log10(padj")+
        #ggtitle(paste0(cellType_name,' excluding ',outlier))+
        #theme(plot.title = element_text(size=8,hjust = 0.5)) + 
        geom_text_repel(aes(label=ifelse(p_val_adj<0.05 & abs(avg_logFC)>2,gene,'')),
                        size=1.5,color='#ca0020',fontface = 'bold', box.padding = unit(0.35, "lines"),
                        point.padding = unit(0.5, "lines"),segment.color = 'grey50')
}

Mac_object = subset(immune,cells=which((immune$Immune_subcell_Ana=='Activated Mac: chemotaxis high') | (immune$Immune_subcell_Ana=='Mac') |
                                       (immune$Immune_subcell_Ana=='Mac Proliferation') | (immune$Immune_subcell_Ana==' Migratory Mac') | 
                                       (immune$Immune_subcell_Ana=='Phagocyting Mac')))
pseudobulk_fun_subcells_2(seu_object = Mac_object,cellType_name  = 'Macrophage',DE_method = 'edgeR',DE_type = 'LRT',
                        num_cols = 1,figwidth = 2,figheight = 4)

Dendritic_object = subset(immune,cells=which((immune$Immune_subcell_Ana=='Dendritic') | (immune$Immune_subcell_Ana=='Dendritic cells ITGAX')))
pseudobulk_fun_subcells_2(seu_object = Dendritic_object,cellType_name  = 'Dendritic',DE_method = 'edgeR',DE_type = 'LRT',
                        num_cols = 1,figwidth = 2,figheight = 4)

immune_noDBL = subset(immune,cells = which((immune$Immune_subcell_Ana != 'Mac-Epithelial') & (immune$Immune_subcell_Ana != 'Mac-Mesenchymal')))
pseudobulk_fun_subcells_2(immune_noDBL,cellType_name  = 'Immune_noDBL',DE_method = 'edgeR',DE_type = 'LRT',
                        num_cols = 1,figwidth = 2,figheight = 4)

mesen_gene_name = read.csv('old/Mesenchymal_gene_name.csv')

mesen <- subset(alldata.filt,features = mesen_gene_name$Mesenchymal_gene,cells = which(alldata.filt$Mesenchymal_subcell_Ana!=alldata.filt$NeuroGlia_subcell_Ana[1]))

pseudobulk_fun(seu_object = mesen,cellType_name  = 'Mesenchymal',DE_method = 'edgeR',DE_type = 'LRT',num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun_subcells(seu_object = mesen,cellType_name  = 'Mesenchymal',DE_method = 'edgeR',DE_type = 'LRT',
                        num_cols = 2,figwidth = 10,figheight = 8,subcell_type='Mesenchymal_subcell_Ana')

FB_object = subset(mesen,cells=which((mesen$Mesenchymal_subcell_Ana=='FB: cell adhesion') | (mesen$Mesenchymal_subcell_Ana=='FB: extracelular matrix') | 
                                     (mesen$Mesenchymal_subcell_Ana=='FB: mobile') | (mesen$Mesenchymal_subcell_Ana=='FB: transmembrane transport')))
pseudobulk_fun_subcells_2(seu_object = FB_object,cellType_name  = 'Fibroblast',DE_method = 'edgeR',DE_type = 'LRT',
                        num_cols = 1,figwidth = 2,figheight = 4)


peri_object = subset(mesen,cells=which((mesen$Mesenchymal_subcell_Ana=='Pericytes') | (mesen$Mesenchymal_subcell_Ana=='vSMC')))
pseudobulk_fun_subcells_2(seu_object = peri_object,cellType_name  = 'Pericytes_vSMC',DE_method = 'edgeR',DE_type = 'LRT',
                        num_cols = 1,figwidth = 2,figheight = 4)

endoth_gene_name = read.csv('old/Endothelial_gene_name.csv')

endoth <- subset(alldata.filt,features = endoth_gene_name$Endothelial_gene,cells = which(alldata.filt$Endothelial_Subcells!=alldata.filt$NeuroGlia_subcell_Ana[1]))

pseudobulk_fun(seu_object = endoth,cellType_name  = 'Endothelial',DE_method = 'edgeR',DE_type = 'LRT',num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun_subcells(seu_object = endoth,cellType_name  = 'Endothelial',DE_method = 'edgeR',DE_type = 'LRT',
                        num_cols = 2,figwidth = 10,figheight = 8,subcell_type='Endothelial_Subcells')

neuroglia_gene_name = read.csv('old/NeuroGlia_gene_name.csv')

neuroglia <- subset(alldata.filt,features = neuroglia_gene_name$NeuroGlia_gene,cells = which(alldata.filt$cell_type=='NeuroGlia' & alldata.filt$NeuroGlia_subcell_Ana!=alldata.filt$NeuroGlia_subcell_Ana[1] ))

pseudobulk_fun(seu_object = neuroglia,cellType_name  = 'NeuroGlia',DE_method = 'edgeR',DE_type = 'LRT',num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun_subcells(seu_object = neuroglia,cellType_name  = 'NeuroGlia',DE_method = 'edgeR',DE_type = 'LRT',
                        num_cols = 3,figwidth = 12,figheight = 8,subcell_type='NeuroGlia_subcell_Ana')

astro_object = subset(neuroglia,cells=which((neuroglia$NeuroGlia_subcell_Ana=='Astrocytes') | (neuroglia$NeuroGlia_subcell_Ana=='Neuroblast: NRCAM, CD44+') | 
                                             (neuroglia$NeuroGlia_subcell_Ana=='NSC: GLI3+')))
pseudobulk_fun_subcells_2(seu_object = astro_object,cellType_name  = 'Astrocyes_neuroblast_nsc',DE_method = 'edgeR',DE_type = 'LRT',
                        num_cols = 1,figwidth = 2,figheight = 4)


ependyma_gene_name = read.csv('old/Ependyma_gene_name.csv')

ependyma <- subset(alldata.filt,features = ependyma_gene_name$Ependyma_gene,cells = which(alldata.filt$cell_type=='Ependyma'))

pseudobulk_fun(seu_object = ependyma,cellType_name  = 'Ependyma',DE_method = 'edgeR',DE_type = 'LRT',num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun_subcells_2(seu_object = ependyma,cellType_name  = 'Ependyma',DE_method = 'edgeR',DE_type = 'LRT',
                        num_cols = 1,figwidth = 2,figheight = 4)

matrices = to_pseudobulk(alldata.filt,min_cells = 5,min_reps = 0)

DE = run_de(alldata.filt,de_method = 'edgeR',de_type = 'QLF')

write.csv(DE,'old/Libra/Libra_edgeR_QLF_allcells_DE.csv')

highlight_DE = DE[(abs(DE$avg_logFC)>2 & DE$p_val_adj<0.05),]
write.csv(highlight_DE,paste0('old/Libra/Libra_edgeR_QLF_allcells_sig_DGEs.csv'))

print(highlight_DE,n=2000)

options(rePercentageFeatureSet.width = 14, repr.plot.height = 6, repr.plot.res = 300,ggrepel.max.overlaps = Inf)
print('Control vs MS: negative avg logFC means down-regulated in Control and positive means up-regulated in Control.')

ggplot(data = DE, aes(x=avg_logFC, y=-log10(p_val_adj)),colour=gene[abs(avg_logFC>2)])+ geom_point(pch=20,size=0.8,color='#cccccc')+geom_point(data=highlight_DE,aes(x=avg_logFC,y=-log10(p_val_adj)),color='#0571b0',size=0.8,pch=20)+
    facet_wrap(~cell_type,ncol = 3)+theme_bw() + theme(legend.position="none")+
    geom_vline(xintercept = 2,size=0.3,color='black') + geom_vline(xintercept = -2,size=0.3,color='black') + geom_hline(yintercept = -log10(0.05),size=0.3,color='black') + xlab("avg logFC") + ylab("-log10(padj")+
    geom_text_repel(aes(label=ifelse(p_val_adj<0.05 & abs(avg_logFC)>2,gene,'')),size=1,color='#ca0020',fontface = 'bold', box.padding = unit(0.35, "lines"),point.padding = unit(0.5, "lines"),segment.color = 'grey50')

pseudobulk_fun(seu_object = epith,cellType_name  = 'Epithelial',DE_method = 'edgeR',DE_type = 'QLF',num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun_subcells(seu_object = epith,cellType_name  = 'Epithelial',DE_method = 'edgeR',DE_type = 'QLF',
                        num_cols = 2,figwidth = 10,figheight = 6,subcell_type='Epithelial_subcell_Ana')

pseudobulk_fun(seu_object = immune,cellType_name  = 'Immune',DE_method = 'edgeR',DE_type = 'QLF',num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun_subcells(seu_object = immune,cellType_name  = 'Immune',DE_method = 'edgeR',DE_type = 'QLF',
                        num_cols = 4,figwidth = 16,figheight = 10,subcell_type='Immune_subcell_Ana')

Mac_object = subset(immune,cells=which((immune$Immune_subcell_Ana=='Activated Mac: chemotaxis high') | (immune$Immune_subcell_Ana=='Mac') |
                                       (immune$Immune_subcell_Ana=='Mac Proliferation') | (immune$Immune_subcell_Ana==' Migratory Mac') | 
                                       (immune$Immune_subcell_Ana=='Phagocyting Mac')))
pseudobulk_fun_subcells_2(seu_object = Mac_object,cellType_name  = 'Macrophage',DE_method = 'edgeR',DE_type = 'QLF',
                        num_cols = 1,figwidth = 2,figheight = 4)

Dendritic_object = subset(immune,cells=which((immune$Immune_subcell_Ana=='Dendritic') | (immune$Immune_subcell_Ana=='Dendritic cells ITGAX')))
pseudobulk_fun_subcells_2(seu_object = Dendritic_object,cellType_name  = 'Dendritic',DE_method = 'edgeR',DE_type = 'QLF',
                        num_cols = 1,figwidth = 2,figheight = 4)

immune_noDBL = subset(immune,cells = which((immune$Immune_subcell_Ana != 'Mac-Epithelial') & (immune$Immune_subcell_Ana != 'Mac-Mesenchymal')))
pseudobulk_fun_subcells_2(immune_noDBL,cellType_name  = 'Immune_noDBL',DE_method = 'edgeR',DE_type = 'QLF',
                        num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun(seu_object = mesen,cellType_name  = 'Mesenchymal',DE_method = 'edgeR',DE_type = 'QLF',num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun_subcells(seu_object = mesen,cellType_name  = 'Mesenchymal',DE_method = 'edgeR',DE_type = 'QLF',
                        num_cols = 2,figwidth = 10,figheight = 8,subcell_type='Mesenchymal_subcell_Ana')

FB_object = subset(mesen,cells=which((mesen$Mesenchymal_subcell_Ana=='FB: cell adhesion') | (mesen$Mesenchymal_subcell_Ana=='FB: extracelular matrix') | 
                                     (mesen$Mesenchymal_subcell_Ana=='FB: mobile') | (mesen$Mesenchymal_subcell_Ana=='FB: transmembrane transport')))
pseudobulk_fun_subcells_2(seu_object = FB_object,cellType_name  = 'Fibroblast',DE_method = 'edgeR',DE_type = 'QLF',
                        num_cols = 1,figwidth = 2,figheight = 4)


peri_object = subset(mesen,cells=which((mesen$Mesenchymal_subcell_Ana=='Pericytes') | (mesen$Mesenchymal_subcell_Ana=='vSMC')))
pseudobulk_fun_subcells_2(seu_object = peri_object,cellType_name  = 'Pericytes_vSMC',DE_method = 'edgeR',DE_type = 'QLF',
                        num_cols = 1,figwidth = 2,figheight = 4)

endoth_gene_name = read.csv('old/Endothelial_gene_name.csv')

endoth <- subset(alldata.filt,features = endoth_gene_name$Endothelial_gene,cells = which(alldata.filt$Endothelial_Subcells!=alldata.filt$NeuroGlia_subcell_Ana[1]))

pseudobulk_fun(seu_object = endoth,cellType_name  = 'Endothelial',DE_method = 'edgeR',DE_type = 'QLF',num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun_subcells(seu_object = endoth,cellType_name  = 'Endothelial',DE_method = 'edgeR',DE_type = 'QLF',
                        num_cols = 2,figwidth = 10,figheight = 8,subcell_type='Endothelial_Subcells')

neuroglia_gene_name = read.csv('old/NeuroGlia_gene_name.csv')

neuroglia <- subset(alldata.filt,features = neuroglia_gene_name$NeuroGlia_gene,cells = which(alldata.filt$cell_type=='NeuroGlia' & alldata.filt$NeuroGlia_subcell_Ana!=alldata.filt$NeuroGlia_subcell_Ana[1] ))

pseudobulk_fun(seu_object = neuroglia,cellType_name  = 'NeuroGlia',DE_method = 'edgeR',DE_type = 'QLF',num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun_subcells(seu_object = neuroglia,cellType_name  = 'NeuroGlia',DE_method = 'edgeR',DE_type = 'QLF',
                        num_cols = 3,figwidth = 12,figheight = 8,subcell_type='NeuroGlia_subcell_Ana')

astro_object = subset(neuroglia,cells=which((neuroglia$NeuroGlia_subcell_Ana=='Astrocytes') | (neuroglia$NeuroGlia_subcell_Ana=='Neuroblast: NRCAM, CD44+') | 
                                             (neuroglia$NeuroGlia_subcell_Ana=='NSC: GLI3+')))
pseudobulk_fun_subcells_2(seu_object = astro_object,cellType_name  = 'Astrocyes_neuroblast_nsc',DE_method = 'edgeR',DE_type = 'QLF',
                        num_cols = 1,figwidth = 2,figheight = 4)


ependyma_gene_name = read.csv('old/Ependyma_gene_name.csv')

ependyma <- subset(alldata.filt,features = ependyma_gene_name$Ependyma_gene,cells = which(alldata.filt$cell_type=='Ependyma'))

pseudobulk_fun(seu_object = ependyma,cellType_name  = 'Ependyma',DE_method = 'edgeR',DE_type = 'QLF',num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun_subcells_2(seu_object = ependyma,cellType_name  = 'Ependyma',DE_method = 'edgeR',DE_type = 'QLF',
                        num_cols = 1,figwidth = 2,figheight = 4)

matrices = to_pseudobulk(alldata.filt,min_cells = 5,min_reps = 0)

DE = run_de(alldata.filt,de_method = 'DESeq2',de_type = 'LRT')

highlight_DE = DE[abs(DE$avg_logFC)>2 & DE$p_val_adj<0.05,]
highlight_DE = na.omit(highlight_DE)
print(highlight_DE,n=2000)

write.csv(DE,paste0('old/Libra/Libra_DEseq2_LRT_allcells_DE.csv'))
write.csv(highlight_DE,'old/Libra/Libra_DEseq2_LRT_allcells_sig_DGEs.csv')

options(rePercentageFeatureSet.width = 14, repr.plot.height = 6, repr.plot.res = 300,ggrepel.max.overlaps = Inf)
print('Control vs MS: negative avg logFC means down-regulated in Control and positive means up-regulated in Control.')

ggplot(data = DE, aes(x=avg_logFC, y=-log10(p_val_adj)),colour=gene[abs(avg_logFC>2)])+ geom_point(pch=20,size=0.8,color='#cccccc')+geom_point(data=highlight_DE,aes(x=avg_logFC,y=-log10(p_val_adj)),color='#0571b0',size=0.8,pch=20)+
    facet_wrap(~cell_type,ncol = 3)+theme_bw() + theme(legend.position="none")+
    geom_vline(xintercept = 2,size=0.3,color='black') + geom_vline(xintercept = -2,size=0.3,color='black') + geom_hline(yintercept = -log10(0.05),size=0.3,color='black') + xlab("avg logFC") + ylab("-log10(padj")+
    geom_text_repel(aes(label=ifelse(p_val_adj<0.05 & abs(avg_logFC)>2,gene,'')),size=1,color='#ca0020',fontface = 'bold', box.padding = unit(0.35, "lines"),point.padding = unit(0.5, "lines"),segment.color = 'grey50')



pseudobulk_fun(seu_object = epith,cellType_name  = 'Epithelial',DE_method = 'DESeq2',DE_type = 'LRT',num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun_subcells(seu_object = epith,cellType_name  = 'Epithelial',DE_method = 'DESeq2',DE_type = 'LRT',
                        num_cols = 2,figwidth = 10,figheight = 6,subcell_type='Epithelial_subcell_Ana')

pseudobulk_fun(seu_object = immune,cellType_name  = 'Immune',DE_method = 'DESeq2',DE_type = 'LRT',num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun_subcells(seu_object = immune,cellType_name  = 'Immune',DE_method = 'DESeq2',DE_type = 'LRT',
                        num_cols = 4,figwidth = 16,figheight = 10,subcell_type='Immune_subcell_Ana')

pseudobulk_fun_subcells_2 <- function(seu_object,cellType_name,DE_method,DE_type,num_cols,figwidth,figheight){
    print(paste0('The total number of cells from ',cellType_name,' is: '))
    print(dim(seu_object)[2])

    seu_object$cell_type=rep(cellType_name,dim(seu_object)[2])
    
    matrices = to_pseudobulk(seu_object,min_cells = 5,min_reps = 0)
    #run DE analysis based on specific method 
    DE = run_de(seu_object,de_method=DE_method,de_type=DE_type)
    write.csv(DE,paste0('old/Libra/Libra_',DE_method,'_',DE_type,'_',cellType_name,'_DE.csv'))
    
    highlight_DE = DE[(abs(DE$avg_logFC)>2 & DE$p_val_adj<0.05),]
    highlight_DE = na.omit(highlight_DE)
    write.csv(highlight_DE,paste0('old/Libra/Libra_',DE_method,'_',DE_type,'_',cellType_name,'_sig_DGEs.csv'))
    print('The significant DEGs:')
    print(highlight_DE,n=2000)
    
    options(rePercentageFeatureSet.width = figwidth, repr.plot.height = figheight, repr.plot.res = 300,ggrepel.max.overlaps = Inf)
    #plot
    print('Control vs MS: negative avg logFC means down-regulated in Control and positive means up-regulated in Control.')
    
    ggplot(data = DE, aes(x=avg_logFC, y=-log10(p_val_adj)),colour=gene[abs(avg_logFC>2)])+ 
        geom_point(pch=20,size=0.8,color='#cccccc')+
        geom_point(data=highlight_DE,aes(x=avg_logFC,y=-log10(p_val_adj)),color='#0571b0',size=0.8,pch=20)+
        facet_wrap(~cell_type,ncol = num_cols)+
        theme_bw() + theme(legend.position="none")+
        geom_vline(xintercept = 2,size=0.3,color='black') + 
        geom_vline(xintercept = -2,size=0.3,color='black') + 
        geom_hline(yintercept = -log10(0.05),size=0.3,color='black') + 
        xlab("avg logFC") + ylab("-log10(padj")+
        #ggtitle(paste0(cellType_name,' excluding ',outlier))+
        #theme(plot.title = element_text(size=8,hjust = 0.5)) + 
        geom_text_repel(aes(label=ifelse(p_val_adj<0.05 & abs(avg_logFC)>2,gene,'')),
                        size=1.5,color='#ca0020',fontface = 'bold', box.padding = unit(0.35, "lines"),
                        point.padding = unit(0.5, "lines"),segment.color = 'grey50')
}

Mac_object = subset(immune,cells=which((immune$Immune_subcell_Ana=='Activated Mac: chemotaxis high') | (immune$Immune_subcell_Ana=='Mac') |
                                       (immune$Immune_subcell_Ana=='Mac Proliferation') | (immune$Immune_subcell_Ana==' Migratory Mac') | 
                                       (immune$Immune_subcell_Ana=='Phagocyting Mac')))
pseudobulk_fun_subcells_2(seu_object = Mac_object,cellType_name  = 'Macrophage',DE_method = 'DESeq2',DE_type = 'LRT',
                        num_cols = 1,figwidth = 2,figheight = 4)

Dendritic_object = subset(immune,cells=which((immune$Immune_subcell_Ana=='Dendritic') | (immune$Immune_subcell_Ana=='Dendritic cells ITGAX')))
pseudobulk_fun_subcells_2(seu_object = Dendritic_object,cellType_name  = 'Dendritic',DE_method = 'DESeq2',DE_type = 'LRT',
                        num_cols = 1,figwidth = 2,figheight = 4)

immune_noDBL = subset(immune,cells = which((immune$Immune_subcell_Ana != 'Mac-Epithelial') & (immune$Immune_subcell_Ana != 'Mac-Mesenchymal')))
pseudobulk_fun_subcells_2(immune_noDBL,cellType_name  = 'Immune_noDBL',DE_method = 'DESeq2',DE_type = 'LRT',
                        num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun(seu_object = mesen,cellType_name  = 'Mesenchymal',DE_method = 'DESeq2',DE_type = 'LRT',num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun_subcells(seu_object = mesen,cellType_name  = 'Mesenchymal',DE_method = 'DESeq2',DE_type = 'LRT',
                        num_cols = 2,figwidth = 10,figheight = 8,subcell_type='Mesenchymal_subcell_Ana')

FB_object = subset(mesen,cells=which((mesen$Mesenchymal_subcell_Ana=='FB: cell adhesion') | (mesen$Mesenchymal_subcell_Ana=='FB: extracelular matrix') | 
                                     (mesen$Mesenchymal_subcell_Ana=='FB: mobile') | (mesen$Mesenchymal_subcell_Ana=='FB: transmembrane transport')))
pseudobulk_fun_subcells_2(seu_object = FB_object,cellType_name  = 'Fibroblast',DE_method = 'DESeq2',DE_type = 'LRT',
                        num_cols = 1,figwidth = 2,figheight = 4)


peri_object = subset(mesen,cells=which((mesen$Mesenchymal_subcell_Ana=='Pericytes') | (mesen$Mesenchymal_subcell_Ana=='vSMC')))
pseudobulk_fun_subcells_2(seu_object = peri_object,cellType_name  = 'Pericytes_vSMC',DE_method = 'DESeq2',DE_type = 'LRT',
                        num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun(seu_object = endoth,cellType_name  = 'Endothelial',DE_method = 'DESeq2',DE_type = 'LRT',num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun_subcells(seu_object = endoth,cellType_name  = 'Endothelial',DE_method = 'DESeq2',DE_type = 'LRT',
                        num_cols = 2,figwidth = 10,figheight = 8,subcell_type='Endothelial_Subcells')

pseudobulk_fun(seu_object = neuroglia,cellType_name  = 'NeuroGlia',DE_method = 'DESeq2',DE_type = 'LRT',num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun_subcells(seu_object = neuroglia,cellType_name  = 'NeuroGlia',DE_method = 'DESeq2',DE_type = 'LRT',
                        num_cols = 3,figwidth = 12,figheight = 8,subcell_type='NeuroGlia_subcell_Ana')

astro_object = subset(neuroglia,cells=which((neuroglia$NeuroGlia_subcell_Ana=='Astrocytes') | (neuroglia$NeuroGlia_subcell_Ana=='Neuroblast: NRCAM, CD44+') | 
                                             (neuroglia$NeuroGlia_subcell_Ana=='NSC: GLI3+')))
pseudobulk_fun_subcells_2(seu_object = astro_object,cellType_name  = 'Astrocyes_neuroblast_nsc',DE_method = 'DESeq2',DE_type = 'LRT',
                        num_cols = 1,figwidth = 2,figheight = 4)


pseudobulk_fun(seu_object = ependyma,cellType_name  = 'Ependyma',DE_method = 'DESeq2',DE_type = 'LRT',num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun_subcells_2(seu_object = ependyma,cellType_name  = 'Ependyma',DE_method = 'DESeq2',DE_type = 'LRT',
                        num_cols = 1,figwidth = 2,figheight = 4)

matrices = to_pseudobulk(alldata.filt,min_cells = 5,min_reps = 0)

DE = run_de(alldata.filt,de_method = 'DESeq2',de_type = 'Wald')

highlight_DE = DE[abs(DE$avg_logFC)>2 & DE$p_val_adj<0.05,]
highlight_DE = na.omit(highlight_DE)
print(highlight_DE,n=2000)

write.csv(DE,paste0('old/Libra/Libra_DEseq2_Wald_allcells_DE.csv'))
write.csv(highlight_DE,'old/Libra/Libra_DEseq2_Wald_allcells_sig_DGEs.csv')

options(rePercentageFeatureSet.width = 14, repr.plot.height = 6, repr.plot.res = 300,ggrepel.max.overlaps = Inf)
print('Control vs MS: negative avg logFC means down-regulated in Control and positive means up-regulated in Control.')

ggplot(data = DE, aes(x=avg_logFC, y=-log10(p_val_adj)),colour=gene[abs(avg_logFC>2)])+ geom_point(pch=20,size=0.8,color='#cccccc')+geom_point(data=highlight_DE,aes(x=avg_logFC,y=-log10(p_val_adj)),color='#0571b0',size=0.8,pch=20)+
    facet_wrap(~cell_type,ncol = 3)+theme_bw() + theme(legend.position="none")+
    geom_vline(xintercept = 2,size=0.3,color='black') + geom_vline(xintercept = -2,size=0.3,color='black') + geom_hline(yintercept = -log10(0.05),size=0.3,color='black') + xlab("avg logFC") + ylab("-log10(padj")+
    geom_text_repel(aes(label=ifelse(p_val_adj<0.05 & abs(avg_logFC)>2,gene,'')),size=1,color='#ca0020',fontface = 'bold', box.padding = unit(0.35, "lines"),point.padding = unit(0.5, "lines"),segment.color = 'grey50')


epith <- subset(alldata.filt,features = epith_gene_name$Epithelia_gene,cells = which(alldata.filt$Epithelial_subcell_Ana!=alldata.filt$NeuroGlia_subcell_Ana[1]))

pseudobulk_fun(seu_object = epith,cellType_name  = 'Epithelial',DE_method = 'DESeq2',DE_type = 'Wald',num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun_subcells(seu_object = epith,cellType_name  = 'Epithelial',DE_method = 'DESeq2',DE_type = 'Wald',
                        num_cols = 2,figwidth = 10,figheight = 6,subcell_type='Epithelial_subcell_Ana')

immune <- subset(alldata.filt,features = immune_gene_name$Immune_gene,cells = which(alldata.filt$Immune_subcell_Ana!=alldata.filt$NeuroGlia_subcell_Ana[1]))

pseudobulk_fun(seu_object = immune,cellType_name  = 'Immune',DE_method = 'DESeq2',DE_type = 'Wald',num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun_subcells(seu_object = immune,cellType_name  = 'Immune',DE_method = 'DESeq2',DE_type = 'Wald',
                        num_cols = 4,figwidth = 16,figheight = 10,subcell_type='Immune_subcell_Ana')

Mac_object = subset(immune,cells=which((immune$Immune_subcell_Ana=='Activated Mac: chemotaxis high') | (immune$Immune_subcell_Ana=='Mac') |
                                       (immune$Immune_subcell_Ana=='Mac Proliferation') | (immune$Immune_subcell_Ana==' Migratory Mac') | 
                                       (immune$Immune_subcell_Ana=='Phagocyting Mac')))

pseudobulk_fun_subcells_2(seu_object = Mac_object,cellType_name  = 'Macrophage',DE_method = 'DESeq2',DE_type = 'Wald',
                        num_cols = 1,figwidth = 2,figheight = 4)

Dendritic_object = subset(immune,cells=which((immune$Immune_subcell_Ana=='Dendritic') | (immune$Immune_subcell_Ana=='Dendritic cells ITGAX')))

pseudobulk_fun_subcells_2(seu_object = Dendritic_object,cellType_name  = 'Dendritic',DE_method = 'DESeq2',DE_type = 'Wald',
                        num_cols = 1,figwidth = 2,figheight = 4)

immune_noDBL = subset(immune,cells = which((immune$Immune_subcell_Ana != 'Mac-Epithelial') & (immune$Immune_subcell_Ana != 'Mac-Mesenchymal')))
pseudobulk_fun_subcells_2(immune_noDBL,cellType_name  = 'Immune_noDBL',DE_method = 'DESeq2',DE_type = 'Wald',
                        num_cols = 1,figwidth = 2,figheight = 4)

mesen <- subset(alldata.filt,features = mesen_gene_name$Mesenchymal_gene,cells = which(alldata.filt$Mesenchymal_subcell_Ana!=alldata.filt$NeuroGlia_subcell_Ana[1]))

pseudobulk_fun(seu_object = mesen,cellType_name  = 'Mesenchymal',DE_method = 'DESeq2',DE_type = 'Wald',num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun_subcells(seu_object = mesen,cellType_name  = 'Mesenchymal',DE_method = 'DESeq2',DE_type = 'Wald',
                        num_cols = 2,figwidth = 10,figheight = 8,subcell_type='Mesenchymal_subcell_Ana')

FB_object = subset(mesen,cells=which((mesen$Mesenchymal_subcell_Ana=='FB: cell adhesion') | (mesen$Mesenchymal_subcell_Ana=='FB: extracelular matrix') | 
                                     (mesen$Mesenchymal_subcell_Ana=='FB: mobile') | (mesen$Mesenchymal_subcell_Ana=='FB: transmembrane transport')))
pseudobulk_fun_subcells_2(seu_object = FB_object,cellType_name  = 'Fibroblast',DE_method = 'DESeq2',DE_type = 'Wald',
                        num_cols = 1,figwidth = 2,figheight = 4)


peri_object = subset(mesen,cells=which((mesen$Mesenchymal_subcell_Ana=='Pericytes') | (mesen$Mesenchymal_subcell_Ana=='vSMC')))
pseudobulk_fun_subcells_2(seu_object = peri_object,cellType_name  = 'Pericytes_vSMC',DE_method = 'DESeq2',DE_type = 'Wald',
                        num_cols = 1,figwidth = 2,figheight = 4)

endoth <- subset(alldata.filt,features = endoth_gene_name$Endothelial_gene,cells = which(alldata.filt$Endothelial_Subcells!=alldata.filt$NeuroGlia_subcell_Ana[1]))

pseudobulk_fun(seu_object = endoth,cellType_name  = 'Endothelial',DE_method = 'DESeq2',DE_type = 'Wald',num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun_subcells(seu_object = endoth,cellType_name  = 'Endothelial',DE_method = 'DESeq2',DE_type = 'Wald',
                        num_cols = 2,figwidth = 10,figheight = 8,subcell_type='Endothelial_Subcells')

neuroglia <- subset(alldata.filt,features = neuroglia_gene_name$NeuroGlia_gene,cells = which(alldata.filt$cell_type=='NeuroGlia' & alldata.filt$NeuroGlia_subcell_Ana!=alldata.filt$NeuroGlia_subcell_Ana[1] ))

pseudobulk_fun(seu_object = neuroglia,cellType_name  = 'NeuroGlia',DE_method = 'DESeq2',DE_type = 'Wald',num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun_subcells(seu_object = neuroglia,cellType_name  = 'NeuroGlia',DE_method = 'DESeq2',DE_type = 'Wald',
                        num_cols = 3,figwidth = 12,figheight = 8,subcell_type='NeuroGlia_subcell_Ana')

astro_object = subset(neuroglia,cells=which((neuroglia$NeuroGlia_subcell_Ana=='Astrocytes') | (neuroglia$NeuroGlia_subcell_Ana=='Neuroblast: NRCAM, CD44+') | 
                                             (neuroglia$NeuroGlia_subcell_Ana=='NSC: GLI3+')))
pseudobulk_fun_subcells_2(seu_object = astro_object,cellType_name  = 'Astrocyes_neuroblast_nsc',DE_method = 'DESeq2',DE_type = 'Wald',
                        num_cols = 1,figwidth = 2,figheight = 4)


ependyma <- subset(alldata.filt,features = ependyma_gene_name$Ependyma_gene,cells = which(alldata.filt$cell_type=='Ependyma'))

pseudobulk_fun(seu_object = ependyma,cellType_name  = 'Ependyma',DE_method = 'DESeq2',DE_type = 'Wald',num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun_subcells_2(seu_object = ependyma,cellType_name  = 'Ependyma',DE_method = 'DESeq2',DE_type = 'Wald',
                        num_cols = 1,figwidth = 2,figheight = 4)

matrices = to_pseudobulk(alldata.filt,min_cells = 5,min_reps = 0)

# transpose matricies
matrices_transpose=lapply(matrices,function(x)t(x))

names(matrices_transpose)

pca_plot_fun <- function(pca_mat,celltype){
    pca <- prcomp(pca_mat, scale = TRUE,center = TRUE, retx = T)
    pca_df = as.data.frame(pca$x)
    pca_df$label = rep(c('Control','MS'),c(length(grep('Control',rownames(pca_df))),length(grep('MS',rownames(pca_df)))))
    pca_df$batch = rownames(pca_df)
    pca_df$batch = gsub(':Control','',pca_df$batch)
    pca_df$batch = gsub(':MS','',pca_df$batch)
    pca_df$batch = factor(pca_df$batch,levels = c('R100_63_Ctrl_F_Age56','R100_64_Ctrl_M_Age64','R100_72_Ctrl_F_Age82','R108_30_Ctrl_M_Age88','R108_31_Ctrl_F_Age82',
                                          'R77_61_Ctrl_F_Age79','Spain_18_006_Ctrl_M_Age51','Spain_19_005_Ctrl_M_Age60','R100_61_MS_F_Age75','R100_62_MS_F_Age62',
                                         'R100_71_MS_F_Age52','R77_60_MS_M_Age46'))
    print(pca_df[,1:3])
    options(rePercentageFeatureSet.width = 4, repr.plot.height = 5, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
    ggplot(pca_df,aes(x=PC1,y=PC2,colour=batch,shape=label))+geom_point(size=3)+ggtitle(celltype)

}

pca_plot_fun(matrices_transpose[[1]],'Endothelial')

pca_plot_fun(matrices_transpose[[2]],'Ependyma')

pca_plot_fun(matrices_transpose[[3]],'Epithelial')

pca_plot_fun(matrices_transpose[[4]],'Immune')

pca_plot_fun(matrices_transpose[[5]],'Mesenchymal')

pca_plot_fun(matrices_transpose[[6]],'NeuroGlia')

pseudobulk_fun_noOutlier <- function(seu_object,cellType_name,DE_method,DE_type,num_cols,figwidth,figheight){
    print(paste0('The total number of cells from ',cellType_name,' is: '))
    print(dim(seu_object)[2])
    seu_object$cell_type=rep(cellType_name,dim(seu_object)[2])

    matrices = to_pseudobulk(seu_object,min_cells = 5,min_reps = 0)
    #run DE analysis based on specific method 
    DE = run_de(seu_object,de_method=DE_method,de_type=DE_type)
    write.csv(DE,paste0('old/Libra/Libra_',DE_method,'_',DE_type,'_',cellType_name,'_noOutlier_allcells_DE.csv'))
    
    highlight_DE = DE[(abs(DE$avg_logFC)>2 & DE$p_val_adj<0.05),]
    highlight_DE = na.omit(highlight_DE)
    write.csv(highlight_DE,paste0('old/Libra/Libra_',DE_method,'_',DE_type,'_',cellType_name,'_noOutlier_allcells_sig_DGEs.csv'))
    print('The significant DEGs:')
    print(highlight_DE,n=2000)
    
    options(rePercentageFeatureSet.width = figwidth, repr.plot.height = figheight, repr.plot.res = 300,ggrepel.max.overlaps = Inf)
    #plot
    print('Control vs MS: negative avg logFC means down-regulated in Control and positive means up-regulated in Control.')
    
    ggplot(data = DE, aes(x=avg_logFC, y=-log10(p_val_adj)),colour=gene[abs(avg_logFC>2)])+ 
        geom_point(pch=20,size=0.8,color='#cccccc')+
        geom_point(data=highlight_DE,aes(x=avg_logFC,y=-log10(p_val_adj)),color='#0571b0',size=0.8,pch=20)+
        facet_wrap(~cell_type,ncol = num_cols)+
        theme_bw() + theme(legend.position="none")+
        geom_vline(xintercept = 2,size=0.3,color='black') + 
        geom_vline(xintercept = -2,size=0.3,color='black') + 
        geom_hline(yintercept = -log10(0.05),size=0.3,color='black') + 
        xlab("avg logFC") + ylab("-log10(padj")+
        #ggtitle(paste0(cellType_name,' excluding ',outlier))+
        #theme(plot.title = element_text(size=8,hjust = 0.5)) + 
        geom_text_repel(aes(label=ifelse(p_val_adj<0.05 & abs(avg_logFC)>2,gene,'')),
                        size=1.5,color='#ca0020',fontface = 'bold', box.padding = unit(0.35, "lines"),
                        point.padding = unit(0.5, "lines"),segment.color = 'grey50')
}

pseudobulk_fun_subcells_noOutlier <- function(seu_object,cellType_name,DE_method,DE_type,subcell_type,num_cols,figwidth,figheight){
    print(paste0('The total number of cells from ',cellType_name,' is: '))
    print(dim(seu_object)[2])
    
    print(table(seu_object@meta.data[subcell_type]))

    matrices = to_pseudobulk(seu_object,min_cells = 5,min_reps = 0, cell_type_col =subcell_type,replicate_col = 'replicate',label_col = 'label')
    #run DE analysis based on specific method 
    DE = run_de(seu_object,de_method=DE_method,de_type=DE_type,cell_type_col =subcell_type)
    write.csv(DE,paste0('old/Libra/Libra_',DE_method,'_',DE_type,'_',cellType_name,'_noOutlier_DE.csv'))
    
    highlight_DE = DE[(abs(DE$avg_logFC)>2 & DE$p_val_adj<0.05),]
    highlight_DE = na.omit(highlight_DE)
    write.csv(highlight_DE,paste0('old/Libra/Libra_',DE_method,'_',DE_type,'_',cellType_name,'_noOutlier_sig_DGEs.csv'))
    print('The significant DEGs:')
    print(highlight_DE,n=2000)
    
    options(rePercentageFeatureSet.width = figwidth, repr.plot.height = figheight, repr.plot.res = 300,ggrepel.max.overlaps = Inf)
    #plot
    print('Control vs MS: negative avg logFC means down-regulated in Control and positive means up-regulated in Control.')
    
    ggplot(data = DE, aes(x=avg_logFC, y=-log10(p_val_adj)),colour=gene[abs(avg_logFC>2)])+ 
        geom_point(pch=20,size=0.8,color='#cccccc')+
        geom_point(data=highlight_DE,aes(x=avg_logFC,y=-log10(p_val_adj)),color='#0571b0',size=0.8,pch=20)+
        facet_wrap(~cell_type,ncol = num_cols)+
        theme_bw() + theme(legend.position="none")+
        geom_vline(xintercept = 2,size=0.3,color='black') + 
        geom_vline(xintercept = -2,size=0.3,color='black') + 
        geom_hline(yintercept = -log10(0.05),size=0.3,color='black') + 
        xlab("avg logFC") + ylab("-log10(padj")+
        #ggtitle(paste0(cellType_name,' excluding ',outlier))+
        #theme(plot.title = element_text(size=8,hjust = 0.5)) + 
        geom_text_repel(aes(label=ifelse(p_val_adj<0.05 & abs(avg_logFC)>2,gene,'')),
                        size=1.5,color='#ca0020',fontface = 'bold', box.padding = unit(0.35, "lines"),
                        point.padding = unit(0.5, "lines"),segment.color = 'grey50')
}

colnames(alldata.filt@meta.data)

endoth_edit = subset(endoth,cells = which(endoth$replicate!='R108_30_Ctrl_M_Age88' & endoth$replicate!='R100_71_MS_F_Age52'))

pseudobulk_fun_noOutlier(seu_object = endoth,cellType_name  = 'Endothelial',DE_method = 'edgeR',DE_type = 'LRT',
                        num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun_subcells_noOutlier(seu_object = endoth,cellType_name  = 'Endothelial',DE_method = 'edgeR',DE_type = 'LRT',
                        num_cols = 2,figwidth = 10,figheight = 6,subcell_type='Endothelial_Subcells')

ependyma_edit = subset(ependyma,cells = which(ependyma$cell_type=='Ependyma' & ependyma$replicate!='R108_30_Ctrl_M_Age88'))

pseudobulk_fun_noOutlier(seu_object = ependyma_edit,cellType_name  = 'Ependyma',DE_method = 'edgeR',DE_type = 'LRT',
                        num_cols = 1,figwidth = 2,figheight = 4)

epith_edit = subset(epith,cells = which(epith$replicate!='Spain_18_006_Ctrl_M_Age51' &  epith$replicate!='Spain_19_005_Ctrl_M_Age60' &  epith$replicate!='R77_60_MS_M_Age46'))

pseudobulk_fun_noOutlier(seu_object = epith,cellType_name  = 'Epithelial',DE_method = 'edgeR',DE_type = 'LRT',
                        num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun_subcells_noOutlier(seu_object = epith,cellType_name  = 'Epithelial',DE_method = 'edgeR',DE_type = 'LRT',
                        num_cols = 2,figwidth = 10,figheight = 6,subcell_type='Epithelial_subcell_Ana')

immune_edit = subset(immune,cells = which(immune$replicate!='R77_61_Ctrl_F_Age79' &  immune$replicate!='R100_71_MS_F_Age52'))

pseudobulk_fun_noOutlier(seu_object = immune_edit,cellType_name  = 'Immune',DE_method = 'edgeR',DE_type = 'LRT',
                        num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun_subcells_noOutlier(seu_object = immune_edit,cellType_name  = 'Immune',DE_method = 'edgeR',DE_type = 'LRT',
                        num_cols = 2,figwidth = 16,figheight = 16,subcell_type='Immune_subcell_Ana')

mesen_edit = subset(mesen,cells = which(mesen$replicate!='Spain_19_005_Ctrl_M_Age60' &  mesen$replicate!='R77_61_Ctrl_F_Age79'))

pseudobulk_fun_noOutlier(seu_object = mesen_edit,cellType_name  = 'Mesenchymal',DE_method = 'edgeR',DE_type = 'LRT',
                        num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun_subcells_noOutlier(seu_object = mesen_edit,cellType_name  = 'Mesenchymal',DE_method = 'edgeR',DE_type = 'LRT',
                        num_cols = 2,figwidth = 16,figheight = 10,subcell_type='Mesenchymal_subcell_Ana')

neuroglia_edit = subset(neuroglia,cells = which(alldata.filt$replicate!='R100_64_Ctrl_M_Age64' & alldata.filt$replicate!='R77_60_MS_M_Age46'))

pseudobulk_fun_noOutlier(seu_object = neuroglia,cellType_name  = 'NeuroGlia',DE_method = 'edgeR',DE_type = 'LRT',
                        num_cols = 1,figwidth = 2,figheight = 4)

pseudobulk_fun_subcells_noOutlier(seu_object = neuroglia,cellType_name  = 'NeuroGlia',DE_method = 'edgeR',DE_type = 'LRT',
                        num_cols = 2,figwidth = 16,figheight = 10,subcell_type='NeuroGlia_subcell_Ana')




