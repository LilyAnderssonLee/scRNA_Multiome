suppressPackageStartupMessages({
    library(Seurat)
    library(CountClust)
    library(Matrix)
    library(dplyr)
    library(SeuratDisk)
    library(loomR)
    library(dplyr)
    library(rafalib)
    library(data.table)
    library(ggplot2)
    library(enrichR)
    library(repr)
    library(scater)
    library(MAST)
    })


#set figures size
library(repr)
options(repr.plot.width=10, repr.plot.height=10)

library(Seurat)

epith<-readRDS('../Seurat/results/Epithelial_seurat.rds')
epith<-SetIdent(epith,value = 'Condition')

seu_file= epith
fdata=as.data.frame(rownames(seu_file),row.names = rownames(seu_file))
colnames(fdata)<-'primerid'
sca<-FromMatrix(log2(as.matrix(seu_file@assays$RNA@counts)+1),seu_file@meta.data,fdata)

#load data 
load('../Seurat/results/Epithelial_MAST.Rdata')
summaryDt <- summaryCond$datatable

# Make a datatable with all the relevant output - combined Hurdle model

fcHurdle <- merge(summaryDt[contrast=='ConditionPD' & component=='H',.(primerid, `Pr(>Chisq)`)], ##hurdle P values
     summaryDt[contrast=='ConditionPD' & component=='logFC', 
               .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

# change name of Pr(>Chisq) to fdr and sort by fdr
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
#install.packages('data.table')
library(data.table)
setorder(fcHurdle, fdr)

write.table(fcHurdle,'../DE_comparisons/markers/Epithelial/Epithelial_MAST_PD_vs_Control.table',sep="\t",quote=F)

fcHurdle<-read.table('../DE_comparisons/markers/Epithelial/Epithelial_MAST_PD_vs_Control.table')
fcHurdleSig<-fcHurdle[fcHurdle$fdr<0.05,]
#fcHurdleSig <- merge(fcHurdle[(fdr<0.05 & abs(coef)>0.25)],as.data.table(mcols(sca)), by='primerid')# 
#setorder(fcHurdleSig, fdr)

dim(fcHurdleSig)

# run all DE methods
methods <- c("wilcox","bimod","t","LR","MAST","roc")
DE <- list()
for (m in methods){
    outfile <- paste0("../DE_comparisons/markers/Epithelial/Epithelial_seurat_",m,"_PD_vs_Control.table")
    DE[[m]]<- FindMarkers(object = epith,`ident.1` = 'PD',`ident.2` = 'Control',
                          test.use = m,verbose=FALSE)
    write.table(DE[[m]],file=outfile,sep="\t",quote=F)
}

DE$MAST_pkg<-fcHurdleSig
rownames(DE$MAST_pkg)<-fcHurdleSig$primerid

methods <- c("negbinom","poisson")
for (m in methods){
    outfile <- paste0("../DE_comparisons/markers/Epithelial/Epithelial_seurat_",m,"_PD_vs_Control.table")
    DE[[m]]<- FindMarkers(object = epith,`ident.1` = 'PD',`ident.2` = 'Control',
                          test.use = m,verbose=FALSE,slot='counts')
    write.table(DE[[m]],file=outfile,sep="\t",quote=F)
}

DE <- list()
files <- c("../DE_comparisons/markers/Epithelial/Epithelial_seurat_wilcox_PD_vs_Control.table",
           "../DE_comparisons/markers/Epithelial/Epithelial_seurat_bimod_PD_vs_Control.table",
           "../DE_comparisons/markers/Epithelial/Epithelial_seurat_t_PD_vs_Control.table",
           "../DE_comparisons/markers/Epithelial/Epithelial_seurat_LR_PD_vs_Control.table",
           "../DE_comparisons/markers/Epithelial/Epithelial_seurat_MAST_PD_vs_Control.table",
           "../DE_comparisons/markers/Epithelial/Epithelial_seurat_roc_PD_vs_Control.table",
           "../DE_comparisons/markers/Epithelial/Epithelial_seurat_negbinom_PD_vs_Control.table",
           "../DE_comparisons/markers/Epithelial/Epithelial_seurat_poisson_PD_vs_Control.table",
           "../DE_comparisons/markers/Epithelial/Epithelial_MAST_PD_vs_Control.table")

for (i in 1:9){ 
  DE[[i]]<-read.table(files[i],header=T)
}
DE[[10]]<-read.csv("../DE_comparisons/markers/Epithelial/Epithelial_Scapy_wilcoxon.csv")

names(DE)<-c("seurat-wilcox", "seurat-bimod","seurat-t","seurat-LR","seurat-MAST","seurat_roc",
             "seurat_negbinom","seurat_poisson","MAST","scanpy-wilcox")

rownames(DE$MAST)<-DE$MAST$primerid

DE$`scanpy-wilcox`<-DE$`scanpy-wilcox`[DE$`scanpy-wilcox`$group=='PD',]
rownames(DE$`scanpy-wilcox`)<-DE$`scanpy-wilcox`$names

overlap_phyper<-function(L,bg=length(unique(unlist(L))),with.unique=TRUE,plot=FALSE){
   # takes a list with indices and creates a matrix with overlap between elements
   # can also plot as a heatmap with coloring according to significance in overlap.
   # phyper test uses all entries in L as background if bg is not specified.

   nL<-length(L)
   M<-mat.or.vec(nL,nL)
   P<-mat.or.vec(nL,nL)
   P[,]<-1
   nU<-mat.or.vec(nL,1)
   for (i in 1:nL){
       nU[i]<-length(setdiff(L[[i]],unlist(L[-i])))
       for (j in  i:nL){
           M[i,j]<-length(intersect(L[[i]],L[[j]]))
           P[i,j]<-1-phyper(M[i,j],length(L[[i]]),bg-length(L[[i]]),length(L[[j]]))
           if (i==j) {P[i,j]<-NA}
      }
   }
   if (with.unique){
      M<-cbind(M,nU)
      colnames(M)<-c(names(L),"unique")
      P<-cbind(P,rep(1,length(nU)))
      colnames(P)<-c(names(L),"unique")
   }else {
      colnames(M)<-names(L)
   }

   rownames(M)<-names(L)
   rownames(P)<-names(L)
   if (plot){
      library(gplots)
      lab<-matrix(as.character(M),nrow=nL)
      lab[is.na(lab)]<-''
      par(oma=c(4,1,2,4),xpd=T,cex=0.1,mfrow=c(1,1))
      cexC = 0.2 + 0.2/log10(nL)
      notecex = 0.2 + 0.2/log10(nL)
      h<-heatmap.2(-log10(P+1e-16),cellnote=lab,scale="none",trace="none",density.info="none",
                   notecex=1,notecol="black",dendrogram="none",Colv=F,Rowv=F,key=T,key.xlab="-log10(p.value)",key.title='')
   }
   return(list(overlap=M,pval=P))
}

lapply(DE,function(x) dim(x))

# get top 200 genes for each test
top.200 <- lapply(DE,function(x) rownames(x)[1:200])

# load a function for plotting overlap
options(repr.plot.width=10, repr.plot.height=10)

o <- overlap_phyper(top.200,plot=T,bg=nrow(DE$MAST))

DE$MAST<-DE$MAST[abs(DE$MAST$coef)>0.25,]
DE$`scanpy-wilcox`<-DE$`scanpy-wilcox`[abs(DE$`scanpy-wilcox`$logfoldchanges)>0.25,]

# not really a p-value for the ROC test, so skip for now 
# (5th entry in the list)
pval.column <- c(5,5,5,5,5,5,5,6,5) # define which column contains the p-value
names(pval.column)<-names(DE)[-6]
sign.genes <- list()
cutP<-0.05
for (n in names(pval.column)){
  sg <- which(DE[[n]][,pval.column[n]] < cutP)
  sign.genes[[n]]<-rownames(DE[[n]])[sg]
}


lapply(sign.genes,length)

options(repr.plot.width=10, repr.plot.height=10)
o <- overlap_phyper(sign.genes,plot=T,bg=nrow(DE$`scanpy-wilcox`))

sign.genes_overlap=Reduce(intersect, sign.genes)
plot_df=DE$`seurat-wilcox`[sign.genes_overlap,]

mypar(mar = c(4, 6, 3, 1))
barplot(sort(setNames(plot_df$avg_log2FC,rownames(plot_df)),F),
        horiz = T,las=1,main = "Up/down-regulated genes in PD Immune",border = "white",yaxs="i",xlab=list('avg_log2FC',cex=1, font=2))

DotPlot(epith, features = names(sort(setNames(plot_df$avg_log2FC,rownames(plot_df)),F)), group.by = 'Condition', assay = "RNA") + coord_flip()

# Perform enrichment
dbs=c('GO_Biological_Process_2017b','GO_Biological_Process_2018','KEGG_2019_Human','WikiPathways_2019_Human')
enrich_results <- enrichr(genes = sign.genes_overlap, databases = dbs)[[1]]

par(mfrow = c(1, 1), mar = c(5, 30, 2, 1))
barplot(height = -log10(enrich_results$Adjusted.P.value)[40:1], names.arg = enrich_results$Term[40:1], 
    horiz = TRUE, las = 1, border = FALSE, cex.names = 0.6,xlab=list('-log10(Adj.P.value)',cex=1,font=2))
abline(v = c(-log10(0.05)), lty = 2)
text(x=2.5,y=49,label = 'Adj.P.value = 0.05',font=2,col='red')

library(msigdbr)

# Download gene sets
msigdbr_species()
all_gene_sets <- msigdbr(species = "Homo sapiens")
all_gene_sets <- as.data.frame(all_gene_sets)
# List available gene sets
unique(all_gene_sets$gs_subcat)


# Subset which gene set you want to use.
all_gene_subset <- all_gene_sets[which(all_gene_sets$gs_subcat == "CP:WIKIPATHWAYS" |all_gene_sets$gs_subcat == "GO:BP" | all_gene_sets$gs_subcat == "'CP:KEGG'" ), ]
gmt <- lapply(unique(all_gene_subset$gs_name), function(x) {
  all_gene_subset[all_gene_subset$gs_name == x, "gene_symbol"]
})
names(gmt) <- unique(paste0(all_gene_subset$gs_name, "_", all_gene_subset$gs_exact_source))


# Create a gene rank based on the gene expression fold change
sign.genes_overlap_df=DE$`seurat-wilcox`[sign.genes_overlap,]
#only select abs(avg_log2FC)>1
#sign.genes_overlap_df<-sign.genes_overlap_df[abs(sign.genes_overlap_df$avg_log2FC)>1,]
#dim(sign.genes_overlap_df)

epith_gene_rank <- setNames(sign.genes_overlap_df$avg_log2FC, casefold(rownames(sign.genes_overlap_df), upper = T))

# Perform enrichemnt analysis
fgseaRes_sig_markers_epith <- fgseaMultilevel(pathways = gmt, stats = epith_gene_rank, minSize = 1, maxSize = Inf,
                                               eps=0,nPermSimple = 5000)

fgseaRes_sig_markers_epith <- fgseaRes_sig_markers_epith[order(fgseaRes_sig_markers_epith$NES, decreasing = T), ]


fgseaRes_sig_markers_epith<-fgseaRes_sig_markers_epith[fgseaRes_sig_markers_epith$padj<0.05,]

names(DE)
genes_overlap=Reduce(intersect,list(rownames(DE[[1]]),rownames(DE[[2]]),rownames(DE[[3]]),rownames(DE[[4]]),rownames(DE[[5]])))
length(genes_overlap)

# Perform enrichment
dbs=c('GO_Biological_Process_2017b','GO_Biological_Process_2018','KEGG_2019_Human','WikiPathways_2019_Human')
enrich_results_seurat_wilcox <- enrichr(genes = rownames(DE$`seurat-wilcox`), databases = dbs)[[1]]

enrich_results_seurat_wilcox_sig<-enrich_results_seurat_wilcox[enrich_results_seurat_wilcox$Adjusted.P.value<0.05,]

par(mfrow = c(1, 1), mar = c(5, 30, 2, 1))
barplot(height = -log10(enrich_results_seurat_wilcox_sig$Adjusted.P.value)[40:1], names.arg = enrich_results_seurat_wilcox_sig$Term[40:1], 
    horiz = TRUE, las = 1, border = FALSE, cex.names = 0.6,xlab=list('-log10(Adj.P.value)',cex=1,font=2))

barplot(height = -log10(enrich_results_seurat_wilcox_sig$Adjusted.P.value)[80:41], names.arg = enrich_results_seurat_wilcox_sig$Term[80:41], 
    horiz = TRUE, las = 1, border = FALSE, cex.names = 0.6,xlab=list('-log10(Adj.P.value)',cex=1,font=2))

# Create a gene rank based on the gene expression fold change
gene_rank_seurat_wilcox <- setNames(DE$`seurat-wilcox`$avg_log2FC, casefold(rownames(DE$`seurat-wilcox`), upper = T))

library(msigdbr)

# Download gene sets
msigdbgmt <- msigdbr::msigdbr("Homo sapiens")
msigdbgmt <- as.data.frame(msigdbgmt)

# List available gene sets
unique(msigdbgmt$gs_subcat)

msigdbgmt_subset <- msigdbgmt[which(msigdbgmt$gs_subcat == "CP:WIKIPATHWAYS"|msigdbgmt$gs_subcat =='CP:KEGG'|msigdbgmt$gs_subcat =='GO:BP'), ]
gmt <- lapply(unique(msigdbgmt_subset$gs_name), function(x) {
    msigdbgmt_subset[msigdbgmt_subset$gs_name == x, "gene_symbol"]
})
names(gmt) <- unique(paste0(msigdbgmt_subset$gs_name, "_", msigdbgmt_subset$gs_exact_source))

library(fgsea)

#https://www.biostars.org/p/479821/
# Perform enrichemnt analysis
fgseaRes_seurat_wilcox <- fgseaMultilevel(pathways = gmt, stats = gene_rank_seurat_wilcox, minSize = 1, maxSize = Inf,
                           eps=0,nPermSimple = 1000)


fgseaRes_seurat_wilcox <- fgseaRes_seurat_wilcox[order(fgseaRes_seurat_wilcox$NES, decreasing = T), ]




fgseaRes_seurat_wilcox_sig<-fgseaRes_seurat_wilcox[fgseaRes_seurat_wilcox$padj<0.05,]

par(mfrow = c(1, 1), mar = c(5, 30, 2, 1))
barplot(height = fgseaRes_seurat_wilcox_sig$NES, names.arg = fgseaRes_seurat_wilcox_sig$pathway, 
    horiz = TRUE, las = 1, border = FALSE, cex.names = 0.5,xlab=list('NES',cex=1,font=2))

# select significant genes for each methods: Seurat_wilcox was filtered by logfc=0.25, MAST was filted by fdr<0.05 and logfc=0.25
# scanpy_wilcox was filtered by pvals_adj and logfoldchanges=0.25
enrich_input=list()
enrich_input$seurat_wilcox<-DE$`seurat-wilcox`
enrich_input$MAST<-DE$MAST[which(abs(DE$MAST$coef)>0.25 & DE$MAST$fdr<0.05),]
enrich_input$scanpy_wilcox<-DE$`scanpy-wilcox`[which(DE$`scanpy-wilcox`$group=='PD' & abs(DE$`scanpy-wilcox`$logfoldchanges)>0.25),]

enrich_out=list()
# Perform enrichment
dbs=c('GO_Biological_Process_2017b','GO_Biological_Process_2018','KEGG_2019_Human','WikiPathways_2019_Human')

enrich_out<-lapply(enrich_input,function(x)enrichr(genes = rownames(x), databases = dbs)[[1]])

# select significant pathways
enrich_out_sig<-lapply(enrich_out,function(x)x<-x[x$Adjusted.P.value<0.05,])
#select pathways were found in all three data set
enrich_out_sig_overlap_pathway<-Reduce(intersect,list(enrich_out_sig$MAST$Term,enrich_out_sig$scanpy_wilcox$Term,enrich_out_sig$seurat_wilcox$Term))
enrich_out_sig_overlap_df<-enrich_out_sig$seurat_wilcox[which(enrich_out_sig$seurat_wilcox$Term %in% enrich_out_sig_overlap_pathway),]   


#plot 
par(mfrow = c(1, 1), mar = c(5, 30, 2, 1))
barplot(height = -log10(enrich_out_sig_overlap_df$Adjusted.P.value)[dim(enrich_out_sig_overlap_df)[1]:1], 
        names.arg = enrich_out_sig_overlap_df$Term[dim(enrich_out_sig_overlap_df)[1]:1], 
        horiz = TRUE, las = 1, border = FALSE, cex.names = 0.6,xlab=list('-log10(Adj.P.value)',cex=1,font=2),
       main='Enrichr: pathways enriched in PD Epith')


library(ReactomePA)

x <- enrichPathway(gene=enrich_input$seurat_wilcox,pvalueCutoff=0.05, readable=T)
head(summary(x))

require(biomaRt)

gene_list<-list()

gene_list$seurat_wilcox<- setNames(enrich_input$seurat_wilcox$avg_log2FC, casefold(rownames(enrich_input$seurat_wilcox), upper = T))
gene_list$MAST<-setNames(enrich_input$MAST$coef,casefold(rownames(enrich_input$MAST), upper = T))
gene_list$scanpy_wilcox<-setNames(enrich_input$scanpy_wilcox$logfoldchanges,casefold(rownames(enrich_input$scanpy_wilcox), upper = T))


mart <- useMart("ensembl")
#head(listDatasets(mart))
mart <- useDataset("hsapiens_gene_ensembl", mart = mart)
#convert gene symbol into Gene id
genes.table <- lapply(gene_list,function(x){
    getBM(filters = 'external_gene_name',
          attributes= c("ensembl_gene_id", "external_gene_name"), 
          values= names(x), mart= mart)})

genes.table<-lapply(genes.table,function(x)x[!duplicated(x$external_gene_name),])

## prepare the genelist for the reactomPA
for (filename in c('seurat_wilcox','MAST','scanpy_wilcox')){
    genes.table[[filename]]<-genes.table[[filename]][genes.table[[filename]]$external_gene_name %in% names(gene_list[[filename]]),]
    logfc<-list()
    for (i in 1:dim(genes.table[[filename]])[1]){
        logfc[[i]]<-gene_list[[filename]][names(gene_list[[filename]])==genes.table[[filename]]$external_gene_name[i]]
    }
    genes.table[[filename]]$logfc<-unlist(logfc)
}

for (f in c('seurat_wilcox','MAST','scanpy_wilcox')){
    genes.table[[f]]$gene_id<-gsub("ENSG00000|ENSG000000|ENSG0000000|ENSG00000000|ENSG000000000|ENSG0000000000",
                                   '',genes.table[[f]]$ensembl_gene_id)
}

gene_list_edited<-list()
for (f in c('seurat_wilcox','MAST','scanpy_wilcox')){
    gene_list_edited[[f]]<-genes.table[[f]]$logfc
    names(gene_list_edited[[f]])<-genes.table[[f]]$gene_id
    gene_list_edited[[f]]<-sort(gene_list_edited[[f]],decreasing = T)
}

de<-lapply(gene_list_edited,names)

library(ReactomePA)

reactom_enrich<-lapply(de,function(x)enrichPathway(gene=x,pvalueCutoff=0.05, readable=T))

lapply(reactom_enrich,function(x)dim(summary(x)))

library(DOSE)
library(igraph)
library(clusterProfiler)

summary(reactom_enrich$scanpy_wilcox)

cnetplot(reactom_enrich$scanpy_wilcox, categorySize="pvalue", 
         foldChange=gene_list_edited$scanpy_wilcox)

viewPathway("SUMOylation of DNA methylation proteins", 
            readable = TRUE, 
            foldChange = gene_list_edited$scanpy_wilcox)

viewPathway("SUMOylation of RNA binding proteins", 
            readable = TRUE, 
            foldChange = gene_list_edited$scanpy_wilcox)

# Create a gene rank based on the gene expression fold change

gene_rank<-list()

gene_rank$seurat_wilcox<- setNames(enrich_input$seurat_wilcox$avg_log2FC, casefold(rownames(enrich_input$seurat_wilcox), upper = T))
gene_rank$MAST<-setNames(enrich_input$MAST$coef,casefold(rownames(enrich_input$MAST), upper = T))
gene_rank$scanpy_wilcox<-setNames(enrich_input$scanpy_wilcox$logfoldchanges,casefold(rownames(enrich_input$scanpy_wilcox), upper = T))

fgseaRes_epith <- lapply(gene_rank,function(x){
    fgseaMultilevel(pathways = gmt, stats = x, minSize = 1, maxSize = Inf,
                           eps=0,nPermSimple = 5000)})



fgseaRes_epith_sorted <- lapply(fgseaRes_epith,function(x){x<-x[order(x$NES, decreasing = T),]})
fgseaRes_epith_overlap_pathway<-Reduce(intersect,list(fgseaRes_epith_sorted$MAST$pathway,
                                                      fgseaRes_epith_sorted$scanpy_wilcox$pathway,
                                                      fgseaRes_epith_sorted$seurat_wilcox$pathway))
# over lap pathways
fgseaRes_epith_overlap_pathway_df<-fgseaRes_epith_sorted$MAST[which(fgseaRes_epith_sorted$MAST$pathway %in% fgseaRes_epith_overlap_pathway),]
#significant pathways
fgseaRes_epith_overlap_pathway_sig<-fgseaRes_epith_overlap_pathway_df[fgseaRes_epith_overlap_pathway_df$padj<0.05,]

par(mfrow = c(1, 1), mar = c(5, 30, 2, 1))
barplot(height = fgseaRes_epith_overlap_pathway_sig$NES, names.arg = fgseaRes_epith_overlap_pathway_sig$pathway, 
    horiz = TRUE, las = 1, border = FALSE, cex.names = 0.6,xlab=list('NES',cex=1,font=2),
       main='RankGenes:up/down-regulated in PD Epith')

save(plot_df,enrich_results,enrich_out,enrich_out_sig_overlap_df,fgseaRes_epith_overlap_pathway,fgseaRes_epith_overlap_pathway_sig,
     file='epith_DE.RData')

immune<-readRDS('../Seurat/results/Immune_seurat.rds')
immune<-SetIdent(immune,value = 'Condition')

seu_file= immune
fdata=as.data.frame(rownames(seu_file),row.names = rownames(seu_file))
colnames(fdata)<-'primerid'
sca<-FromMatrix(log2(as.matrix(seu_file@assays$RNA@counts)+1),seu_file@meta.data,fdata)

#load data 
load('../Seurat/results/Immune_MAST.Rdata') # 
summaryDt <- summaryCond$datatable

# Make a datatable with all the relevant output - combined Hurdle model

fcHurdle <- merge(summaryDt[contrast=='ConditionPD' & component=='H',.(primerid, `Pr(>Chisq)`)], ##hurdle P values
     summaryDt[contrast=='ConditionPD' & component=='logFC', 
               .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

# change name of Pr(>Chisq) to fdr and sort by fdr
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
#install.packages('data.table')
library(data.table)
setorder(fcHurdle, fdr)
write.table(fcHurdle,'../DE_comparisons/markers/Immune/Immune_MAST_PD_vs_Control.table',sep="\t",quote=F)
rownames(fcHurdle)<-fcHurdle$primerid


# run all DE methods
methods <- c("wilcox","bimod","t","LR","MAST","roc")
DE <- list()
for (m in methods){
    outfile <- paste0("../DE_comparisons/markers/Immune/Immune_seurat_",m,"_PD_vs_Control.table")
    DE[[m]]<- FindMarkers(object = immune,`ident.1` = 'PD',`ident.2` = 'Control',
                          test.use = m,verbose=FALSE)
    write.table(DE[[m]],file=outfile,sep="\t",quote=F)
}

methods <- c("negbinom","poisson")
for (m in methods){
    outfile <- paste0("../DE_comparisons/markers/Immune/Immune_seurat_",m,"_PD_vs_Control.table")
    DE[[m]]<- FindMarkers(object = immune,`ident.1` = 'PD',`ident.2` = 'Control',
                          test.use = m,verbose=FALSE,slot='counts')
    write.table(DE[[m]],file=outfile,sep="\t",quote=F)
}



DE <- list()
files <- c("../DE_comparisons/markers/Immune//Immune_seurat_wilcox_PD_vs_Control.table",
           "../DE_comparisons/markers/Immune/Immune_seurat_bimod_PD_vs_Control.table",
           "../DE_comparisons/markers/Immune//Immune_seurat_t_PD_vs_Control.table",
           "../DE_comparisons/markers/Immune//Immune_seurat_LR_PD_vs_Control.table",
           "../DE_comparisons/markers/Immune//Immune_seurat_MAST_PD_vs_Control.table",
           "../DE_comparisons/markers/Immune//Immune_seurat_roc_PD_vs_Control.table",
           "../DE_comparisons/markers/Immune//Immune_seurat_negbinom_PD_vs_Control.table",
           "../DE_comparisons/markers//Immune//Immune_seurat_poisson_PD_vs_Control.table",
           "../DE_comparisons/markers/Immune/Immune_MAST_PD_vs_Control.table")

for (i in 1:9){ 
  DE[[i]]<-read.table(files[i],header=T)
}

DE[[10]]<-read.csv( "../DE_comparisons/markers/Immune/Immune_Scapy_wilcoxon.csv",header=T)
names(DE)<-c("seurat_wilcox", "seurat_bimod","seurat_t","seurat_LR","seurat_MAST","seurat_roc",
             "seurat_negbinom","seurat_poisson","MAST","scanpy_wilcox")

DE$scanpy_wilcox<-DE$scanpy_wilcox[DE$scanpy_wilcox$group=='PD',]
rownames(DE$scanpy_wilcox)<-DE$scanpy_wilcox$names
rownames(DE$MAST)<-DE$MAST$primerid

# get top 200 genes for each test
top.200 <- lapply(DE,function(x) rownames(x)[1:200])
o <- overlap_phyper(top.200,plot=T,bg=nrow(DE$MAST))

names(DE)

pval.column <- c(5,5,5,5,5,5,5) # define which column contains the p-value
names(pval.column)<-names(DE)[-c(6,9,10)]
sign.genes <- list()
cutP<-0.05
for (n in names(pval.column)){
  sg <- which(DE[[n]][,pval.column[n]] < cutP)
  sign.genes[[n]]<-rownames(DE[[n]])[sg]
}

sign.genes$MAST<-DE$MAST$primerid[which(abs(DE$MAST$coef)>0.25 & DE$MAST$fdr<0.05)]
sign.genes$scanpy_wilcox<-DE$scanpy_wilcox$names[which(abs(DE$scanpy_wilcox$logfoldchanges)>0.25 & DE$scanpy_wilcox$pvals_adj<0.05)]

o <- overlap_phyper(sign.genes,plot=T,bg=nrow(DE$scanpy_wilcox))

sign.genes_overlap=Reduce(intersect, sign.genes)
sign.genes_overlap_df=DE$seurat_wilcox[sign.genes_overlap,]
#only select abs(avg_log2FC)>1
sign.genes_overlap_df<-sign.genes_overlap_df[abs(sign.genes_overlap_df$avg_log2FC)>1,]
dim(sign.genes_overlap_df)

#up-regulated in PD
sign.genes_overlap_df_up<-sign.genes_overlap_df[sign.genes_overlap_df$avg_log2FC>1,]
mypar(mar = c(4, 6, 3, 1))
barplot(sort(setNames(sign.genes_overlap_df_up$avg_log2FC,rownames(sign.genes_overlap_df_up)),F),
        horiz = T,las=1,border = "white",yaxs="i",xlab=list('avg_log2FC',cex=1, font=2),
       main='Up-regulated genes in PD Immune')

sign.genes_overlap_df_down<-sign.genes_overlap_df[sign.genes_overlap_df$avg_log2FC<1,]
mypar(mar = c(4, 6, 3, 1))
barplot(sort(setNames(sign.genes_overlap_df_down$avg_log2FC,
                      rownames(sign.genes_overlap_df_down)),T),
        horiz = T,las=1,border = "white",yaxs="i",xlab=list('avg_log2FC',cex=1, font=2),
       main='Down-regulated genes in Immune')


DotPlot(immune, features = names(sort(setNames(sign.genes_overlap_df_up$avg_log2FC,rownames(sign.genes_overlap_df_up)),F)), 
        group.by = 'Condition', assay = "RNA") + coord_flip()
DotPlot(immune, features = names(sort(setNames(sign.genes_overlap_df_down$avg_log2FC,rownames(sign.genes_overlap_df_down)),F)), 
        group.by = 'Condition', assay = "RNA") + coord_flip()


# Perform enrichment
dbs=c('GO_Biological_Process_2017b','GO_Biological_Process_2018','KEGG_2019_Human','WikiPathways_2019_Human')
enrich_results <- enrichr(genes = sign.genes_overlap, databases = dbs)[[1]]

#significant pathways
enrich_results_sig<-enrich_results[enrich_results$Adjusted.P.value<0.05,]
par(mfrow = c(1, 1), mar = c(5, 30, 2, 1))
barplot(height = -log10(enrich_results_sig$Adjusted.P.value)[40:1], names.arg = enrich_results_sig$Term[40:1], 
    horiz = TRUE, las = 1, border = FALSE, cex.names = 0.7,xlab=list('-log10(Adj.P.value)',cex=1,font=2),
       main='Enrichr for sig_markers in PD Immune')

par(mfrow = c(1, 1), mar = c(5, 30, 2, 1))
barplot(height = -log10(enrich_results_sig$Adjusted.P.value)[80:41], names.arg = enrich_results_sig$Term[80:41], 
    horiz = TRUE, las = 1, border = FALSE, cex.names = 0.5,xlab=list('-log10(Adj.P.value)',cex=1,font=2),
       main='Enrichr for sig_markers in PD Immune')


# Create a gene rank based on the gene expression fold change
immune_gene_rank <- setNames(sign.genes_overlap_df$avg_log2FC, casefold(rownames(sign.genes_overlap_df), upper = T))

# Perform enrichemnt analysis
fgseaRes_sig_markers_immune <- fgseaMultilevel(pathways = gmt, stats = immune_gene_rank, minSize = 1, maxSize = Inf,
                                               eps=0,nPermSimple = 5000)

fgseaRes_sig_markers_immune <- fgseaRes_sig_markers_immune[order(fgseaRes_sig_markers_immune$NES, decreasing = T), ]


fgseaRes_sig_markers_immune<-fgseaRes_sig_markers_immune[fgseaRes_sig_markers_immune$padj<0.05,]

par(mfrow = c(1, 1), mar = c(5, 30, 2, 1))
barplot(height = fgseaRes_sig_markers_immune$NES[7:1], names.arg = fgseaRes_sig_markers_immune$pathway[7:1], 
    horiz = TRUE, las = 1, border = FALSE, cex.names =0.8,xlab=list('NES',cex=1,font=2),
       main='Up/down-regulated in PD Immune')

# select significant genes for each methods: Seurat_wilcox was filtered by logfc=0.25, MAST was filted by fdr<0.05 and logfc=0.25
# scanpy_wilcox was filtered by pvals_adj and logfoldchanges=0.25
enrich_input=list()
enrich_input$seurat_wilcox<-DE$seurat_wilcox
enrich_input$MAST<-DE$MAST[which(abs(DE$MAST$coef)>0.25 & DE$MAST$fdr<0.05),]
enrich_input$scanpy_wilcox<-DE$scanpy_wilcox[abs(DE$scanpy_wilcox$logfoldchanges)>0.25,]
rownames(enrich_input$MAST)<-enrich_input$MAST$primerid

enrich_out_immune=list()
# Perform enrichment
dbs=c('GO_Biological_Process_2017b','GO_Biological_Process_2018','KEGG_2019_Human','WikiPathways_2019_Human')

enrich_out_immune<-lapply(enrich_input,function(x)enrichr(genes = rownames(x), databases = dbs)[[1]])

# select significant pathways
enrich_out_immune_sig<-lapply(enrich_out_immune,function(x)x<-x[x$Adjusted.P.value<0.05,])
#select pathways were found in all three data set
enrich_out_immune_sig_overlap_pathway<-Reduce(intersect,list(enrich_out_immune_sig$MAST$Term,enrich_out_immune_sig$scanpy_wilcox$Term,enrich_out_immune_sig$seurat_wilcox$Term))
enrich_out_immune_sig_overlap_df<-enrich_out_immune_sig$seurat_wilcox[which(enrich_out_immune_sig$seurat_wilcox$Term %in% enrich_out_immune_sig_overlap_pathway),]   



#plot 
par(mfrow = c(1, 1), mar = c(5, 30, 2, 1))
barplot(height = -log10(enrich_out_immune_sig_overlap_df$Adjusted.P.value)[40:1], 
        names.arg = enrich_out_immune_sig_overlap_df$Term[40:1], 
        horiz = TRUE, las = 1, border = FALSE, cex.names = 0.6,xlab=list('-log10(Adj.P.value)',cex=1,font=2),
       main='Enrichr:pathways enriched in PD Immune')

barplot(height = -log10(enrich_out_immune_sig_overlap_df$Adjusted.P.value)[80:41], 
        names.arg = enrich_out_immune_sig_overlap_df$Term[80:41], 
        horiz = TRUE, las = 1, border = FALSE, cex.names = 0.5,xlab=list('-log10(Adj.P.value)',cex=1,font=2),
       main='Enrichr:pathways enriched in PD Immune')


require(biomaRt)
library(ReactomePA)
library(DOSE)
library(igraph)
library(clusterProfiler)

gene_list<-list()

gene_list$seurat_wilcox<- setNames(enrich_input$seurat_wilcox$avg_log2FC, casefold(rownames(enrich_input$seurat_wilcox), upper = T))
gene_list$MAST<-setNames(enrich_input$MAST$coef,casefold(rownames(enrich_input$MAST), upper = T))
gene_list$scanpy_wilcox<-setNames(enrich_input$scanpy_wilcox$logfoldchanges,casefold(rownames(enrich_input$scanpy_wilcox), upper = T))

mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart = mart)
#convert gene symbol into Gene id
genes.table <- lapply(gene_list,function(x){
    getBM(filters = 'external_gene_name',
          attributes= c("ensembl_gene_id", "external_gene_name"), 
          values= names(x), mart= mart)})

genes.table<-lapply(genes.table,function(x)x[!duplicated(x$external_gene_name),])

## prepare the genelist for the reactomPA
for (filename in c('seurat_wilcox','MAST','scanpy_wilcox')){
    genes.table[[filename]]<-genes.table[[filename]][genes.table[[filename]]$external_gene_name %in% names(gene_list[[filename]]),]
    logfc<-list()
    for (i in 1:dim(genes.table[[filename]])[1]){
        logfc[[i]]<-gene_list[[filename]][names(gene_list[[filename]])==genes.table[[filename]]$external_gene_name[i]]
    }
    genes.table[[filename]]$logfc<-unlist(logfc)
}

for (f in c('seurat_wilcox','MAST','scanpy_wilcox')){
    genes.table[[f]]$gene_id<-gsub("ENSG00000|ENSG000000|ENSG0000000|ENSG00000000|ENSG000000000|ENSG0000000000",
                                   '',genes.table[[f]]$ensembl_gene_id)
}

gene_list_edited<-list()
for (f in c('seurat_wilcox','MAST','scanpy_wilcox')){
    gene_list_edited[[f]]<-genes.table[[f]]$logfc
    names(gene_list_edited[[f]])<-genes.table[[f]]$gene_id
    gene_list_edited[[f]]<-sort(gene_list_edited[[f]],decreasing = T)
}

de<-lapply(gene_list_edited,names)

reactom_enrich<-lapply(de,function(x)enrichPathway(gene=x,pvalueCutoff=0.05, readable=T))
lapply(reactom_enrich,function(x)dim(summary(x)))


summary(reactom_enrich$seurat_wilcox)

library(enrichplot)

viewPathway("Signaling by FGFR3 fusions in cancer", 
            readable = TRUE, 
            foldChange = gene_list_edited$seurat_wilcox)

cnetplot(reactom_enrich$seurat_wilcox, categorySize="pvalue", 
         foldChange=gene_list_edited$seurat_wilcox)

#Reactome pathway gene set enrichment analysis
gse <- lapply(gene_list_edited,function(x)gsePathway(x, pvalueCutoff = 0.5,pAdjustMethod = "BH", verbose = FALSE))
lapply(gse,dim)

 dotplot(gse$scanpy_wilcox, font.size=14)

# Create a gene rank based on the gene expression fold change

gene_rank<-list()

gene_rank$seurat_wilcox<- setNames(enrich_input$seurat_wilcox$avg_log2FC, casefold(rownames(enrich_input$seurat_wilcox), upper = T))
gene_rank$MAST<-setNames(enrich_input$MAST$coef,casefold(rownames(enrich_input$MAST), upper = T))
gene_rank$scanpy_wilcox<-setNames(enrich_input$scanpy_wilcox$logfoldchanges,casefold(rownames(enrich_input$scanpy_wilcox), upper = T))


fgseaRes_immune <- lapply(gene_rank,function(x){
    fgseaMultilevel(pathways = gmt, stats = x, minSize = 1, maxSize = Inf,
                           eps=0,nPermSimple = 5000)})


fgseaRes_immune_sorted <- lapply(fgseaRes_immune,function(x){x<-x[order(x$NES, decreasing = T),]})
fgseaRes_immune_overlap_pathway<-Reduce(intersect,list(fgseaRes_immune_sorted$MAST$pathway,
                                                      fgseaRes_immune_sorted$scanpy_wilcox$pathway,
                                                      fgseaRes_immune_sorted$seurat_wilcox$pathway))
# over lap pathways
fgseaRes_immune_overlap_pathway_df<-fgseaRes_immune_sorted$MAST[which(fgseaRes_immune_sorted$MAST$pathway %in% fgseaRes_immune_overlap_pathway),]
#significant pathways
fgseaRes_immune_overlap_pathway_sig<-fgseaRes_immune_overlap_pathway_df[fgseaRes_immune_overlap_pathway_df$padj<0.05,]

par(mfrow = c(1, 1), mar = c(5, 30, 2, 1))
barplot(height = fgseaRes_immune_overlap_pathway_sig$NES, names.arg = fgseaRes_immune_overlap_pathway_sig$pathway, 
    horiz = TRUE, las = 1, border = FALSE, cex.names = 0.6,xlab=list('NES',cex=1,font=2),
       main='Up/down-regulated in PD Immune')

save(fgseaRes_immune,enrich_out_immune,enrich_results_sig,sign.genes_overlap_df,fgseaRes_immune_overlap_pathway_sig,file = 'immune_DE.Rdata')

mesen<-readRDS('../Seurat/results/Mesenchymal_seurat.rds')
mesen<-SetIdent(mesen,value = 'Condition')

seu_file= mesen
fdata=as.data.frame(rownames(seu_file),row.names = rownames(seu_file))
colnames(fdata)<-'primerid'
sca<-FromMatrix(log2(as.matrix(seu_file@assays$RNA@counts)+1),seu_file@meta.data,fdata)

#load data 
load('../Seurat/results/Mesenchymal_MAST.Rdata') # 
summaryDt <- summaryCond$datatable

# Make a datatable with all the relevant output - combined Hurdle model

fcHurdle <- merge(summaryDt[contrast=='ConditionPD' & component=='H',.(primerid, `Pr(>Chisq)`)], ##hurdle P values
     summaryDt[contrast=='ConditionPD' & component=='logFC', 
               .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

# change name of Pr(>Chisq) to fdr and sort by fdr
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
#install.packages('data.table')
library(data.table)
setorder(fcHurdle, fdr)
rownames(fcHurdle)<-fcHurdle$primerid
write.table(fcHurdle,'../DE_comparisons/markers/Mesenchymal/Mesenchymal_MAST_PD_vs_Control.table',sep="\t",quote=F)



# run all DE methods
methods <- c("wilcox","bimod","t","LR","MAST","roc")
DE <- list()
for (m in methods){
    outfile <- paste0("../DE_comparisons/markers/Mesenchymal/Mesenchymal_seurat_",m,"_PD_vs_Control.table")
    DE[[m]]<- FindMarkers(object = immune,`ident.1` = 'PD',`ident.2` = 'Control',
                          test.use = m,verbose=FALSE)
    write.table(DE[[m]],file=outfile,sep="\t",quote=F)
}

methods <- c("negbinom","poisson")
for (m in methods){
    outfile <- paste0("../DE_comparisons/markers/Mesenchymal/Mesenchymal_seurat_",m,"_PD_vs_Control.table")
    DE[[m]]<- FindMarkers(object = immune,`ident.1` = 'PD',`ident.2` = 'Control',
                          test.use = m,verbose=FALSE,slot='counts')
    write.table(DE[[m]],file=outfile,sep="\t",quote=F)
}



DE[[9]]<-fcHurdle
DE[[10]]<-read.csv("../DE_comparisons/markers/Mesenchymal/Mesenchymal_Scapy_wilcoxon.csv",header = T)
names(DE)<-c("seurat_wilcox", "seurat_bimod","seurat_t","seurat_LR","seurat_MAST","seurat_roc",
             "seurat_negbinom","seurat_poisson","MAST","scanpy_wilcox")

DE$scanpy_wilcox<-DE$scanpy_wilcox[DE$scanpy_wilcox$group=='PD',]
rownames(DE$scanpy_wilcox)<-DE$scanpy_wilcox$names
rownames(DE$MAST)<-fcHurdle$primerid

names(DE)

# get top 200 genes for each test
top.200 <- lapply(DE,function(x) rownames(x)[1:200])
o <- overlap_phyper(top.200,plot=T,bg=nrow(DE$MAST))

pval.column <- c(5,5,5,5,5,5,5) # define which column contains the p-value
names(pval.column)<-names(DE)[-c(6,9,10)]
sign.genes <- list()
cutP<-0.05
for (n in names(pval.column)){
  sg <- which(DE[[n]][,pval.column[n]] < cutP)
  sign.genes[[n]]<-rownames(DE[[n]])[sg]
}

sign.genes$MAST<-DE$MAST$primerid[which(abs(DE$MAST$coef)>0.25 & DE$MAST$fdr<0.05)]
sign.genes$scanpy_wilcox<-DE$scanpy_wilcox$names[which(abs(DE$scanpy_wilcox$logfoldchanges)>0.25 & DE$scanpy_wilcox$pvals_adj<0.05)]

o <- overlap_phyper(sign.genes,plot=T,bg=nrow(DE$scanpy_wilcox))

sign.genes_overlap=Reduce(intersect, sign.genes)
sign.genes_overlap_df=DE$seurat_wilcox[sign.genes_overlap,]
#only select abs(avg_log2FC)>1
sign.genes_overlap_df_logfc1<-sign.genes_overlap_df[abs(sign.genes_overlap_df$avg_log2FC)>1,]

#up-regulated in PD
mypar(mar = c(4, 6, 3, 1))
barplot(sort(setNames(sign.genes_overlap_df_logfc1$avg_log2FC,rownames(sign.genes_overlap_df_logfc1)),F),
        horiz = T,las=1,border = "white",yaxs="i",xlab=list('avg_log2FC',cex=1, font=2),
       main='Up/down regulated genes in PD Mesen')

DotPlot(mesen, features = names(sort(setNames(sign.genes_overlap_df_logfc1$avg_log2FC,rownames(sign.genes_overlap_df_logfc1)),F)), 
        group.by = 'Condition', assay = "RNA") + coord_flip()

dbs=c('GO_Biological_Process_2017b','GO_Biological_Process_2018','KEGG_2019_Human','WikiPathways_2019_Human')
enrich_results <- enrichr(genes = sign.genes_overlap, databases = dbs)[[1]]

#significant pathways
enrich_results_sig<-enrich_results[enrich_results$Adjusted.P.value<0.05,]

par(mfrow = c(1, 1), mar = c(5, 30, 2, 1))
barplot(height = -log10(enrich_results$Adjusted.P.value)[40:1], names.arg = enrich_results$Term[40:1], 
    horiz = TRUE, las = 1, border = FALSE, cex.names = 0.7,xlab=list('-log10(Adj.P.value)',cex=1,font=2),
       main='Enrichr for sig_markers in PD Mesen')

par(mfrow = c(1, 1), mar = c(5, 30, 2, 1))
barplot(height = -log10(enrich_results$Adjusted.P.value)[80:41], names.arg = enrich_results$Term[80:41], 
    horiz = TRUE, las = 1, border = FALSE, cex.names = 0.6,xlab=list('-log10(Adj.P.value)',cex=1,font=2),
       main='Enrichr for sig_markers in PD Mesen')


# Create a gene rank based on the gene expression fold change
mesen_gene_rank <- setNames(sign.genes_overlap_df$avg_log2FC, casefold(rownames(sign.genes_overlap_df), upper = T))

# Perform enrichemnt analysis
fgseaRes_sig_markers_mesen <- fgseaMultilevel(pathways = gmt, stats = mesen_gene_rank, minSize = 1, maxSize = Inf,
                                               eps=0,nPermSimple = 5000)

fgseaRes_sig_markers_mesen <- fgseaRes_sig_markers_mesen[order(fgseaRes_sig_markers_mesen$NES, decreasing = T), ]


fgseaRes_sig_markers_mesen<-fgseaRes_sig_markers_mesen[fgseaRes_sig_markers_mesen$padj<0.05,]

par(mfrow = c(1, 1), mar = c(5, 30, 2, 1))
barplot(height = fgseaRes_sig_markers_mesen$NES[dim(fgseaRes_sig_markers_mesen)[1]:1], 
        names.arg = fgseaRes_sig_markers_mesen$pathway[dim(fgseaRes_sig_markers_mesen)[1]:1], 
    horiz = TRUE, las = 1, border = FALSE, cex.names =0.6,xlab=list('NES',cex=1,font=2),
       main='Up/down-regulated in PD Mesen')

# select significant genes for each methods: Seurat_wilcox was filtered by logfc=0.25, MAST was filted by fdr<0.05 and logfc=0.25
# scanpy_wilcox was filtered by pvals_adj and logfoldchanges=0.25
enrich_input=list()
enrich_input$seurat_wilcox<-DE$seurat_wilcox
enrich_input$MAST<-DE$MAST[which(abs(DE$MAST$coef)>0.25 & DE$MAST$fdr<0.05),]
enrich_input$scanpy_wilcox<-DE$scanpy_wilcox[abs(DE$scanpy_wilcox$logfoldchanges)>0.25,]
rownames(enrich_input$MAST)<-enrich_input$MAST$primerid

enrich_out_mesen=list()
# Perform enrichment
dbs=c('GO_Biological_Process_2017b','GO_Biological_Process_2018','KEGG_2019_Human','WikiPathways_2019_Human')

enrich_out_mesen<-lapply(enrich_input,function(x)enrichr(genes = rownames(x), databases = dbs)[[1]])

# select significant pathways
enrich_out_mesen_sig<-lapply(enrich_out_mesen,function(x)x<-x[x$Adjusted.P.value<0.05,])
#select pathways were found in all three data set
enrich_out_mesen_sig_overlap_pathway<-Reduce(intersect,list(enrich_out_mesen_sig$MAST$Term,enrich_out_mesen_sig$scanpy_wilcox$Term,enrich_out_mesen_sig$seurat_wilcox$Term))
enrich_out_mesen_sig_overlap_df<-enrich_out_mesen_sig$seurat_wilcox[which(enrich_out_mesen_sig$seurat_wilcox$Term %in% enrich_out_mesen_sig_overlap_pathway),]   


#plot 
par(mfrow = c(1, 1), mar = c(5, 30, 2, 1))
barplot(height = -log10(enrich_out_mesen_sig_overlap_df$Adjusted.P.value)[30:1], 
        names.arg = enrich_out_mesen_sig_overlap_df$Term[30:1], 
        horiz = TRUE, las = 1, border = FALSE, cex.names = 0.6,xlab=list('-log10(Adj.P.value)',cex=1,font=2),
       main='Enrichr:pathways enriched in PD Mesen')

barplot(height = -log10(enrich_out_mesen_sig_overlap_df$Adjusted.P.value)[58:31], 
        names.arg = enrich_out_mesen_sig_overlap_df$Term[58:31], 
        horiz = TRUE, las = 1, border = FALSE, cex.names = 0.5,xlab=list('-log10(Adj.P.value)',cex=1,font=2),
       main='Enrichr:pathways enriched in PD Mesen')


# Create a gene rank based on the gene expression fold change

gene_rank<-list()

gene_rank$seurat_wilcox<- setNames(enrich_input$seurat_wilcox$avg_log2FC, casefold(rownames(enrich_input$seurat_wilcox), upper = T))
gene_rank$MAST<-setNames(enrich_input$MAST$coef,casefold(rownames(enrich_input$MAST), upper = T))
gene_rank$scanpy_wilcox<-setNames(enrich_input$scanpy_wilcox$logfoldchanges,casefold(rownames(enrich_input$scanpy_wilcox), upper = T))


fgseaRes_mesen <- lapply(gene_rank,function(x){
    fgseaMultilevel(pathways = gmt, stats = x, minSize = 1, maxSize = Inf,
                           eps=0,nPermSimple = 5000)})

fgseaRes_mesen_sorted <- lapply(fgseaRes_mesen,function(x){x<-x[order(x$NES, decreasing = T),]})
fgseaRes_mesen_overlap_pathway<-Reduce(intersect,list(fgseaRes_mesen_sorted$MAST$pathway,
                                                      fgseaRes_mesen_sorted$scanpy_wilcox$pathway,
                                                      fgseaRes_mesen_sorted$seurat_wilcox$pathway))
# over lap pathways
fgseaRes_mesen_overlap_pathway_df<-fgseaRes_mesen_sorted$MAST[which(fgseaRes_mesen_sorted$MAST$pathway %in% fgseaRes_mesen_overlap_pathway),]
#significant pathways
fgseaRes_mesen_overlap_pathway_sig<-fgseaRes_mesen_overlap_pathway_df[fgseaRes_mesen_overlap_pathway_df$padj<0.05,]

par(mfrow = c(1, 1), mar = c(5, 30, 2, 1))
barplot(height = fgseaRes_mesen_overlap_pathway_sig$NES[20:1], names.arg = fgseaRes_mesen_overlap_pathway_sig$pathway[20:1], 
    horiz = TRUE, las = 1, border = FALSE, cex.names = 0.6,xlab=list('NES',cex=1,font=2),
       main='Up/down-regulated in PD Mesen')

save(fgseaRes_mesen,enrich_out_mesen,enrich_results_sig,sign.genes_overlap_df,fgseaRes_mesen_overlap_pathway_sig,file = 'mesenchymal_DE.Rdata')

neuroglial<-readRDS('../Seurat/results/NeuronGlials_seurat.rds')
neuroglial<-SetIdent(neuroglial,value = 'Condition')

seu_file= neuroglial
fdata=as.data.frame(rownames(seu_file),row.names = rownames(seu_file))
colnames(fdata)<-'primerid'
sca<-FromMatrix(log2(as.matrix(seu_file@assays$RNA@counts)+1),seu_file@meta.data,fdata)

#load data 
load('../Seurat/results/NeuroGlial_MAST.Rdata') # 
summaryDt <- summaryCond$datatable

# Make a datatable with all the relevant output - combined Hurdle model

fcHurdle <- merge(summaryDt[contrast=='ConditionPD' & component=='H',.(primerid, `Pr(>Chisq)`)], ##hurdle P values
     summaryDt[contrast=='ConditionPD' & component=='logFC', 
               .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

# change name of Pr(>Chisq) to fdr and sort by fdr
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
#install.packages('data.table')
library(data.table)
setorder(fcHurdle, fdr)
rownames(fcHurdle)<-fcHurdle$primerid
write.table(fcHurdle,'../DE_comparisons/markers/NeuronGlial/NeuronGlial_MAST_PD_vs_Control.table',sep="\t",quote=F)


# run all DE methods
methods <- c("wilcox","bimod","t","LR","MAST","roc")
DE <- list()
for (m in methods){
    outfile <- paste0("../DE_comparisons/markers/NeuronGlial/NeuronGlial_seurat_",m,"_PD_vs_Control.table")
    DE[[m]]<- FindMarkers(object = immune,`ident.1` = 'PD',`ident.2` = 'Control',
                          test.use = m,verbose=FALSE)
    write.table(DE[[m]],file=outfile,sep="\t",quote=F)
}

methods <- c("negbinom","poisson")
for (m in methods){
    outfile <- paste0("../DE_comparisons/markers/NeuronGlial/NeuronGlial_seurat_",m,"_PD_vs_Control.table")
    DE[[m]]<- FindMarkers(object = immune,`ident.1` = 'PD',`ident.2` = 'Control',
                          test.use = m,verbose=FALSE,slot='counts')
    write.table(DE[[m]],file=outfile,sep="\t",quote=F)
}


DE[[9]]<-fcHurdle
DE[[10]]<-read.csv("../DE_comparisons/markers/NeuronGlial/NeuronGlial_Scapy_wilcoxon.csv",header = T)
names(DE)<-c("seurat_wilcox", "seurat_bimod","seurat_t","seurat_LR","seurat_MAST","seurat_roc",
             "seurat_negbinom","seurat_poisson","MAST","scanpy_wilcox")

DE$scanpy_wilcox<-DE$scanpy_wilcox[DE$scanpy_wilcox$group=='PD',]
rownames(DE$scanpy_wilcox)<-DE$scanpy_wilcox$names
rownames(DE$MAST)<-fcHurdle$primerid

# get top 200 genes for each test
top.200 <- lapply(DE,function(x) rownames(x)[1:200])
o <- overlap_phyper(top.200,plot=T,bg=nrow(DE$MAST))

pval.column <- c(5,5,5,5,5,5,5) # define which column contains the p-value
names(pval.column)<-names(DE)[-c(6,9,10)]
sign.genes <- list()
cutP<-0.05
for (n in names(pval.column)){
  sg <- which(DE[[n]][,pval.column[n]] < cutP)
  sign.genes[[n]]<-rownames(DE[[n]])[sg]
}

sign.genes$MAST<-DE$MAST$primerid[which(abs(DE$MAST$coef)>0.25 & DE$MAST$fdr<0.05)]
sign.genes$scanpy_wilcox<-DE$scanpy_wilcox$names[which(abs(DE$scanpy_wilcox$logfoldchanges)>0.25 & DE$scanpy_wilcox$pvals_adj<0.05)]

o <- overlap_phyper(sign.genes,plot=T,bg=nrow(DE$MAST))

sign.genes_overlap=Reduce(intersect, sign.genes)
sign.genes_overlap_df=DE$seurat_wilcox[sign.genes_overlap,]
#only select abs(avg_log2FC)>1
sign.genes_overlap_df_logfc1<-sign.genes_overlap_df[abs(sign.genes_overlap_df$avg_log2FC)>1,]

mypar(mar = c(4, 6, 3, 1))
barplot(sort(setNames(sign.genes_overlap_df_logfc1$avg_log2FC,rownames(sign.genes_overlap_df_logfc1)),F),
        horiz = T,las=1,border = "white",yaxs="i",xlab=list('avg_log2FC',cex=1, font=2),
       main='Up/down regulated genes in PD NeuroGlial')

DotPlot(mesen, features = names(sort(setNames(sign.genes_overlap_df_logfc1$avg_log2FC,rownames(sign.genes_overlap_df_logfc1)),F)), 
        group.by = 'Condition', assay = "RNA") + coord_flip()

dbs=c('GO_Biological_Process_2017b','GO_Biological_Process_2018','KEGG_2019_Human','WikiPathways_2019_Human')
enrich_results <- enrichr(genes = sign.genes_overlap, databases = dbs)[[1]]

#significant pathways
enrich_results_sig<-enrich_results[enrich_results$Adjusted.P.value<0.05,]

par(mfrow = c(1, 1), mar = c(5, 30, 2, 1))
barplot(height = -log10(enrich_results$Adjusted.P.value)[11:1], names.arg = enrich_results$Term[11:1], 
    horiz = TRUE, las = 1, border = FALSE, cex.names =0.9,xlab=list('-log10(Adj.P.value)',cex=1,font=2),
       main='Enrichr for sig_markers in PD NeuroGlial')


# Create a gene rank based on the gene expression fold change
neuroglial_gene_rank <- setNames(sign.genes_overlap_df$avg_log2FC, casefold(rownames(sign.genes_overlap_df), upper = T))

# Perform enrichemnt analysis
fgseaRes_sig_markers_neuroglial <- fgseaMultilevel(pathways = gmt, stats = neuroglial_gene_rank, minSize = 1, maxSize = Inf,
                                               eps=0,nPermSimple = 5000)

fgseaRes_sig_markers_neuroglial <- fgseaRes_sig_markers_neuroglial[order(fgseaRes_sig_markers_neuroglial$NES, decreasing = T), ]


fgseaRes_sig_markers_neuroglial<-fgseaRes_sig_markers_neuroglial[fgseaRes_sig_markers_neuroglial$padj<0.05,]

par(mfrow = c(1, 1), mar = c(5, 30, 2, 1))
barplot(height = fgseaRes_sig_markers_neuroglial$NES[dim(fgseaRes_sig_markers_neuroglial)[1]:1], 
        names.arg = fgseaRes_sig_markers_neuroglial$pathway[dim(fgseaRes_sig_markers_neuroglial)[1]:1], 
    horiz = TRUE, las = 1, border = FALSE, cex.names =0.9,xlab=list('NES',cex=1,font=2),
       main='Up/down-regulated in PD NeuroGlial')

# select significant genes for each methods: Seurat_wilcox was filtered by logfc=0.25, MAST was filted by fdr<0.05 and logfc=0.25
# scanpy_wilcox was filtered by pvals_adj and logfoldchanges=0.25
enrich_input=list()
enrich_input$seurat_wilcox<-DE$seurat_wilcox
enrich_input$MAST<-DE$MAST[which(abs(DE$MAST$coef)>0.25 & DE$MAST$fdr<0.05),]
enrich_input$scanpy_wilcox<-DE$scanpy_wilcox[abs(DE$scanpy_wilcox$logfoldchanges)>0.25,]
rownames(enrich_input$MAST)<-enrich_input$MAST$primerid

enrich_out_neuroglial=list()
# Perform enrichment
dbs=c('GO_Biological_Process_2017b','GO_Biological_Process_2018','KEGG_2019_Human','WikiPathways_2019_Human')

enrich_out_neuroglial<-lapply(enrich_input,function(x)enrichr(genes = rownames(x), databases = dbs)[[1]])

# select significant pathways
enrich_out_neuroglial_sig<-lapply(enrich_out_neuroglial,function(x)x<-x[x$Adjusted.P.value<0.05,])
#select pathways were found in all three data set
enrich_out_neuroglial_sig_overlap_pathway<-Reduce(intersect,list(enrich_out_neuroglial_sig$MAST$Term,enrich_out_neuroglial_sig$scanpy_wilcox$Term,enrich_out_neuroglial_sig$seurat_wilcox$Term))
enrich_out_neuroglial_sig_overlap_df<-enrich_out_neuroglial_sig$seurat_wilcox[which(enrich_out_neuroglial_sig$seurat_wilcox$Term %in% enrich_out_neuroglial_sig_overlap_pathway),]   



enrich_out_neuroglial_sig_overlap_df

# Create a gene rank based on the gene expression fold change

gene_rank<-list()

gene_rank$seurat_wilcox<- setNames(enrich_input$seurat_wilcox$avg_log2FC, casefold(rownames(enrich_input$seurat_wilcox), upper = T))
gene_rank$MAST<-setNames(enrich_input$MAST$coef,casefold(rownames(enrich_input$MAST), upper = T))
gene_rank$scanpy_wilcox<-setNames(enrich_input$scanpy_wilcox$logfoldchanges,casefold(rownames(enrich_input$scanpy_wilcox), upper = T))


fgseaRes_neuroglial <- lapply(gene_rank,function(x){
    fgseaMultilevel(pathways = gmt, stats = x, minSize = 1, maxSize = Inf,
                           eps=0,nPermSimple = 50000)})

fgseaRes_neuroglial_sorted <- lapply(fgseaRes_neuroglial,function(x){x<-x[order(x$NES, decreasing = T),]})
fgseaRes_neuroglial_overlap_pathway<-Reduce(intersect,list(fgseaRes_neuroglial_sorted$MAST$pathway,
                                                      fgseaRes_neuroglial_sorted$scanpy_wilcox$pathway,
                                                      fgseaRes_neuroglial_sorted$seurat_wilcox$pathway))
# over lap pathways
fgseaRes_neuroglial_overlap_pathway_df<-fgseaRes_neuroglial_sorted$MAST[which(fgseaRes_neuroglial_sorted$MAST$pathway %in% fgseaRes_neuroglial_overlap_pathway),]
#significant pathways
fgseaRes_neuroglial_overlap_pathway_sig<-fgseaRes_neuroglial_overlap_pathway_df[fgseaRes_neuroglial_overlap_pathway_df$padj<0.05,]

plot_df1<-head(fgseaRes_neuroglial_overlap_pathway_sig,n=40)
plot_df2<-tail(fgseaRes_neuroglial_overlap_pathway_sig,n=40)

par(mfrow = c(1, 1), mar = c(5, 30, 2, 1))
barplot(height = plot_df1$NES[40:1], names.arg = plot_df1$pathway[40:1], 
    horiz = TRUE, las = 1, border = FALSE, cex.names = 0.6,xlab=list('NES',cex=1,font=2),
       main='Upregulated in PD NeuroGlial')

par(mfrow = c(1, 1), mar = c(5, 30, 2, 1))
barplot(height = plot_df2$NES, names.arg = plot_df2$pathway, 
    horiz = TRUE, las = 1, border = FALSE, cex.names = 0.6,xlab=list('NES',cex=1,font=2),
       main='Downregulate in PD NeuroGlial')


save(fgseaRes_neuroglial,enrich_out_neuroglial,enrich_results_sig,sign.genes_overlap_df,fgseaRes_neuroglial_overlap_pathway_sig,file = 'neuroglial_DE.Rdata')

library(ReactomePA)

data(geneList)










































