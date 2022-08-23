suppressPackageStartupMessages({
  library(Seurat)
  library(venn)
  library(dplyr)
  library(cowplot)
  library(ggplot2)
  library(gridExtra)
  library(CellChat)
  library(patchwork)
  library(NMF)
  library(ggalluvial)
  library(Matrix)
})


cellchat.ms<-readRDS('../output/LV_MS.rds')
cellchat.ctl<-readRDS('../output/LV_control.rds')
object.list<-list(LV_Control=cellchat.ctl,LV_MS=cellchat.ms)
cellchat<-mergeCellChat(object.list,add.names = names(object.list))

rm(cellchat.ctl,cellchat.ms)
gc()

options(rePercentageFeatureSet.width = 4, repr.plot.height = 4, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
gg1<-compareInteractions(cellchat,show.legend = F,group = c(1,2),size.text = 14)
gg2<-compareInteractions(cellchat,show.legend = F,group = c(1,2),measure = "weight",size.text = 14)
print(gg1+gg2)

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
par(mfrow = c(1,1))
netVisual_diffInteraction(cellchat,weight.scale = T,measure = "count",
                          title.name = "Differential No. interactions in MS vs Control")
netVisual_diffInteraction(cellchat,weight.scale = T,measure = "weight",
                          title.name = "Differential No. interactions strength in MS vs Control")

gg1<-netVisual_heatmap(cellchat,title.name = "Differential No.interactions in MS vs Control ",font.size.title = 8)
gg2<-netVisual_heatmap(cellchat,measure = "weight",title.name = "Differential No.interactions strength in MS vs Control",font.size.title = 8)
gg1

gg2

#total counts of interactions
weight.max<-getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow=c(1,2))
options(rePercentageFeatureSet.width = 10, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count,weight.scale = T,label.edge = F,
                   edge.weight.max = weight.max[2],edge.width.max = 12,
                   vertex.weight =as.numeric(table(object.list[[i]]@idents)),
                   vertex.weight.max = weight.max[1],vertex.size.max = 12,
                   title.name = paste0("Number of interactions-",names(object.list)[i]))
}


weight.max<-getMaxWeight(object.list, attribute = c("idents","weight"))
options(rePercentageFeatureSet.width = 10, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
par(mfrow=c(1,2))
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight,weight.scale = T,label.edge = F,
                   edge.weight.max = weight.max[2],edge.width.max = 12,
                   vertex.weight =as.numeric(table(object.list[[i]]@idents)),
                   vertex.weight.max = weight.max[1],vertex.size.max = 12,
                   title.name = paste0("Number of strength-",names(object.list)[i]))
}

#interaction strength from each cell group to others
mat_ctl <- cellchat@net$LV_Control$weight
mat_ms<-cellchat@net$LV_MS$weight
groupSize_ctl <- as.numeric(table(cellchat@idents$LV_Control))
groupSize_ms <- as.numeric(table(cellchat@idents$LV_MS))

weight.max<-getMaxWeight(object.list, attribute = c("idents","weight"))

options(rePercentageFeatureSet.width = 8, repr.plot.height = 4, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
for (i in 1:nrow(mat_ctl)) {
  par(mfrow = c(1,2))# xpd=TRUE
  mat1 <- matrix(0, nrow = nrow(mat_ctl), ncol = ncol(mat_ctl), dimnames = dimnames(mat_ctl))
  mat1[i, ] <- mat_ctl[i, ]
  netVisual_circle(mat1, 
                         weight.scale = T, edge.weight.max = weight.max[2], 
                         edge.width.max = 12,label.edge = T,margin = 0,
                         vertex.weight =as.numeric(table(object.list[[1]]@idents)),
                         vertex.weight.max = weight.max[1],vertex.size.max = 12,
                         title.name = paste0(rownames(mat_ctl)[i],"_LV_Control"))
  mat2 <- matrix(0, nrow = nrow(mat_ms), ncol = ncol(mat_ms), dimnames = dimnames(mat_ms))
  mat2[i, ] <- mat_ms[i, ]
  netVisual_circle(mat2,
                         weight.scale = T, edge.weight.max = weight.max[2], 
                         edge.width.max = 12,label.edge = T,margin = 0,
                         vertex.weight =as.numeric(table(object.list[[2]]@idents)),
                         vertex.weight.max = weight.max[1],vertex.size.max = 12,
                         title.name = paste0(rownames(mat_ms)[i],"_LV_MS"))
}

options(rePercentageFeatureSet.width = 12, repr.plot.height = 5, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
num.link <- sapply(object.list, function(x) {
  rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)
})

weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], 
                                               title = names(object.list)[i], 
                                               weight.MinMax = weight.MinMax)
}

patchwork::wrap_plots(plots = gg)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Endothelial")
gg1

gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Mesenchymal")
gg2


gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Immune")
gg3



gg4 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Epithelial")
gg4

# Identify signaling groups based on their functional similarity
cellchat<-computeNetSimilarityPairwise(cellchat,type = "functional")
cellchat<-netEmbedding(cellchat,type = "functional")
cellchat<-netClustering(cellchat,type = "functional")
options(rePercentageFeatureSet.width = 5, repr.plot.height = 5, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
netVisual_embeddingPairwise(cellchat,type = "functional",label.size = 2)

options(rePercentageFeatureSet.width = 8, repr.plot.height = 8, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2)

# Identify signaling groups based on structure similarity
cellchat<-computeNetSimilarityPairwise(cellchat,type = "structural")
cellchat<-netEmbedding(cellchat,type = "structural")
cellchat<-netClustering(cellchat,type = "structural")

options(rePercentageFeatureSet.width = 8, repr.plot.height = 8, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
netVisual_embeddingPairwise(cellchat,type = "structural",label.size = 2)

options(rePercentageFeatureSet.width = 20, repr.plot.height = 8, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 3)

options(rePercentageFeatureSet.width = 14, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
par(mfrow=c(1,2))
gg1<-rankSimilarity(cellchat,type = "functional",font.size = 10,
               title ="Based on functional similarity")
gg2<-rankSimilarity(cellchat,type = "structural",font.size = 10,
               title ="Based on structure similarity")
gg1+gg2

par(mfrow=c(1,2))
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

library(ComplexHeatmap)
# combining all the identified signaling pathways from different datasets
i=1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)

# outgoing signaling
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union,
                                        title = names(object.list)[i], width = 6, height = 12,font.size = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, 
                                        title = names(object.list)[i+1], width = 6, height = 12,font.size = 6)
options(rePercentageFeatureSet.width = 8, repr.plot.height = 8, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
ht1

ht2

# incoming signaling
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union,
                                        title = names(object.list)[i], width = 6, height = 12, font.size = 6,
                                        color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, 
                                        title = names(object.list)[i+1], width = 6, height = 12, font.size = 6,
                                        color.heatmap = "GnBu")
ht1


ht2

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, 
                                        title = names(object.list)[i], width = 6, height = 12, 
                                        color.heatmap = "OrRd",font.size = 8)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, 
                                        title = names(object.list)[i+1], width = 6, height = 12, 
                                        color.heatmap = "OrRd",font.size = 8)
ht1

ht2

options(rePercentageFeatureSet.width = 6, repr.plot.height = 10, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
netVisual_bubble(cellchat,sources.use = c(1),targets.use = c(2:4),comparison = c(1,2),angle.x = 45)
netVisual_bubble(cellchat,sources.use = c(2),targets.use = c(1,3,4),comparison = c(1,2),angle.x = 45)
netVisual_bubble(cellchat,sources.use = c(4),targets.use = c(1,2,3),comparison = c(1,2),angle.x = 45)
options(rePercentageFeatureSet.width = 5, repr.plot.height = 5, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
netVisual_bubble(cellchat,sources.use = c(3),targets.use = c(1,2,4),comparison = c(1,2),angle.x = 45)

netVisual_bubble(cellchat, sources.use = c(1), targets.use = c(2:4), comparison = c(1, 2), max.dataset = 2, 
                        title.name = "Increased signaling in MS", angle.x = 45, remove.isolate = T)
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2:4),  comparison = c(1,2), max.dataset = 1, 
                        title.name = "Decreased signaling in MS", angle.x = 45, remove.isolate = T)

options(rePercentageFeatureSet.width = 5, repr.plot.height = 8, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1,3,4), comparison = c(1, 2), max.dataset = 2, 
                        title.name = "Increased signaling in MS", angle.x = 45, remove.isolate = T)
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1,3,4),  comparison = c(1,2), max.dataset = 1, 
                        title.name = "Decreased signaling in MS", angle.x = 45, remove.isolate = T)

options(rePercentageFeatureSet.width = 5, repr.plot.height = 3, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1,2,4), comparison = c(1, 2), max.dataset = 2, 
                        title.name = "Increased signaling in MS", angle.x = 45, remove.isolate = T)
netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1,2,4),  comparison = c(1,2), max.dataset = 1, 
                        title.name = "Decreased signaling in MS", angle.x = 45, remove.isolate = T)

options(rePercentageFeatureSet.width = 5, repr.plot.height = 8, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1,2,3), comparison = c(1, 2), max.dataset = 2, 
                        title.name = "Increased signaling in MS", angle.x = 45, remove.isolate = T)
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1,2,3),  comparison = c(1,2), max.dataset = 1, 
                        title.name = "Decreased signaling in MS", angle.x = 45, remove.isolate = T)

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "LV_MS"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, 
                                       features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, 
                                       thresh.fc = 0.1, thresh.p = 0.05)
# map the results of differential expression analysis onto the inferred cell-cell communications
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in MS
net.up <- subsetCommunication(cellchat, net = net, datasets = "LV_MS",ligand.logFC = 0.2, receptor.logFC = 0.2)
# extract the ligand-receptor pairs with down-regulated ligands and upregulated recetptors in MS
net.down <- subsetCommunication(cellchat, net = net, datasets = "LV_Control",ligand.logFC = (-0.2), receptor.logFC =(-0.2))

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)


#up-regulated in MS
options(rePercentageFeatureSet.width = 5, repr.plot.height = 3, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pairLR.use.up = net.up[, "interaction_name", drop = F]
netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(1,2), 
                 targets.use = c(3,4), comparison = c(1,2), angle.x = 45, remove.isolate = T,
                        title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

#down-regulated in MS
pairLR.use.down = net.down[, "interaction_name", drop = F]
netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(1,2), 
                 targets.use = c(3,4), comparison = c(1, 2), angle.x = 45, remove.isolate = T,
                 title.name = paste0("Down-regulated signaling in ",names(object.list)[2]))

#identify the same pathways between two datasets. 
pathway.commom<-intersect(cellchat@netP$LV_Control$pathways,cellchat@netP$LV_MS$pathways)
options(rePercentageFeatureSet.width = 12, repr.plot.height = 5, repr.plot.res = 200,ggrepel.max.overlaps = Inf)

for (i in 1:length(pathway.commom)) {
    pathways.show<-pathway.commom[i]
    weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute =pathways.show)
    par(mfrow = c(1,2),xpd=TRUE)
    #circle plot
    for (i in 1:length(object.list)) {
        netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle",
                            edge.weight.max = weight.max[1], edge.width.max = 10,
                            signaling.name = paste0(pathways.show,' ', names(object.list)[i]))
    }
}


pathway.commom<-intersect(cellchat@netP$LV_Control$pathways,cellchat@netP$LV_MS$pathways)
pathway.commom

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("COLLAGEN") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("MK") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("APP") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("PTN") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("VEGF") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("LAMININ") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("PECAM1") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("ESAM") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("CDH5") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("ANGPTL") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("THBS") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("MIF") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("CXCL") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("HSPG") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("SPP1") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("PDGF") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("GAS") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("SEMA3") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("JAM") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("SEMA7") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("CD45") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("VTN") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("IGF") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("BMP") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("ANGPT") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("CALCR") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("GRN") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("PROS") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("CLDN") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("LIFR") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("TWEAK") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("CD39") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

options(rePercentageFeatureSet.width = 6, repr.plot.height = 6, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
pathways.show <- c("OCLN") 
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),font.size = 12)
}
ht[[1]]

ht[[2]]

# Chord diagram
options(rePercentageFeatureSet.width = 10, repr.plot.height = 8, repr.plot.res = 200,ggrepel.max.overlaps = Inf)

for (i in 1:length(levels(cellchat@idents$LV_Control))) {
  for (j in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[j]], sources.use = i, 
                       targets.use = seq(1:length(levels(cellchat@idents$LV_Control)))[-i], lab.cex = 0.4, 
                      title.name = paste0("Signaling from ",levels(cellchat@idents$LV_MS)[i],' ',names(object.list)[j]),legend.pos.x=8,legend.pos.y=30)
  }
}

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("LV_Control", "LV_MS"))


pathway.all<-unique(c(cellchat@netP$LV_MS$pathways,cellchat@netP$LV_Control$pathways))
for (i in 1:length(pathway.all)) {
  print(plotGeneExpression(cellchat, signaling = pathway.all[i],split.by='datasets',colors.ggplot=TRUE)) #default only plot significant communications
}

saveRDS(cellchat,'../output/LV_control_MS.rds')

sessionInfo()












