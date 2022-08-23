setwd("~/OneDrive - Karolinska Institutet/Mac/Documents/snRNA/MS/mouse/ms_model_LV_SLV/cellChat/")
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

lv <-readRDS('input/lv_cellchat_input.rds')


lv_MS <- subset(lv,cells = which(lv$Condition=='MS'))

dim(lv_MS)

#normalize data
lv_MS<-NormalizeData(lv_MS)
#extract normalized count matrix from seu object
data.input<-GetAssayData(lv_MS,assay = "RNA",slot="data")
groups<-lv_MS$CellType_correct_edit1
meta<-data.frame(group=groups,row.names = names(groups))

rm(lv)
gc()

cellchat <- createCellChat(object = data.input, meta = lv_MS@meta.data, group.by = "CellType_correct_edit1")
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents)) 
rm(lv_MS,data.input)
gc()

options(rePercentageFeatureSet.width = 4, repr.plot.height = 4, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use

unique(cellchat@DB$interaction$pathway_name)

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) 

# remove H2-BI and H2-Ea-ps
print(which(CellChatDB.use[["interaction"]]$ligand == "H2-BI")) # 1887
CellChatDB.use[["interaction"]] <- CellChatDB.use[["interaction"]][-1887,]
print(which(CellChatDB.use[["interaction"]]$ligand == "H2-Ea-ps")) #1900
CellChatDB.use[["interaction"]] <- CellChatDB.use[["interaction"]][-1900,]

#future::plan("multiprocess", workers = 4) # do parallel
#Identify over-expressed signaling genes associated with each cell group
cellchat <- identifyOverExpressedGenes(cellchat)
#Identify over-expressed ligand-receptor interactions (pairs)
cellchat <- identifyOverExpressedInteractions(cellchat)

head(cellchat@var.features$features.info) 
#pct.1: percentage of average feature expression across cluster1, pct.2: percentage average expression across other clusters

cellchat <- computeCommunProb(cellchat, raw.use = TRUE) # with updated slot 'net'
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)
par(mfrow = c(1,2),xpd=TRUE)#,
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

options(rePercentageFeatureSet.width = 8, repr.plot.height = 8, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
mat <- cellchat@net$weight
par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

options(rePercentageFeatureSet.width = 4, repr.plot.height = 4, repr.plot.res = 200,ggrepel.max.overlaps = Inf)

for (i in 1:length(pathways.show.all)){
    netVisual_aggregate(cellchat, signaling = pathways.show.all[i],  vertex.receiver = vertex.receiver,layout='circle')
}

options(rePercentageFeatureSet.width = 5, repr.plot.height = 5, repr.plot.res = 200,ggrepel.max.overlaps = Inf)

for (i in 1:length(pathways.show.all)){
    netVisual_aggregate(cellchat, signaling = pathways.show.all[i],  vertex.receiver = vertex.receiver,layout='chord')
}


options(rePercentageFeatureSet.width = 5, repr.plot.height = 5, repr.plot.res = 200,ggrepel.max.overlaps = Inf)

for (i in 1:length(pathways.show.all)){
    gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i],width = 0.07)
    print(gg)
}

options(rePercentageFeatureSet.width = 8, repr.plot.height = 10, repr.plot.res = 200,ggrepel.max.overlaps = Inf)

cell_cluster<-levels(cellchat@idents)
vec<-1:length(cell_cluster)
for (i in 1:length(cell_cluster)) {
  print(netVisual_bubble(cellchat, sources.use = i, targets.use =vec[-i], remove.isolate = FALSE))
}

saveRDS(cellchat,'output/LV_MS.rds')

cellchat <-readRDS('output/LV_MS.rds')

options(rePercentageFeatureSet.width = 8, repr.plot.height = 8, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
for (i in 1:length(levels(cellchat@idents))) {
  print(netVisual_chord_gene(cellchat, sources.use = i, targets.use = seq(1,4)[-i], lab.cex = 0.4,legend.pos.y = 30))
}

for (i in 1:length(levels(cellchat@idents))) {
  print(netVisual_chord_gene(cellchat, sources.use = i, targets.use = seq(1,4)[-i],
                             slot.name = "netP",legend.pos.x = 10,lab.cex = 0.6))
}

pathways.show.all <- cellchat@netP$pathways
for (i in 1:length(pathways.show.all)) {
  print(plotGeneExpression(cellchat, signaling = pathways.show.all[i],
                           group.by = cellchat@idents,type = "violin")) #default only plot significant communications
}


cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
print(netAnalysis_signalingRole_network(cellchat,  width = 8, height = 6, font.size = 10))


options(rePercentageFeatureSet.width = 5, repr.plot.height = 5, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
netAnalysis_signalingRole_scatter(cellchat,title = "all_signal_pathway")

options(rePercentageFeatureSet.width = 5, repr.plot.height = 5, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
for (i in 1:length(pathways.show.all)) {
  print(netAnalysis_signalingRole_scatter(cellchat, signaling = pathways.show.all[i],
                                          title = pathways.show.all[i]))
}

options(rePercentageFeatureSet.width = 5, repr.plot.height = 8, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",font.size=8,width = 6,height = 14)
ht1

ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",font.size=8,width = 6,height = 14)
ht2

library(ggalluvial)
options(rePercentageFeatureSet.width = 5, repr.plot.height = 5, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
selectK(cellchat, pattern = "outgoing")

nPatterns = 4 #Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 3.
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns,font.size = 7,width = 6,height = 10)
netAnalysis_river(cellchat, pattern = "outgoing")
netAnalysis_dot(cellchat, pattern = "outgoing")

options(rePercentageFeatureSet.width = 5, repr.plot.height = 5, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
selectK(cellchat, pattern = "incoming")

nPatterns = 4 #Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 3.
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns,font.size = 7,width = 6,height = 10)
netAnalysis_river(cellchat, pattern = "incoming")
netAnalysis_dot(cellchat, pattern = "incoming")


cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
(grep("NaN",cellchat@netP$similarity))

cellchat <- netClustering(cellchat, type = "functional")
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)

cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")

options(rePercentageFeatureSet.width = 8, repr.plot.height = 8, repr.plot.res = 200,ggrepel.max.overlaps = Inf)
netVisual_embedding(cellchat, type = "structural", label.size = 2.5,
                    title = "Structure similarity")
netVisual_embeddingZoomIn(cellchat, type = "structural",nCol = 2)



saveRDS(cellchat, file = "output/LV_MS.rds")

sessionInfo()






