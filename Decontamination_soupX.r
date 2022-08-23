library(SoupX)
library(Matrix)

Epithelial_marker = c("TTR","CLIC6","FOLR1","IGFBP2","PDGFA")
Endothelial_marker = c("PLVAP","CD34","ICAM1","CLDN5")
Immune_marker = c("AIF1","TYROBP","FCER1G","CD74",'CD68')
Mesenchyma_marker = c("LUM","COL1A1","COL1A2","COL3A1","CD44")
T_cells_NK = c("SKAP1","PTPRC","MBNL1","B2M","CDC42SE2")
Astrocyte = c('GFAP')
Oligos = c("FA2H","CNP","UGT8")
Neuros = c("KCNIP4","NRXN1","KCND2","RIMS1","RIMS2","DSCAM","TNR","VCAN")
OPCs = c('PDGFRA', 'SOX6','SOX10')

toc = Seurat::Read10X('../../../../raw_counts_scRNA/xxxx/outs/filtered_feature_bc_matrix/')
tod = Seurat::Read10X('../../../../raw_counts_scRNA/xxxx/outs/raw_feature_bc_matrix/')
sc = SoupChannel(tod, toc)

clt = read.csv('../../../../raw_counts_scRNA/xxxx/outs/analysis/clustering/kmeans_6_clusters/clusters.csv',header = T)
tsne = read.csv('../../../../raw_counts_scRNA/xxxx/outs/analysis/tsne/2_components/projection.csv',header = T)
rownames(tsne) = tsne$Barcode

sc = setClusters(sc, setNames(clt$Cluster, rownames(clt)))
sc = setDR(sc, tsne[colnames(sc$toc), c("TSNE.1", "TSNE.2")])

#clustering K-mean=6
library(ggplot2)
dd = sc$metaData[colnames(sc$toc), ]
mids = aggregate(cbind(TSNE.1, TSNE.2) ~ clusters, data = dd, FUN = mean)
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = clusters), size = 0.2) + 
    geom_label(data = mids, aes(label = clusters)) + ggtitle("Spain_18_006_Ctrl_M") + 
    guides(colour = guide_legend(override.aes = list(size = 1)))
plot(gg)

#epithelial
dd$CLIC6 = sc$toc["CLIC6", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = CLIC6 > 0))
plot(gg)

#mesenchyal
dd$PLVAP = sc$toc["PLVAP", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = PLVAP > 0))
plot(gg)
#endothelial
dd$COL3A1 = sc$toc["COL3A1", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = COL3A1 > 0))
plot(gg)
#immune C1qa
dd$SKAP1 = sc$toc["SKAP1", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = SKAP1 > 0))
plot(gg)
#neuroglial Cnp
dd$CNP = sc$toc["CNP", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = CNP > 0))
plot(gg)

sc = autoEstCont(sc)

out = adjustCounts(sc)

plotChangeMap(sc, out, "CLIC6")
plotChangeMap(sc, out, "PLVAP")
plotChangeMap(sc, out, "COL3A1")
plotChangeMap(sc, out, "SKAP1")
plotChangeMap(sc, out, "CNP")

DropletUtils:::write10xCounts("../../../../raw_counts_scRNA/spain_human/AU7003-18-006-SI-GA-C12/outs/strainedCounts", out)

toc = Seurat::Read10X('../../../../raw_counts_scRNA/spain_human/AU7006-19-005-SI-GA-D3/outs/filtered_feature_bc_matrix/')
tod = Seurat::Read10X('../../../../raw_counts_scRNA/spain_human/AU7006-19-005-SI-GA-D3/outs/raw_feature_bc_matrix/')
sc = SoupChannel(tod, toc)

clt = read.csv('../../../../raw_counts_scRNA/spain_human/AU7006-19-005-SI-GA-D3/outs/analysis/clustering/kmeans_6_clusters/clusters.csv',header = T)
tsne = read.csv('../../../../raw_counts_scRNA/spain_human/AU7006-19-005-SI-GA-D3/outs/analysis/tsne/2_components/projection.csv',header = T)
rownames(tsne) = tsne$Barcode

sc = setClusters(sc, setNames(clt$Cluster, rownames(clt)))
sc = setDR(sc, tsne[colnames(sc$toc), c("TSNE.1", "TSNE.2")])

#clustering K-mean=6
library(ggplot2)
dd = sc$metaData[colnames(sc$toc), ]
mids = aggregate(cbind(TSNE.1, TSNE.2) ~ clusters, data = dd, FUN = mean)
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = clusters), size = 0.2) + 
    geom_label(data = mids, aes(label = clusters)) + ggtitle("Spain_19_005_Ctrl_M") + 
    guides(colour = guide_legend(override.aes = list(size = 1)))
plot(gg)

#epithelial
dd$CLIC6 = sc$toc["CLIC6", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = CLIC6 > 0))
plot(gg)

#mesenchyal
dd$PLVAP = sc$toc["PLVAP", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = PLVAP > 0))
plot(gg)
#endothelial
dd$COL3A1 = sc$toc["COL3A1", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = COL3A1 > 0))
plot(gg)
#immune C1qa
dd$SKAP1 = sc$toc["SKAP1", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = SKAP1 > 0))
plot(gg)
#neuroglial Cnp
dd$CNP = sc$toc["CNP", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = CNP > 0))
plot(gg)

sc = autoEstCont(sc)

out = adjustCounts(sc)

plotChangeMap(sc, out, "CLIC6")
plotChangeMap(sc, out, "PLVAP")
plotChangeMap(sc, out, "COL3A1")
plotChangeMap(sc, out, "SKAP1")
plotChangeMap(sc, out, "CNP")

DropletUtils:::write10xCounts("../../../../raw_counts_scRNA/spain_human/AU7006-19-005-SI-GA-D3/outs/strainedCounts", out)

toc = Seurat::Read10X('../../../../raw_counts_scRNA/10X_20_AF_09_R77/10X_20_061_premrna/outs/filtered_feature_bc_matrix/')
tod = Seurat::Read10X('../../../../raw_counts_scRNA/10X_20_AF_09_R77/10X_20_061_premrna/outs/raw_feature_bc_matrix/')
sc = SoupChannel(tod, toc)

clt = read.csv('../../../../raw_counts_scRNA/10X_20_AF_09_R77/10X_20_061_premrna/outs/analysis/clustering/kmeans_6_clusters/clusters.csv',header = T)
tsne = read.csv('../../../../raw_counts_scRNA/10X_20_AF_09_R77/10X_20_061_premrna/outs/analysis/tsne/2_components/projection.csv',header = T)
rownames(tsne) = tsne$Barcode

sc = setClusters(sc, setNames(clt$Cluster, rownames(clt)))
sc = setDR(sc, tsne[colnames(sc$toc), c("TSNE.1", "TSNE.2")])

#clustering K-mean=6
library(ggplot2)
dd = sc$metaData[colnames(sc$toc), ]
mids = aggregate(cbind(TSNE.1, TSNE.2) ~ clusters, data = dd, FUN = mean)
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = clusters), size = 0.2) + 
    geom_label(data = mids, aes(label = clusters)) + ggtitle("R77_61_Ctrl_F") + 
    guides(colour = guide_legend(override.aes = list(size = 1)))
plot(gg)

#epithelial
dd$CLIC6 = sc$toc["CLIC6", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = CLIC6 > 0))
plot(gg)

#mesenchyal
dd$CLDN5 = sc$toc["PLVAP", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = CLDN5 > 0))
plot(gg)
#endothelial
dd$COL3A1 = sc$toc["COL3A1", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = COL3A1 > 0))
plot(gg)
#immune C1qa
dd$SKAP1 = sc$toc["SKAP1", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = SKAP1 > 0))
plot(gg)
#neuroglial Cnp
dd$CNP = sc$toc["CNP", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = CNP > 0))
plot(gg)

sc = autoEstCont(sc)

out = adjustCounts(sc)

plotChangeMap(sc, out, "CLIC6")
plotChangeMap(sc, out, "PLVAP")
plotChangeMap(sc, out, "COL3A1")
plotChangeMap(sc, out, "SKAP1")
plotChangeMap(sc, out, "CNP")

DropletUtils:::write10xCounts("../../../../raw_counts_scRNA/10X_20_AF_09_R77/10X_20_061_premrna/outs/strainedCounts", out)

toc = Seurat::Read10X('../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_063_new/outs/filtered_feature_bc_matrix/')
tod = Seurat::Read10X('../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_063_new/outs/raw_feature_bc_matrix/')
sc = SoupChannel(tod, toc)

clt = read.csv('../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_063_new/outs/analysis/clustering/kmeans_6_clusters/clusters.csv',header = T)
tsne = read.csv('../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_063_new/outs/analysis/tsne/2_components/projection.csv',header = T)
rownames(tsne) = tsne$Barcode

sc = setClusters(sc, setNames(clt$Cluster, rownames(clt)))
sc = setDR(sc, tsne[colnames(sc$toc), c("TSNE.1", "TSNE.2")])

#clustering K-mean=6
library(ggplot2)
dd = sc$metaData[colnames(sc$toc), ]
mids = aggregate(cbind(TSNE.1, TSNE.2) ~ clusters, data = dd, FUN = mean)
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = clusters), size = 0.2) + 
    geom_label(data = mids, aes(label = clusters)) + ggtitle("R100_63_Ctrl_F") + 
    guides(colour = guide_legend(override.aes = list(size = 1)))
plot(gg)

#epithelial
dd$CLIC6 = sc$toc["CLIC6", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = CLIC6 > 0))
plot(gg)

#mesenchyal
dd$CLDN5 = sc$toc["CLDN5", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = CLDN5 > 0))
plot(gg)
#endothelial
dd$COL3A1 = sc$toc["COL3A1", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = COL3A1 > 0))
plot(gg)
#immune C1qa
dd$SKAP1 = sc$toc["SKAP1", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = SKAP1 > 0))
plot(gg)
#neuroglial Cnp
dd$CNP = sc$toc["CNP", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = CNP > 0))
plot(gg)

clt$CellType=clt$Cluster
clt$CellType = gsub('1','Epihelial',clt$CellType)
clt$CellType = gsub('2','Epihelial',clt$CellType)
clt$CellType = gsub('3','Epihelial',clt$CellType)
clt$CellType = gsub('4','Mesenchymal',clt$CellType)
clt$CellType = gsub('5','Mesenchymal',clt$CellType)
clt$CellType = gsub('6','Immune',clt$CellType)

dd$CellType=clt$CellType
mids = aggregate(cbind(TSNE.1, TSNE.2) ~ CellType, data = dd, FUN = mean)
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = CellType), size = 0.2) + 
    geom_label(data = mids, aes(label = CellType)) + ggtitle("R100_63_Ctrl_F") + 
    guides(colour = guide_legend(override.aes = list(size = 1)))
plot(gg)


sc = autoEstCont(sc)

out = adjustCounts(sc)

cntSoggy = rowSums(sc$toc > 0)
cntStrained = rowSums(out > 0)
mostZeroed = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
mostZeroed

#If on the other hand we focus on genes for which there is a quantitative difference
tail(sort(rowSums(sc$toc > out)/rowSums(sc$toc > 0)), n = 20)

plotChangeMap(sc, out, "CLIC6")
plotChangeMap(sc, out, "CLDN5")
plotChangeMap(sc, out, "COL3A1")
plotChangeMap(sc, out, "SKAP1")
plotChangeMap(sc, out, "CNP")

DropletUtils:::write10xCounts("../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_063_new/outs/strainedCounts", out)

toc = Seurat::Read10X('../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_064_new/outs/filtered_feature_bc_matrix/')
tod = Seurat::Read10X('../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_064_new/outs/raw_feature_bc_matrix/')
sc = SoupChannel(tod, toc)

clt = read.csv('../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_064_new/outs/analysis/clustering/kmeans_6_clusters/clusters.csv',header = T)
tsne = read.csv('../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_064_new/outs/analysis/tsne/2_components/projection.csv',header = T)
rownames(tsne) = tsne$Barcode

sc = setClusters(sc, setNames(clt$Cluster, rownames(clt)))
sc = setDR(sc, tsne[colnames(sc$toc), c("TSNE.1", "TSNE.2")])

#clustering K-mean=6
library(ggplot2)
dd = sc$metaData[colnames(sc$toc), ]
mids = aggregate(cbind(TSNE.1, TSNE.2) ~ clusters, data = dd, FUN = mean)
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = clusters), size = 0.2) + 
    geom_label(data = mids, aes(label = clusters)) + ggtitle("R100_64_Ctrl_M") + 
    guides(colour = guide_legend(override.aes = list(size = 1)))
plot(gg)

#epithelial
dd$CLIC6 = sc$toc["CLIC6", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = CLIC6 > 0))
plot(gg)

#mesenchyal
dd$CLDN5 = sc$toc["CLDN5", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = CLDN5 > 0))
plot(gg)
#endothelial
dd$COL3A1 = sc$toc["COL3A1", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = COL3A1 > 0))
plot(gg)
#immune C1qa
dd$SKAP1 = sc$toc["SKAP1", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = SKAP1 > 0))
plot(gg)
#neuroglial Cnp
dd$CNP = sc$toc["CNP", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = CNP > 0))
plot(gg)

sc = autoEstCont(sc)

out = adjustCounts(sc)

plotChangeMap(sc, out, "CLIC6")
plotChangeMap(sc, out, "CLDN5")
plotChangeMap(sc, out, "COL3A1")
plotChangeMap(sc, out, "SKAP1")
plotChangeMap(sc, out, "CNP")

DropletUtils:::write10xCounts("../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_064_new/outs/strainedCounts", out)

toc = Seurat::Read10X('../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_072_new/outs/filtered_feature_bc_matrix/')
tod = Seurat::Read10X('../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_072_new/outs/raw_feature_bc_matrix/')
sc = SoupChannel(tod, toc)

clt = read.csv('../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_072_new/outs/analysis/clustering/kmeans_6_clusters/clusters.csv',header = T)
tsne = read.csv('../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_072_new/outs/analysis/tsne/2_components/projection.csv',header = T)
rownames(tsne) = tsne$Barcode

sc = setClusters(sc, setNames(clt$Cluster, rownames(clt)))
sc = setDR(sc, tsne[colnames(sc$toc), c("TSNE.1", "TSNE.2")])

#clustering K-mean=6
library(ggplot2)
dd = sc$metaData[colnames(sc$toc), ]
mids = aggregate(cbind(TSNE.1, TSNE.2) ~ clusters, data = dd, FUN = mean)
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = clusters), size = 0.2) + 
    geom_label(data = mids, aes(label = clusters)) + ggtitle("R100_72_Ctrl_F") + 
    guides(colour = guide_legend(override.aes = list(size = 1)))
plot(gg)

#epithelial
dd$CLIC6 = sc$toc["CLIC6", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = CLIC6 > 0))
plot(gg)

#mesenchyal
dd$PLVAP = sc$toc["PLVAP", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = PLVAP > 0))
plot(gg)
#endothelial
dd$COL3A1 = sc$toc["COL3A1", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = COL3A1 > 0))
plot(gg)
#immune C1qa
dd$SKAP1 = sc$toc["SKAP1", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = SKAP1 > 0))
plot(gg)
#neuroglial Cnp
dd$CNP = sc$toc["CNP", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = CNP > 0))
plot(gg)

sc = autoEstCont(sc)

out = adjustCounts(sc)

plotChangeMap(sc, out, "CLIC6")
plotChangeMap(sc, out, "PLVAP")
plotChangeMap(sc, out, "COL3A1")
plotChangeMap(sc, out, "SKAP1")
plotChangeMap(sc, out, "CNP")

DropletUtils:::write10xCounts("../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_072_new/outs/strainedCounts", out)

toc = Seurat::Read10X('../../../../raw_counts_scRNA/10X_20_AF_09_R77/10X_20_060_premrna/outs/filtered_feature_bc_matrix/')
tod = Seurat::Read10X('../../../../raw_counts_scRNA/10X_20_AF_09_R77/10X_20_060_premrna/outs/raw_feature_bc_matrix/')
sc = SoupChannel(tod, toc)

clt = read.csv('../../../../raw_counts_scRNA/10X_20_AF_09_R77/10X_20_060_premrna/outs/analysis/clustering/kmeans_6_clusters/clusters.csv',header = T)
tsne = read.csv('../../../../raw_counts_scRNA/10X_20_AF_09_R77/10X_20_060_premrna/outs/analysis/tsne/2_components/projection.csv',header = T)
rownames(tsne) = tsne$Barcode

sc = setClusters(sc, setNames(clt$Cluster, rownames(clt)))
sc = setDR(sc, tsne[colnames(sc$toc), c("TSNE.1", "TSNE.2")])

#clustering K-mean=6
library(ggplot2)
dd = sc$metaData[colnames(sc$toc), ]
mids = aggregate(cbind(TSNE.1, TSNE.2) ~ clusters, data = dd, FUN = mean)
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = clusters), size = 0.2) + 
    geom_label(data = mids, aes(label = clusters)) + ggtitle("R77_60_MS_M") + 
    guides(colour = guide_legend(override.aes = list(size = 1)))
plot(gg)

#epithelial
dd$CLIC6 = sc$toc["CLIC6", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = CLIC6 > 0))
plot(gg)

#mesenchyal
dd$PLVAP = sc$toc["PLVAP", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = PLVAP > 0))
plot(gg)
#endothelial
dd$COL3A1 = sc$toc["COL3A1", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = COL3A1 > 0))
plot(gg)
#immune C1qa
dd$SKAP1 = sc$toc["SKAP1", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = SKAP1 > 0))
plot(gg)
#neuroglial Cnp
dd$CNP = sc$toc["CNP", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = CNP > 0))
plot(gg)

sc = autoEstCont(sc)

out = adjustCounts(sc)

plotChangeMap(sc, out, "CLIC6")
plotChangeMap(sc, out, "PLVAP")
plotChangeMap(sc, out, "COL3A1")
plotChangeMap(sc, out, "SKAP1")
plotChangeMap(sc, out, "CNP")

DropletUtils:::write10xCounts("../../../../raw_counts_scRNA/10X_20_AF_09_R77/10X_20_060_premrna/outs/strainedCounts", out)

toc = Seurat::Read10X('../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_061_new/outs/filtered_feature_bc_matrix/')
tod = Seurat::Read10X('../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_061_new/outs/raw_feature_bc_matrix/')
sc = SoupChannel(tod, toc)

clt = read.csv('../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_061_new/outs/analysis/clustering/kmeans_6_clusters/clusters.csv',header = T)
tsne = read.csv('../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_061_new/outs/analysis/tsne/2_components/projection.csv',header = T)
rownames(tsne) = tsne$Barcode

sc = setClusters(sc, setNames(clt$Cluster, rownames(clt)))
sc = setDR(sc, tsne[colnames(sc$toc), c("TSNE.1", "TSNE.2")])

#clustering K-mean=6
library(ggplot2)
dd = sc$metaData[colnames(sc$toc), ]
mids = aggregate(cbind(TSNE.1, TSNE.2) ~ clusters, data = dd, FUN = mean)
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = clusters), size = 0.2) + 
    geom_label(data = mids, aes(label = clusters)) + ggtitle("R100_61_MS_F") + 
    guides(colour = guide_legend(override.aes = list(size = 1)))
plot(gg)

#epithelial
dd$CLIC6 = sc$toc["CLIC6", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = CLIC6 > 0))
plot(gg)

#mesenchyal
dd$PLVAP = sc$toc["PLVAP", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = PLVAP > 0))
plot(gg)
#endothelial
dd$COL3A1 = sc$toc["COL3A1", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = COL3A1 > 0))
plot(gg)
#immune C1qa
dd$SKAP1 = sc$toc["SKAP1", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = SKAP1 > 0))
plot(gg)
#neuroglial Cnp
dd$CNP = sc$toc["CNP", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = CNP > 0))
plot(gg)

sc = autoEstCont(sc)

out = adjustCounts(sc)

plotChangeMap(sc, out, "CLIC6")
plotChangeMap(sc, out, "PLVAP")
plotChangeMap(sc, out, "COL3A1")
plotChangeMap(sc, out, "SKAP1")
plotChangeMap(sc, out, "CNP")

DropletUtils:::write10xCounts("../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_061_new/outs/strainedCounts", out)

toc = Seurat::Read10X('../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_062_new/outs/filtered_feature_bc_matrix/')
tod = Seurat::Read10X('../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_062_new/outs/raw_feature_bc_matrix/')
sc = SoupChannel(tod, toc)

clt = read.csv('../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_062_new/outs/analysis/clustering/kmeans_6_clusters/clusters.csv',header = T)
tsne = read.csv('../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_062_new/outs/analysis/tsne/2_components/projection.csv',header = T)
rownames(tsne) = tsne$Barcode

sc = setClusters(sc, setNames(clt$Cluster, rownames(clt)))
sc = setDR(sc, tsne[colnames(sc$toc), c("TSNE.1", "TSNE.2")])

#clustering K-mean=6
library(ggplot2)
dd = sc$metaData[colnames(sc$toc), ]
mids = aggregate(cbind(TSNE.1, TSNE.2) ~ clusters, data = dd, FUN = mean)
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = clusters), size = 0.2) + 
    geom_label(data = mids, aes(label = clusters)) + ggtitle("R100_62_MS_F") + 
    guides(colour = guide_legend(override.aes = list(size = 1)))
plot(gg)

#epithelial
dd$CLIC6 = sc$toc["CLIC6", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = CLIC6 > 0))
plot(gg)

#mesenchyal
dd$PLVAP = sc$toc["PLVAP", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = PLVAP > 0))
plot(gg)
#endothelial
dd$COL3A1 = sc$toc["COL3A1", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = COL3A1 > 0))
plot(gg)
#immune C1qa
dd$SKAP1 = sc$toc["SKAP1", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = SKAP1 > 0))
plot(gg)
#neuroglial Cnp
dd$CNP = sc$toc["CNP", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = CNP > 0))
plot(gg)

sc = autoEstCont(sc)

out = adjustCounts(sc)

plotChangeMap(sc, out, "CLIC6")
plotChangeMap(sc, out, "PLVAP")
plotChangeMap(sc, out, "COL3A1")
plotChangeMap(sc, out, "SKAP1")
plotChangeMap(sc, out, "CNP")

DropletUtils:::write10xCounts("../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_062_new/outs/strainedCounts", out)

toc = Seurat::Read10X('../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_071_new/outs/filtered_feature_bc_matrix/')
tod = Seurat::Read10X('../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_071_new/outs/raw_feature_bc_matrix/')
sc = SoupChannel(tod, toc)

clt = read.csv('../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_071_new/outs/analysis/clustering/kmeans_6_clusters/clusters.csv',header = T)
tsne = read.csv('../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_071_new/outs/analysis/tsne/2_components/projection.csv',header = T)
rownames(tsne) = tsne$Barcode

sc = setClusters(sc, setNames(clt$Cluster, rownames(clt)))
sc = setDR(sc, tsne[colnames(sc$toc), c("TSNE.1", "TSNE.2")])

#clustering K-mean=6
library(ggplot2)
dd = sc$metaData[colnames(sc$toc), ]
mids = aggregate(cbind(TSNE.1, TSNE.2) ~ clusters, data = dd, FUN = mean)
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = clusters), size = 0.2) + 
    geom_label(data = mids, aes(label = clusters)) + ggtitle("R100_71_MS_F") + 
    guides(colour = guide_legend(override.aes = list(size = 1)))
plot(gg)

#epithelial
dd$CLIC6 = sc$toc["CLIC6", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = CLIC6 > 0))
plot(gg)

#mesenchyal
dd$PLVAP = sc$toc["PLVAP", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = PLVAP > 0))
plot(gg)
#endothelial
dd$COL3A1 = sc$toc["COL3A1", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = COL3A1 > 0))
plot(gg)
#immune C1qa
dd$SKAP1 = sc$toc["SKAP1", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = SKAP1 > 0))
plot(gg)
#neuroglial Cnp
dd$CNP = sc$toc["CNP", ]
gg = ggplot(dd, aes(TSNE.1, TSNE.2)) + geom_point(aes(colour = CNP > 0))
plot(gg)

sc = autoEstCont(sc)

out = adjustCounts(sc)

plotChangeMap(sc, out, "CLIC6")
plotChangeMap(sc, out, "PLVAP")
plotChangeMap(sc, out, "COL3A1")
plotChangeMap(sc, out, "SKAP1")
plotChangeMap(sc, out, "CNP")

DropletUtils:::write10xCounts("../../../../raw_counts_scRNA/10X_21_AF_16_R100/10X_21_071_new/outs/strainedCounts", out)










