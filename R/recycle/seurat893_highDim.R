#The update includes getting higher dimensions in PCA; Include all samples; Include meta data that has patients' ID and gene descriptions


# LIBRARIES
library(tidyverse)
library(plyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(Seurat)


# LOAD DATA
data_h5ad <- ReadH5AD('/home/hdd/alex/ncbi/public/sra/output/Trimmed/aligned/countMatrix_aggregate3.h5ad')
count_matrix <- GetAssayData(data_h5ad, slot = 'data') %>% as.matrix
mydata0 <- CreateSeuratObject(counts = count_matrix, meta.data = data_h5ad@meta.data)

# INSPECT DATA
View(rownames(mydata0))
View(head(mydata0@meta.data))


# FROM https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html

#mydata0[["percent.mt"]] <- PercentageFeatureSet(object = mydata0, pattern = "^MT-")
#VlnPlot(mydata0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#plot1 <- FeatureScatter(mydata0, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(mydata0, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#CombinePlots(plots = list(plot1, plot2))

#FILTER CELLS
#mydata1 <- subset(mydata0, subset = nFeature_RNA > 2000 & nFeature_RNA < 9000 & percent.mt < 18)

#NORMALIZATION
mydata1 <- NormalizeData(mydata0, normalization.method = "LogNormalize", scale.factor = 10000)

#FEATURE SELECTION
mydata1 <- FindVariableFeatures(mydata1, selection.method = "vst", nfeatures = 15000)
View(rownames(mydata1))


# Identify the 10 most highly variable genes
top10_variable <- head(VariableFeatures(mydata1), 10)
top20_variable <- head(VariableFeatures(mydata1), 20)
top100_variable <- head(VariableFeatures(mydata1), 100)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(mydata1)
LabelPoints(plot = plot1, points = top10_variable, repel = TRUE)

# SCALE GENE EXPRESSION
all.genes <- rownames(mydata1)
mydata1 <- ScaleData(mydata1, features = all.genes)

#PCA
mydata1 <- RunPCA(mydata1, features = VariableFeatures(object = mydata1))
print(mydata1[["pca"]], dims = 1:5, nfeatures = 20)
VizDimLoadings(mydata1, dims = 1:2, reduction = "pca")
DimPlot(mydata1, reduction = "pca")

DimHeatmap(mydata1, dims = 1:3, cells = 500, balanced = TRUE)

#Determine the ‘dimensionality’ of the dataset (https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html)
mydata1 <- JackStraw(mydata1, num.replicate = 100)
mydata1 <- ScoreJackStraw(mydata1, dims = 1:20) #dims 20 is maximum
JackStrawPlot(mydata1, dims = 1:20)
ElbowPlot(mydata1)

#cluster cells
mydata1 <- FindNeighbors(mydata1, dims = 1:50) #first 10 PCs
mydata1 <- FindClusters(mydata1, resolution = 0.5) # "increased values leading to a greater number of clusters. We find that setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells"
#head(Idents(pbmc), 5) #the cluster ID of the first 5 cells. The clusters can be found using the Idents function

#VISUALIZATION
mydata1 <- RunUMAP(mydata1, dims = 1:50)
DimPlot(mydata1, reduction = "umap", cols.highlight = "#SOX2")

#PLOT WITH FEATURES
features_stemlike <- c("CD44", "CD34", "SOX2", "KLF4", "NOTCH1", "POU2F1", "PROM1")
feature_differentiation <- c("S100B", "KCNJ10", "ID4", "METRN", "GFAP", "PLP1")
feature_proliferation <- c("TIMP1", "EGFR")
features_metabolic <- c("GAPDH", "BPGM", "ENO1", "ENO2", "PGK1", "GPI", "ALDOA", "PGAM1", "PGAM2")
feature_metastasis <- c("VIM", "TGFBI", "FAS", "RHOA")
feature_angiogenesis <- c("SERPINE1", "LOX", "EPAS1", "NRP1", "VEGFA", "ITGA5") #, "MAP2K1", "ANG", "CAV1", "ANXA2", "FN1", "PTGS2", "ITGB1", "ADM")
feature_HIF1 <- c("SERPINE1", "GAPDH", "EGLN3", "HK2", "MAP2K1", "EIF4EBP1")#, "PGK1", "ALDOA", "SLC2A1", "ENO2", "LDHA", "VEGFA")
feature_others <- c("UBB", "F3")
feature_cellcycle <- c("RHOA", "SPP1")

FeaturePlot(mydata1, features = features_stemlike)

#SAVE
saveRDS(mydata1, file = "/home/hdd/alex/ncbi/public/sra/output/Trimmed/aligned/mydata1_aggregated.rds")

############################################################################################
#FIND CLUSTER BIOMARKERS

cluster0.markers <- FindMarkers(mydata1, ident.1 = 1, min.pct = 0.25)
cluster0.top100 <- head(cluster1.markers, n = 100)
write.table(cluster0.top100, file="cluster0.top100")

cluster1.markers <- FindMarkers(mydata1, ident.1 = 1, min.pct = 0.25)
cluster1.top100 <- head(cluster1.markers, n = 100)
write.table(cluster1.top100, file="cluster1.top100")

cluster2.markers <- FindMarkers(mydata1, ident.1 = 2, min.pct = 0.25)
cluster2.top100 <- head(cluster2.markers, n = 100)
write.table(cluster2.top100, file="cluster2.top100")

cluster3.markers <- FindMarkers(mydata1, ident.1 = 3, min.pct = 0.25)
cluster3.top100 <- head(cluster3.markers, n = 100)
write.table(cluster3.top100, file="cluster3.top100")

cluster4.markers <- FindMarkers(mydata1, ident.1 = 4, min.pct = 0.25)
cluster4.top100 <- head(cluster4.markers, n = 100)
write.table(cluster4.top100, file="cluster4.top100")

cluster5.markers <- FindMarkers(mydata1, ident.1 = 5, min.pct = 0.25)
cluster5.top100 <- head(cluster5.markers, n = 100)
write.table(cluster5.top100, file="cluster5.top100")

cluster6.markers <- FindMarkers(mydata1, ident.1 = 6, min.pct = 0.25)
cluster6.top100 <- head(cluster6.markers, n = 100)
write.table(cluster6.top100, file="cluster6.top100")

cluster7.markers <- FindMarkers(mydata1, ident.1 = 7, min.pct = 0.25)
cluster7.top100 <- head(cluster7.markers, n = 100)
write.table(cluster7.top100, file="cluster7.top100")

cluster8.markers <- FindMarkers(mydata1, ident.1 = 8, min.pct = 0.25)
cluster8.top100 <- head(cluster8.markers, n = 100)
write.table(cluster8.top100, file="cluster8.top100")

cluster9.markers <- FindMarkers(mydata1, ident.1 = 9, min.pct = 0.25)
cluster9.top100 <- head(cluster9.markers, n = 100)
write.table(cluster9.top100, file="cluster9.top100")

cluster10.markers <- FindMarkers(mydata1, ident.1 = 10, min.pct = 0.25)
cluster10.top100 <- head(cluster10.markers, n = 100)
write.table(cluster10.top100, file="cluster10.top100")


cluster0.markers <- FindMarkers(mydata1, ident.1 = 1, min.pct = 0.25)
cluster0.top500 <- head(cluster1.markers, n = 500)
write.table(cluster0.top500, file="cluster0.top500")

cluster1.markers <- FindMarkers(mydata1, ident.1 = 1, min.pct = 0.25)
cluster1.top500 <- head(cluster1.markers, n = 500)
write.table(cluster1.top500, file="cluster1.top500")

cluster2.markers <- FindMarkers(mydata2, ident.1 = 1, min.pct = 0.25)
cluster2.top500 <- head(cluster2.markers, n = 500)
write.table(cluster2.top500, file="cluster2.top500")

cluster3.markers <- FindMarkers(mydata1, ident.1 = 1, min.pct = 0.25)
cluster3.top500 <- head(cluster3.markers, n = 500)
write.table(cluster3.top500, file="cluster3.top500")

cluster4.markers <- FindMarkers(mydata1, ident.1 = 1, min.pct = 0.25)
cluster4.top500 <- head(cluster4.markers, n = 500)
write.table(cluster4.top500, file="cluster4.top500")

cluster5.markers <- FindMarkers(mydata1, ident.1 = 1, min.pct = 0.25)
cluster5.top500 <- head(cluster5.markers, n = 500)
write.table(cluster5.top500, file="cluster5.top500")

cluster6.markers <- FindMarkers(mydata1, ident.1 = 1, min.pct = 0.25)
cluster6.top500 <- head(cluster6.markers, n = 500)
write.table(cluster6.top500, file="cluster6.top500")

cluster7.markers <- FindMarkers(mydata1, ident.1 = 1, min.pct = 0.25)
cluster7.top500 <- head(cluster7.markers, n = 500)
write.table(cluster7.top500, file="cluster7.top500")

cluster8.markers <- FindMarkers(mydata1, ident.1 = 1, min.pct = 0.25)
cluster8.top500 <- head(cluster8.markers, n = 500)
write.table(cluster8.top500, file="cluster8.top500")

cluster9.markers <- FindMarkers(mydata1, ident.1 = 1, min.pct = 0.25)
cluster9.top500 <- head(cluster9.markers, n = 500)
write.table(cluster9.top500, file="cluster9.top500")

cluster10.markers <- FindMarkers(mydata1, ident.1 = 1, min.pct = 0.25)
cluster10.top500 <- head(cluster10.markers, n = 500)
write.table(cluster10.top500, file="cluster10.top500")

# find markers for every cluster compared to all remaining cells, report only the positive ones
mydata1.markers <- FindAllMarkers(mydata1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
view5Markers <- mydata1.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
VlnPlot(mydata1, features = c("COL1A2", "RPL21", "SULF1", "PLP1", "TSPAN31", "KCNE5", "SLC2A1"))
# heat map of the first 10 marker genes
top5 <- mydata1.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(mydata1, features = top5$gene, size = 5) + NoLegend()
