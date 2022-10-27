library(SingleCellExperiment)
library(scater)
library(scran)
library(uwot)
library(Rtsne)
library(tidyverse)
library(plyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(msigdbr)
library(PCAtools)
library("RColorBrewer")

setwd('/home/hdd/alex/ncbi/public/sra/output/Trimmed/aligned')

countMatrix_aggregate <- read.table("countMatrix_aggregate", header = TRUE) #countMatrix_aggregate was made by collecting est_counts from each cell and transcripts with the same gene names were aggregated. None of the cells or genes were omited
data1 <- data.frame(countMatrix_aggregate[,-1], row.names=countMatrix_aggregate[,1])
meta_data1 <- read.csv("neo_metadata.csv")
sce <- SingleCellExperiment(
  assays = list(counts = as.matrix(data1)), 
  colData = meta_data1
)

#QC matrix
mito <- grepl("^MT-", rownames(sce))
qc <- perCellQCMetrics(sce, subsets=list(Mito=mito))

#fixed criteria
criteria.lib <- qc$sum < 5e4
criteria.nexprs <- qc$detected < 3e3
criteria.mito <- qc$subsets_Mito_percent > 10
discard <- criteria.lib | criteria.nexprs | criteria.mito

#adaptive criteria
criteria.lib2 <- isOutlier(qc$sum, log=TRUE, type="lower")
criteria.nexprs2 <- isOutlier(qc$detected, log=TRUE, type="lower")
criteria.mito2 <- isOutlier(qc$subsets_Mito_percent, type="higher")
discard2 <- criteria.lib2 | criteria.nexprs2 | criteria.mito2

# Summarize the number of cells removed for each reason. https://osca.bioconductor.org/quality-control.html#identifying-low-quality-cells
DataFrame(LibSize=sum(criteria.lib2), NExprs=sum(criteria.nexprs2), MitoProp=sum(criteria.mito2), Total=sum(discard2))

#Filter out cells
filtered <- sce[,!discard2]
sce.filtered <- filtered[!mito, ]

#NORMALIZATION
set.seed(100)
clust.zeisel <- quickCluster(sce.filtered) # This is needed to calculate size factors
#table(clust.zeisel)

#deconv.sf.zeisel <- calculateSumFactors(sce.final, cluster=clust.zeisel) # This is just to calculate a value called deconv.sf.zeisel, didn't apply on matrix
#summary(deconv.sf.zeisel) # To see it

sce.zeisel <- computeSumFactors(sce.filtered, cluster=clust.zeisel, min.mean=0.1) #Normalize the countMatrix
sce.zeisel <- logNormCounts(sce.zeisel)
assayNames(sce.zeisel)

#FIND HVG (feature selection)
###################################find the trend (technical noise, not considering uninterested )
dec.zeisel <- modelGeneVar(sce.zeisel)
fit.zeisel <- metadata(dec.zeisel)
plot(fit.zeisel$mean, fit.zeisel$var, xlab="Mean", ylab="Variance")
curve(fit.zeisel$trend(x), col="dodgerblue", add=TRUE, lwd=2) # visualizing the fit
dec.zeisel[order(dec.zeisal$bio, decreasing=TRUE),] # Ordering by most interesting genes for inspection.
###################################find the trend (technical noise+biological noise)
#set.seed(0010101)
#dec.pois.zeisel <- modelGeneVarByPoisson(sce.zeisel)
#dec.pois.zeisel <- dec.pois.zeisel[order(dec.pois.zeisel$bio, decreasing=TRUE),]
#head(dec.pois.zeisel)
#plot(dec.pois.zeisel$mean, dec.pois.zeisel$total, pch=16, xlab="Mean of log-expression",
#     ylab="Variance of log-expression")
#curve(metadata(dec.pois.zeisel)$trend(x), col="dodgerblue", add=TRUE)
######################################batch effect
######################################select HVG
hvg <- getTopHVGs(dec.zeisel, var.threshold=0)
sce.hvg <- sce.zeisel[hvg,]

#DIMENSION REDUCION
set.seed(100) # See below.
sce.zeisel <- runPCA(sce.zeisel, subset_row=hvg) 
reducedDimNames(sce.zeisel)
ncol(reducedDim(sce.zeisel))

########################################choosing amount of PCs
percent.var <- attr(reducedDim(sce.zeisel), "percentVar")
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
chosen.elbow
par(mar = c(5,5,3,1))
plot(percent.var, xlab="PC", ylab="Variance explained (%)")
abline(v=chosen.elbow, col="red")
#######################################Choose top 25 PCs
reducedDim(sce.zeisel, "PCA_25") <- reducedDim(sce.zeisel, "PCA")[,1:25]
reducedDimNames(sce.zeisel)
#######################################Visualizing reduced dimension

plotReducedDim(sce.zeisel, dimred="PCA", colour_by = "patient_id")

set.seed(00101001101)
sce.zeisel <- runTSNE(sce.zeisel, dimred="PCA")
plotReducedDim(sce.zeisel, dimred="TSNE", colour_by = "SOX2")

sce.zeisel <- runUMAP(sce.zeisel, dimred="PCA")
plotReducedDim(sce.zeisel, dimred="UMAP", colour_by="patient_id")

#CLUSTERING

g <- buildSNNGraph(sce.zeisel, k=8, use.dimred = 'PCA') #changing the number of k will change the number of clusters, like the "resolution" in FindClusters (seurat)
clust <- igraph::cluster_walktrap(g)$membership
table(clust)
sce.zeisel$cluster <- factor(clust)
plotReducedDim(sce.zeisel, "TSNE", colour_by="cluster")
###################################################assessing clustering
library(pheatmap)
ratio <- clusterModularity(g, clust, as.ratio=TRUE)
dim(ratio)
pheatmap(log2(ratio+1), cluster_rows=FALSE, cluster_cols=FALSE,
         color=colorRampPalette(c("white", "blue"))(100))


#FIND MARKER GENES FOR CLUSTERS
markers.zeisel <- findMarkers(sce.zeisel, sce.zeisel$cluster, pval.type="some")#The summary.logFC field provides a convenient summary of the direction and effect size for each gene, and is defined here as the log-fold change from the comparison with the lowest p-value
markers.zeisel

markers_cluster9 <- markers.zeisel[["1"]]
colnames(markers_cluster9)
markers_cluster9[1:20,1:3]


markers <- combineMarkers(markers.zeisel@listData , pairs, pval.field = "p.value",effect.field = "logFC",pval.type = c("any"),)

cluster9_top6 <- markers_cluster9[markers_cluster9$Top == 1,]
logFCs_cluster9 <- getMarkerEffects(cluster9_top6)
library(pheatmap)
pheatmap(logFCs_cluster9, breaks=seq(-5, 5, length.out=101))


#TRY DB
c2.sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP")
head(unique(c2.sets$gs_name))
csc.sets <- c2.sets[grep("BEIER", c2.sets$gs_name),]#this gave the desired gene sets

