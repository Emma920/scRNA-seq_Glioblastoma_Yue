library(tidyverse)
library(plyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(msigdbr)
library(SingleCellExperiment)
library(scater)
library(scran)
library(uwot)
library(Rtsne)
library(PCAtools)
library("RColorBrewer")
library(d3heatmap)


cluster.marker.gene <- function(sce.cscM, igraph, cluster_louvain) {
  g <- buildSNNGraph(sce.cscM, k=10, use.dimred = 'PCA')
  clust.louvain <- igraph::cluster_louvain(g)$membership
  sce.cscM$label <- factor(clust.louvain)
  
  markers.cscM <- findMarkers(sce.cscM, groups=sce.cscM$label)
}