library(tidyverse)
library(plyr)
library(dplyr)
library(data.table)
library(SingleCellExperiment)
library(scater)
library(scran)
library(Rtsne)
library(PCAtools)

clusting <- function(sce) {
  g <- buildSNNGraph(sce, k=10, use.dimred = NULL)
  clust <- igraph::cluster_louvain(g)$membership
  return(clust)
}

# g <- buildSNNGraph(sce.BT_S2.total, k=10, use.dimred = NULL)
# clust <- igraph::cluster_louvain(g)$membership
# sce.BT_S2.total@colData@listData[["label"]] <- clust
# markers.BT_S2 <- findMarkers(sce.BT_S2.total, sce.BT_S2.total@colData@listData[["label"]])