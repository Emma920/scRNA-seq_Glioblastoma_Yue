library(tidyverse)
library(plyr)
library(dplyr)
library(data.table)
library(SingleCellExperiment)
library(scater)
library(scran)
library(Rtsne)
library(PCAtools)

clusting.change <- function(sce) {
  m <- read.csv("/home/yue/hdd/yue/data/text/prom1/clustering_parameters_2nd.csv", row.names = 1)
  g <- buildSNNGraph(sce, k=m[deparse(substitute(sce)), "k"], use.dimred = NULL)
  clust <- igraph::cluster_leiden(g, resolution_parameter = m[deparse(substitute(sce)), "resolution"])$membership

  return(clust)
}

# g <- buildSNNGraph(sce.BT_S2.total, k=10, use.dimred = NULL)
# clust <- igraph::cluster_louvain(g)$membership
# sce.BT_S2.total@colData@listData[["label"]] <- clust
# markers.BT_S2 <- findMarkers(sce.BT_S2.total, sce.BT_S2.total@colData@listData[["label"]])