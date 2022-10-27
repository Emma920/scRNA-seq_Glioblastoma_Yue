library(tidyverse)
library(plyr)
library(dplyr)
library(data.table)
library(SingleCellExperiment)
library(scater)
library(scran)
library(Rtsne)
library(PCAtools)

clusting.tpm <- function(sce) {
  g <- buildSNNGraph(sce, assay.type = "tpm", k=10, use.dimred = NULL)
  clust <- igraph::cluster_louvain(g)$membership
  return(clust)
}