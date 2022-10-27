library(tidyverse)
library(plyr)
library(dplyr)
library(data.table)
library(SingleCellExperiment)
library(scater)
library(scran)
library(Rtsne)
library(PCAtools)

clusting.tpm.change <- function(sce) {
  m <- read.csv("/home/yue/hdd/yue/data/text/prom1/clustering_parameters_2nd.csv", row.names = 1)
  g <- buildSNNGraph(sce, assay.type = "tpm", k=m[deparse(substitute(sce)), "k"], use.dimred = NULL)
  clust <- igraph::cluster_leiden(g, resolution_parameter = m[deparse(substitute(sce)), "resolution"])$membership
  #g <- buildSNNGraph(sce, assay.type = "tpm", k=10, use.dimred = NULL)
  #clust <- igraph::cluster_leiden(g, resolution_parameter = 1)$membership


  return(clust)
}