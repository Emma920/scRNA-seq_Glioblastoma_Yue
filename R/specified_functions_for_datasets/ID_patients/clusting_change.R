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
  m <- read.csv("/home/yue/hdd/yue/data/text/ID_patients/clustering_parameters_0.csv", row.names = 1)
  g <- buildSNNGraph(sce, k=m[str_sub(deparse(substitute(sce)), 5, 9), "k"], use.dimred = NULL)
  clust <- igraph::cluster_leiden(g, resolution_parameter = m[str_sub(deparse(substitute(sce)), 5, 9), "resolution"])$membership
  #g <- buildSNNGraph(sce, k=11, use.dimred = NULL)
  #clust <- igraph::cluster_leiden(g, resolution_parameter = 0.8)$membership
  
  
  return(clust)
}

