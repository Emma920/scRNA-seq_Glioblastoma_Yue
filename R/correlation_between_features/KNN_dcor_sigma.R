library(tidyverse)
library(plyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(msigdbr)
library(SingleCellExperiment)
library(scater)
library(scran)
library(BiocNeighbors)


source("/home/hdd/yue/code/R/Darmanis/dcor_bootstrapping.R")
kNN.dcor.sigma <- function(gene.set1, gene.set2, genes, sce.object) {
  k <- findKNN(t(logcounts(sce.object)), k=30)
  sigma_cells_obj <- vector()
  for (i in 1:(ncol(sce.object))){
    index <- k[["index"]][i, ]
    sce_cell <- sce.object[, index]
    dcor_cells <- distance.correlation(gene.set1, gene.set2, sce_cell)
    KNN_set1_rand.dcor <- dcor.bootstrapping(gene.set1, gene.set2, genes, sce_cell)
    sigma <- (dcor_cells-mean(KNN_set1_rand.dcor))/sd(KNN_set1_rand.dcor)
    sigma_cells_obj <- c(sigma_cells_obj, sigma) 
    print(i)
  }
  names(sigma_cells_obj) <- colnames(sce.object)
  return(sigma_cells_obj)
}