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


kNN.dcor <- function(gene.set1, gene.set2, sce.object) {
  k <- findKNN(t(logcounts(sce.object)), k=30)
  dcor_cells_obj <- vector()
  for (i in 1:(ncol(sce.object))){
    index <- k[["index"]][i, ]
    sce_cell <- sce.object[, index]
    dcor_cells <- distance.correlation(gene.set1, gene.set2, sce_cell)
    dcor_cells_obj <- c(dcor_cells_obj, dcor_cells) 
  }
  names(dcor_cells_obj) <- colnames(sce.object)
  return(dcor_cells_obj)
}