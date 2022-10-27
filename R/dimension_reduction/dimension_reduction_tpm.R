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


dimension.reduction.tpm <- function(selected.set, sce.objects) {
  sce.selected <- sce.objects[selected.set,]
  set.seed(100) # See below.
  sce.selected <- runPCA(sce.selected, subset_row=selected.set,exprs_values = "tpm")
  reducedDim(sce.selected, "PCA_25") <- reducedDim(sce.selected, "PCA")[,1:25]
  set.seed(00101001101)
  sce.selected <- runTSNE(sce.selected, dimred="PCA",exprs_values = "logcounts")
  sce.selected <- runUMAP(sce.selected, dimred="PCA",exprs_values = "logcounts")
  
  return(sce.selected)
  #
}