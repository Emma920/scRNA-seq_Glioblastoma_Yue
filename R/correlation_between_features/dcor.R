library("energy")
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

distance.correlation <- function(gene.set1, gene.set2, object) {
  set1.counts <- logcounts(object[gene.set1,])
  set2.counts <- logcounts(object[gene.set2,])
  
  dcor.info <- dcor(t(set1.counts), t(set2.counts))
  return(dcor.info)
}

distance.correlation.ttest <- function(gene.set1, gene.set2, object) {
  set1.counts <- logcounts(object[gene.set1,])
  set2.counts <- logcounts(object[gene.set2,])
  
  dcor.info.ttest <- dcor.ttest(t(set1.counts), t(set2.counts), distance=FALSE)
  return(dcor.info.ttest)
}


#gene.set1 and gene.set2 are gene sets that are interested in looking for correlations in, and object is the 
#overall object where the gene sets are taken from, normally is sce.hvg.
