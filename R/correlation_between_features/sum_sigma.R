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
sum.sigma <- function(gene.set, sce.object) {
  sum.feature <- colSums(logcounts(sce.object[gene.set, ]))

  rand.sum <- vector()
  set.seed(100)
  for (i in 1:100) {
    rand <- sample(nrow(sce.object), length(gene.set))
    sce.rand <- sce.object[rand,]
    sum <- colSums(logcounts(sce.rand))
    rand.sum <- c(rand.sum, sum)
  }

  sigma <- (sum.feature-mean(rand.sum))/sd(rand.sum)

  return(sigma)

}