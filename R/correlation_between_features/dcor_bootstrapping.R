library(tidyverse)
library(plyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(msigdbr)
library(SingleCellExperiment)
library(scater)
library(scran)


# dcor.bootstrapping <- function(sets1, sets2, genes, sce.object) {
#   source("/home/hdd/yue/code/R/Darmanis/dcor.R")
#   dcor.bt <- data.frame(dcor = numeric(0))
#   sce.minus.sets1_sets2 <- sce.object[!(genes %in% c(sets1, sets2)), ]
#   
#   set.seed(100)
#   for (i in 1:500) {
#     rand <- sample(nrow(sce.minus.sets1_sets2), length(sets2))
#     dis <- distance.correlation(sets1, rand, sce.object)
#     dcor.bt <- rbind.data.frame(dcor.bt, dis) 
#   }
#   colnames(dcor.bt) <- c("dcor")
#   return(dcor.bt)
# }

dcor.bootstrapping <- function(sets1, sets2, genes, sce.object) {
  source("/home/hdd/yue/code/R/Darmanis/dcor.R")
  dcor.bt <- vector()
  sce.minus.sets1_sets2 <- sce.object[!(genes %in% c(sets1, sets2)), ]
  
  set.seed(100)
  for (i in 1:100) {
    rand <- sample(nrow(sce.minus.sets1_sets2), length(sets2))
    dis <- distance.correlation(sets1, rand, sce.object)
    dcor.bt <- c(dcor.bt, dis) 
  }
  return(dcor.bt)
}