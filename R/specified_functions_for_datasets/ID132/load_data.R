library(tidyverse)
library(plyr)
library(dplyr)
library(data.table)
library(SingleCellExperiment)

load.data.GT <- function(cm.original, meta.data) {
  countMatrix_aggregate <- read.table(cm.original, header = TRUE) #countMatrix_aggregate was made by collecting est_counts from each cell and transcripts with the same gene names were aggregated. None of the cells or genes were omited
  countMatrix_aggregate1 <- countMatrix_aggregate[!duplicated(countMatrix_aggregate[1]), ]
  data1 <- data.frame(countMatrix_aggregate1[,-1], row.names=countMatrix_aggregate1[,1])
  meta_data1 <- read.csv(meta.data)
  sce <- SingleCellExperiment(assays = list(counts = as.matrix(data1)), colData = meta_data1)
  return(sce)
}