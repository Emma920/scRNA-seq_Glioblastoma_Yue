
library(tidyverse)
library(plyr)
library(dplyr)
library(data.table)
library(SingleCellExperiment)

load.data.Tirosh <- function(cm.original, meta.data) {
  countMatrix_aggregate <- read.table(cm.original, header = TRUE) #countMatrix_aggregate was made by collecting est_counts from each cell and transcripts with the same gene names were aggregated. None of the cells or genes were omited
  countMatrix_aggregate1 <- countMatrix_aggregate[!duplicated(countMatrix_aggregate[1]), ]
  data1 <- data.frame(countMatrix_aggregate1[,-1], row.names=countMatrix_aggregate1[,1])
  library(readr)
  meta_data1 <- read_delim(meta.data, 
                            "\t", escape_double = FALSE, col_types = cols(NAME = col_factor(levels = c()), 
                                                                          Sample = col_factor(levels = c()), 
                                                                          Type = col_factor(levels = c())), 
                                                                          trim_ws = TRUE)
  sce <- SingleCellExperiment(assays = list(tpm = as.matrix(data1)), colData = meta_data1)
  return(sce)
}







