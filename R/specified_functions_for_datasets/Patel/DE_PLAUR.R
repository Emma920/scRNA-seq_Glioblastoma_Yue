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


DE.PLAUR <- function(sceP.object) {
  logcountsP_PLAUR <- logcounts(sceP.object["PLAUR",])
  factor_test2 <- factor(logcounts(sceP.object["PLAUR",]) > max(logcounts_PLAUR["PLAUR",])-0.75*(max(logcounts_PLAUR["PLAUR",])-min(logcounts_PLAUR["PLAUR",])))
  
  library(DEsingle)
  DE_result <- DEsingle(sce.object, factor_test2, parallel = FALSE, BPPARAM = bpparam())
  DE_rankByP <- DE_result[order(DE_result$pvalue),] 
  return(DE_rankByP)
}