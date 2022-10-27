library(tidyverse)
library(plyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(msigdbr)
library(SingleCellExperiment)
library(scater)
library(scran)


edgeR <- function(sce.object) {
  factor_test <- factor(sce.object@colData@listData[["neoplastic"]]=="Neoplastic")
  library(edgeR)
  dge <- DGEList(
    counts = logcounts(sce.object),
    lib.size = NULL,
    norm.factors = NULL, 
    group = factor_test
  )
  group_edgeR <- factor_test
  design <- model.matrix(~ group_edgeR)
  dge <- estimateDisp(dge)
  fit <- glmFit(dge, design) 
  res <- glmLRT(fit)
  pVals <- res$table[,4]
  #result <- res$table
  result2 <- topTags(res, n = 10000, adjust.method = "BH", sort.by = "PValue", p.value = 1)$table
  #names(pVals) <- rownames(res$table)
  #pVals_rank <- order(pVals) 
  #pVals_ranked <- as.matrix(pVals[pVals_rank])
  return(result2)
}

##########################################################################################




