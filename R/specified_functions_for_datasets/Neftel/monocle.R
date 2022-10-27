library(tidyverse)
library(plyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(msigdbr)
library(SingleCellExperiment)
library(scater)
library(scran)


monocle <- function(sce.object, gene) {
  logcounts_gene <- logcounts(sce.object[gene,])
  factor_test <- factor(logcounts(sce.object[gene,]) > max(logcounts_gene[gene,])-0.75*(max(logcounts_gene[gene,])-min(logcounts_gene[gene,])))
  library(monocle)
  pd <- data.frame(group = factor_test)
  rownames(pd) <- colnames(logcounts(sce.object))
  pd <- new("AnnotatedDataFrame", data = pd)
  
  Obj <- newCellDataSet(
    as.matrix(logcounts(sce.object)), 
    phenoData = pd, 
    expressionFamily = negbinomial.size()
  )
  Obj <- estimateSizeFactors(Obj)
  Obj <- estimateDispersions(Obj)
  res <- differentialGeneTest(Obj, fullModelFormulaStr = "~group")
  return(res)
}