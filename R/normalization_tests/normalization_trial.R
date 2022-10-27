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


  
  # QC MATRIX
  mito <- grepl("^MT-", rownames(sce.BT_S2))
  qc <- perCellQCMetrics(sce.BT_S2, subsets=list(Mito=mito))
  
  
  # FIXED CRITERIA
  criteria.lib <- qc$sum < 5e4
  criteria.nexprs <- qc$detected < 3e3
  criteria.mito <- qc$subsets_Mito_percent > 10
  discard <- criteria.lib | criteria.nexprs | criteria.mito
  
  
  # ADAPTIVE CRITERIA
  criteria.lib2 <- isOutlier(qc$sum, log=TRUE, type="lower")
  criteria.nexprs2 <- isOutlier(qc$detected, log=TRUE, type="lower")
  criteria.mito2 <- isOutlier(qc$subsets_Mito_percent, type="higher")
  discard2 <- criteria.lib2 | criteria.nexprs2 | criteria.mito2
  # Summarize the number of cells removed for each reason. https://osca.bioconductor.org/quality-control.html#identifying-low-quality-cells
  #DataFrame(LibSize=sum(criteria.lib2), NExprs=sum(criteria.nexprs2), MitoProp=sum(criteria.mito2), Total=sum(discard2))
  
  
  #FILTER OUT CELLS
  filtered <- sce.BT_S2[,!discard2]
  sce.BT_S2_filtered <- filtered[!mito, ]
  
  #FILTER OUT GENES
  sce.BT_S2_filtered <- sce.BT_S2_filtered[rowSums(counts(sce.filtered))> 100, ]
  
  
  #NORMALIZATION
  set.seed(100)
  clust <- quickCluster(sce.BT_S2_filtered) # This is needed to calculate size factors
  #table(clust)
  #summary(clust) # To see it
  deconv.sf <- calculateSumFactors(sce.BT_S2_filtered, cluster=clust)
  sum.factor<- computeSumFactors(sce.filtered, cluster=clust, min.mean=0.1) #Normalize the countMatrix
  lib.sf <- librarySizeFactors(sce.filtered)
  sce.total <- logNormCounts(sce.total)
  cpm(sce.total) <- calculateCPM(sce.filtered)
  
