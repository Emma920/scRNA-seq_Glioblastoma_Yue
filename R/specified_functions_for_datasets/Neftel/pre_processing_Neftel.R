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

pre.processing.Neftel <- function(sce) {
  
  # QC MATRIX
  mito <- grepl("^MT-", rownames(sce))
  qc <- perCellQCMetrics(sce, subsets=list(Mito=mito))

  
  # ADAPTIVE CRITERIA
  criteria.lib2 <- isOutlier(qc$sum, log=TRUE, type="lower")
  criteria.nexprs2 <- isOutlier(qc$detected, log=TRUE, type="lower")
  criteria.mito2 <- isOutlier(qc$subsets_Mito_percent, type="higher")
  discard2 <- criteria.lib2 | criteria.nexprs2 | criteria.mito2
  # Summarize the number of cells removed for each reason. https://osca.bioconductor.org/quality-control.html#identifying-low-quality-cells
  #DataFrame(LibSize=sum(criteria.lib2), NExprs=sum(criteria.nexprs2), MitoProp=sum(criteria.mito2), Total=sum(discard2))
  
  
  #FILTER OUT CELLS
  filtered <- sce[,!discard2]
  sce.filtered <- filtered[!mito, ]
  
  #FILTER OUT GENES
  #discard_gene <- isOutlier(rowSums(tpm(sce.filtered)), type="lower")
  sce.total <- sce.filtered#[!discard_gene,]
  
  #NORMALIZATION
  # lib.sf <- librarySizeFactors(sce.filtered)
  # #sce.total <- logNormCounts(sce.filtered, size_factors = lib.sf, pseudo_count = max(1, 1/min(lib.sf)-1/max(lib.sf)))
  # sce.total <- logNormCounts(sce.filtered, size_factors = lib.sf, pseudo_count = 1)
  # cpm(sce.total) <- calculateCPM(sce.filtered)
  
  return(sce.total)
}


