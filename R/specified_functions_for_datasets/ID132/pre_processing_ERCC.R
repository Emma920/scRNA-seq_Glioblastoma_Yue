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

pre.processing.ERCC <- function(sce) {
  
  # QC MATRIX
  mito <- grepl("^MT-", rownames(sce))
  qc <- perCellQCMetrics(sce, subsets=list(Mito=mito))
  
  
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
  filtered <- sce[,!discard2]
  sce.filtered <- filtered[!mito, ]
  
  #FILTER OUT GENES
  sce.filtered <- sce.filtered[rowSums(counts(sce.filtered))> 100, ]
  
  
  #NORMALIZATION
  #lib.sf <- librarySizeFactors(sce.filtered)
  #sce.total <- logNormCounts(sce.filtered, size_factors = lib.sf, pseudo_count = max(1, 1/min(lib.sf)-1/max(lib.sf)))
  #cpm(sce.total) <- calculateCPM(sce.filtered)
  
  isSpike(sce.filtered, "ERCC") <- grepl("^ERCC-", rownames(sce.filtered))
  sce.total <- splitAltExps(sce.filtered, ifelse(isSpike(sce.filtered), "ERCC", "gene"))
  sce.total <- computeSpikeFactors(sce.total, "ERCC")
  to.plot <- data.frame(SpikeFactor=sizeFactors(sce.total))
  sce.total <- logNormCounts(sce.total, size_factors=to.plot$SpikeFactor)
  
  #summary(sizeFactors(sce.ID132))
  # to.plot <- data.frame(
  #   DeconvFactor=calculateSumFactors(sce.ID132),
  #   SpikeFactor=sizeFactors(sce.ID132)
  # )
  # ggplot(to.plot, aes(x=DeconvFactor, y=SpikeFactor)) +
  #   geom_point() + scale_x_log10() + 
  #   scale_y_log10() + geom_abline(intercept=0, slope=1, color="red")
  
  return(sce.total)
}

