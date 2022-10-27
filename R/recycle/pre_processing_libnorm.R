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

pre.processing.libnorm <- function(sce) {
  
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
  
  
  #NORMALIZATION
  set.seed(100)
  clust <- quickCluster(sce.filtered) # This is needed to calculate size factors
  #table(clust.zeisel)
  #deconv.sf.zeisel <- calculateSumFactors(sce.final, cluster=clust.zeisel) # This is just to calculate a value called deconv.sf.zeisel, didn't apply on matrix
  #summary(deconv.sf.zeisel) # To see it
  sce.total <- computeSumFactors(sce.filtered, cluster=clust, min.mean=0.1) #Normalize the countMatrix
  sce.total <- logNormCounts(sce.total)
  cpm(sce.total) <- calculateCPM(sce.filtered)
  
  #SELECT HIGHLY VARIABLE GENES
  dec.total <- modelGeneVar(sce.total) #original method
  head(dec.total[order(dec.total$bio, decreasing=TRUE),1:6]) # check original method
  #fit.zeisel <- metadata(dec.total)
  #plot(fit.zeisel$mean, fit.zeisel$var, xlab="Mean", ylab="Variance")
  #curve(fit.zeisel$trend(x), col="dodgerblue", add=TRUE, lwd=2) # visualizing the fit
  #dec.total[order(dec.total$bio, decreasing=TRUE),]
  hvg <- getTopHVGs(dec.total, var.threshold=0)  #original method
  sce.hvg <- sce.total[hvg,] #original method
  
  #par(mfrow=c(1,2)) # check dec.total
  #blocked.original.stats <- dec.total$per.block #check dec.total
  return(sce.hvg)
}
