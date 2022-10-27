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

hvg.standard <- function(sce.total) {
  
  #SELECT HIGHLY VARIABLE GENES
  dec.total <- modelGeneVar(sce.total)
  #fit.zeisel <- metadata(dec.total)
  #plot(fit.zeisel$mean, fit.zeisel$var, xlab="Mean", ylab="Variance")
  #curve(fit.zeisel$trend(x), col="dodgerblue", add=TRUE, lwd=2) # visualizing the fit
  #dec.total[order(dec.total$bio, decreasing=TRUE),]
  hvg <- getTopHVGs(dec.total, var.threshold=0)
  sce.hvg <- sce.total[hvg,]
  return(sce.hvg)
}
