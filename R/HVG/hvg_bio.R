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

hvg.bio <- function(sce.total) {

  #HIGHLY VARIABLE GENES CONSIDERING BATCH EFFECT
  dec.bio.total <- modelGeneVarByPoisson(sce.total) #dec.total considering plate_id
  hvg.bio <- getTopHVGs(dec.bio.total, var.threshold=0) # hvg considering plate_id
  sce.bio.hvg <- sce.total[hvg.bio,] #sce.hvg considering plate_id
  #head(dec.bio.total[order(dec.bio.total$bio, decreasing=TRUE),1:6]) # check dec.bio.total
  
  # PLOT TREND OF DEC.BLOCK.TOTAL
  # par(mfrow=c(1,2)) # check dec.block.total
  # blocked.stats <- dec.block.total$per.block # check dec.block.total
  # blocked.stats@listData[["X1001000028"]] # check dec.block.total by each batch
  # for (i in colnames(blocked.stats)) {
  #   current <- blocked.stats[[i]]
  #   plot(current$mean, current$total, main=i, pch=16, cex=0.5, xlab="Mean of log-expression", ylab="Variance of log-expression")
  #   xlim <- c(0,20)
  #   ylim <- c(0,40)
  #   plot.window(xlim,ylim )
  #   curfit <- metadata(current)
  #   points(curfit$mean, curfit$var, col="red", pch=16)
  #   curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2) 
  # }
  
  # plot batch trend one by one
  #current <- blocked.stats[["X1001000174"]]
  #xlim <- c(0,20)
  #ylim <- c(0,40)
  #plot.window(xlim,ylim )
  #plot(current$mean, current$total, main="X1001000174", pch=16, cex=0.5, xlab="Mean of log-expression", ylab="Variance of log-expression")
  #curfit <- metadata(current)
  #points(curfit$mean, curfit$var, col="red", pch=16)
  #curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2) 
  
  return(sce.bio.hvg)
}