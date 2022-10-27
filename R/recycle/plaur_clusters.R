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


dataList <- mget(ls(pattern = '*.total'))
plaur.clusters <- list()
for (i in 1:(length(dataList))){
  library(rlist)
  df <- data.frame(clusters =c(), 
                   cells_expressing=c(),
                   cells_total=c(),
                   percentages=c())
  for (j in 1: length(sort(unique(dataList[[i]]@colData@listData[["label"]])))){
    df0 <- data.frame(clusters =c(sort(unique(dataList[i][[names(dataList[i])]]@colData@listData[["label"]]))[j]), 
                      cells_expressing=length(which(logcounts(dataList[i][[names(dataList[i])]]["PLAUR",dataList[i][[names(dataList[i])]]@colData@listData[["label"]]==j])!=0)),
                      cells_total=ncol(logcounts(dataList[i][[names(dataList[i])]]["PLAUR",dataList[i][[names(dataList[i])]]@colData@listData[["label"]]==j])),
                      percentages=length(which(logcounts(dataList[i][[names(dataList[i])]]["PLAUR",dataList[i][[names(dataList[i])]]@colData@listData[["label"]]==j])!=0))*100/ncol(logcounts(dataList[i][[names(dataList[i])]]["PLAUR",dataList[i][[names(dataList[i])]]@colData@listData[["label"]]==j])))
    df <- rbind(df,df0)
    #assign(str_sub(names(dataList[i]), -11, -7), df)
  }
  plaur.clusters <- list.append(plaur.clusters, df)
  names(plaur.clusters)[i]<- str_sub(names(dataList[i]), 5, -7)
  
}



