library(tidyverse)
library(plyr)
library(dplyr)
library(data.table)

setwd("yue/data/aligned/20200313")
dataList <- list.files(pattern= "*ReadsPerGene.out.tab")
rawMatrix <- as.data.frame(read.table(dataList[1],header=TRUE)[-1:-4, -2:-3])
colnames(rawMatrix) <- c("gene_id", str_sub(dataList[1], 1, 2))
for(i in 1:(length(dataList)-1)){
  rawMatrix1 <- as.data.frame(read.table(dataList[i+1],header=TRUE)[-1:-4, -2:-3])
  colnames(rawMatrix1) <- c("gene_id", str_sub(dataList[i+1], 1, 2))
  rawMatrix<- join(rawMatrix, rawMatrix1, by=c("gene_id"))
  print(i)
}