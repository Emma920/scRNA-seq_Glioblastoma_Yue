library(tidyverse)
library(plyr)
library(dplyr)
library(data.table)
# 
# 
# knit_aggregate_GT <- function(dir, id) {
#   setwd(dir)
#   dataList <- list.files(pattern= id)
#   rawMatrix <- as.data.frame(read.table(dataList[1],header=TRUE)[-1:-4, -2:-3])
#   colnames(rawMatrix) <- c("gene_id", str_sub(dataList[1], 1, 2))
#   for(i in 1:(length(dataList)-1)){
#     rawMatrix1 <- as.data.frame(read.table(dataList[i+1],header=TRUE)[-1:-4, -2:-3])
#     colnames(rawMatrix1) <- c("gene_id", str_sub(dataList[i+1], 1, 2))
#     rawMatrix<- join(rawMatrix, rawMatrix1, by=c("gene_id"))
#     print(i)
#   }
#   
#   
assign_names <- function(dir) {
  geneNames <- read.csv("data/text/gene_name_ERCC.csv") 
  rawMatrix <- read.csv("data/CM/ID136/ID136.csv")
  raw1Matrix <- left_join(geneNames, rawMatrix,  by = c("ensembl_gene_id"="gene_id")) 
  cMatrix <- na.omit(raw1Matrix)
  countMatrix <- cMatrix[, -1]
  
  
  #aggregate rows that have the same gene names
  #countMatrix <- aggregate(x = countMatrix[-1], by = list(countMatrix$gene), FUN = sum)
  return(countMatrix)
}



