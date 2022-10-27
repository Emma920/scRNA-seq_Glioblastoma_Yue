library(tidyverse)
library(plyr)
library(dplyr)
library(data.table)
#setwd("~/hdd")


knit_aggregate_GT <- function(dir, id) {
  #setwd(dir)
  dataList <- list.files(pattern= id)
  rawMatrix <- as.data.frame(read.table(dataList[1],header=TRUE)[-1:-4, -2:-3])
  colnames(rawMatrix) <- c("gene_id", str_sub(dataList[1], 1, 2))
  for(i in 1:(length(dataList)-1)){
    rawMatrix1 <- as.data.frame(read.table(dataList[i+1],header=TRUE)[-1:-4, -2:-3])
    colnames(rawMatrix1) <- c("gene_id", str_sub(dataList[i+1], 1, 2))
    rawMatrix<- join(rawMatrix, rawMatrix1, by=c("gene_id"))
    print(i)
  }
  
  
  # geneNames <- read.table("/yue/data/text/20200313/gene_name") 
  # raw1Matrix <- left_join(geneNames, rawMatrix,  by = c("ensembl_gene_id"="gene_id")) 
  # cMatrix <- na.omit(raw1Matrix)
  # countMatrix <- cMatrix[, -1]
  
  #aggregate rows that have the same gene names
  #countMatrix <- aggregate(x = countMatrix[-1], by = list(countMatrix$gene), FUN = sum)
  #return(countMatrix)
  return(rawMatrix)
}

#countMatrix <- knit_aggregate("/home/hdd/yue/data/aligned/Darmanis/kallisto/", "*.tsv")
#write.table(countMatrix, file="/home/hdd/yue/data/CM/Darmanis/countMatrix")





#make a evaluation on genes ######################################################################################
#QC <- data.frame(cell_counts = rowSums(countMatrix_aggregate[, -1]!=0), sum_genes = rowSums(countMatrix_aggregate[, -1]))
#ggplot(QC, 
#       aes(x=cell_counts,
#           y=sum_genes)) +
#  geom_point() +
#  scale_y_continuous(trans = "log10") +
#  scale_x_continuous(trans = "log10")

#ggplot(QC,
#       aes(x=cell_counts)) +
#  geom_histogram(bins = 50) +
#  scale_x_log10()

#filter out all genes that have less than log10(x)=8 cells expressed
#Filter0 <- QC$cell_counts > 10
#countMatrix00 <- countMatrix[Filter0, ]
#write.table(countMatrix00, file="countMatrix00")

#filter out genes that have less than 100 estimated read
#row_sum <- apply(countMatrix_aggregate[, -1], 1, sum)
#countMatrix_aggregate1 <- countMatrix_aggregate[row_sum >= 100, ]
####################################################################################################################

#writing table######################################################################################################
#write.table(rawMatrix, file="rawMatrix")
#l <- read_tsv('SRR3934349_quant.tsv')
#s <- read_tsv('SRR3934351_quant.tsv')
#d <- read_tsv('hdd/alex/ncbi/public/sra/output/Trimmed/aligned/SRR3934352_quant.tsv')
#counts <- pull(l, est_counts)
#names <- pull(l, target_id)
#matrix <- bind_cols(names, counts)
#total <- merge(l,s,by.x="target_id", by.y="target_id") 
#read_tsv('transcripts_to_genes.txt', col_names = FALSE) %>% set_names(c('transcript', 'gene', 'name'))->t
####################################################################################################################

#other comments#####################################################################################################
#row_sum<-apply(head(rawMatrix, n=1000)[, -1], 1, sum)
#X <- head(rawMatrix, n=1000)
#row_sum<-apply(X[, -1], 1, sum)
#X[row_sum != 0, ]
#X[row_sum != 0, ] %>% View()
#row_sum<-apply(X[, -1], 1, sum)
#length(row_sum)
#dim(X)
#head(row_sum[row_sum != 0])
#hist(row_sum)
#summary(row_sum)
#sort(row_sum)
#order(row_sum)
#rev(order(row_sum))
#largest.sum.first <- rev(order(row_sum))
#X[largest.sum.first, ] %>% View
