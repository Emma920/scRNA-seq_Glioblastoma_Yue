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
library(RColorBrewer)
library(SCnorm)


source("/home/hdd/yue/code/R/Darmanis/correlation_matrix.R")

#LOAD FILES####################################################################################
source("/home/hdd/yue/code/R/Darmanis/load_files.R")
#countMatrix <- knit_aggregate("/home/hdd/yue/data/aligned/Darmanis/kallisto/", "*.tsv")
#write.table(countMatrix, file="/home/hdd/yue/data/CM/Darmanis/countMatrix")
#source("/home/hdd/yue/code/R/Darmanis/load_STARfiles.R")
#countMatrix_S <- knit_STARaggregate("/home/hdd/yue/data/aligned/Darmanis/STAR/", "*ReadsPerGene.out.tab")
#write.table(countMatrix_S, file="/home/hdd/yue/data/CM/Darmanis/countMatrix_S")

#LOAD DATA#####################################################################################
source("/home/hdd/yue/code/R/Darmanis/load_data.R")
source("/home/hdd/yue/code/R/Darmanis/pre_processing.R")
source("/home/hdd/yue/code/R/Darmanis/pre_processing_lib.R")


sce <-load.data("/home/hdd/yue/data/CM/Darmanis/countMatrix_S", "/home/hdd/yue/data/text/ncbi_Acc_list/Darmanis/neo_metadata.csv" )
sce.BT_S2 <- subset(sce, , patient_id=="BT_S2")

#FILTER CELLS AND GENES###############################################################################
source("/home/hdd/yue/code/R/Darmanis/filter.R")
sce.BT_S2.total <- filter(sce.BT_S2)
source("/home/hdd/yue/code/R/Darmanis/pre_processing.R")
sce.BT_S2.pre <- pre.processing(sce.BT_S2)
########################################################################################################

Conditions = rep(c(1), each= 431)
par(mfrow=c(2,2))
DataNorm <- SCnorm(Data = sce.BT_S2.total, Conditions = Conditions,PrintProgressPlots = TRUE,
                   FilterCellNum = 10, K = 4,NCores=10, reportSF = TRUE)

NormalizedData <- SingleCellExperiment::normcounts(DataNorm)
plotCountDepth(Data = sce.BT_S2.total,
               NormalizedData = logNormalizedData,
               Conditions = Conditions,FilterCellProportion = .1, NCores=3)

logNormalizedData <- log(NormalizedData+5.66)

genes.BT_S2 <- row.names(counts(sce.BT_S2.total))

source("/home/hdd/yue/code/R/Darmanis/interested_genes.R")
cscM_2 <- interested.genes("cscM", sce.BT_S2.total)
diffA_2 <- interested.genes("diffA", sce.BT_S2.total)
angioA_2 <- interested.genes("angioA", sce.BT_S2.total)
cycleA_2 <- interested.genes("cycleA", sce.BT_S2.total)
migrHM_2 <- interested.genes("migrHM", sce.BT_S2.total)
#nscA <- interested.genes("angioM", sce.total)
#nscHM <- interested.genes("angioM", sce.total)
srHM_2 <- interested.genes("srHM", sce.BT_S2.total)

set.seed(100)
rand <- sample(nrow(NormalizedData), 200)
sum.all <- as.matrix(data.frame(csc=colSums(cpm(sce.BT_S2.pre[cscM_2, ])),#/sumall, 
                                diff=colSums(cpm(sce.BT_S2.pre[diffA_2, ])),#/sumall,
                                migr=colSums(cpm(sce.BT_S2.pre[migrHM_2, ])),#/sumall, 
                                angio=colSums(cpm(sce.BT_S2.pre[angioA_2, ]))#/sumall
))
#                                cycle = colSums(logcounts(sce.BT_S2.pre[cycleA_2, ]))

#                                random =colSums(logcounts(sce.BT_S2.pre[rand, ])),

sum.all <- as.matrix(data.frame(csc=colSums(NormalizedData[cscM_2, ]),#/sumall, 
                                diff=colSums(NormalizedData[diffA_2, ]),#/sumall,
                                migr=colSums(NormalizedData[migrHM_2, ]),#/sumall, 
                                angio=colSums(NormalizedData[angioA_2, ])
#/sumall
))

source("/home/hdd/yue/code/R/Darmanis/sum_sigma.R")
sumsigma <- as.matrix(data.frame(csc=sum.sigma(cscM_2, sce.BT_S2.pre),#/sumall, 
                                 diff=sum.sigma(diffA_2, sce.BT_S2.pre),#/sumall,
                                 migr=sum.sigma(migrHM_2, sce.BT_S2.pre),#/sumall, 
                                 angio=sum.sigma(angioA_2, sce.BT_S2.pre),#/sumall
                                 random =sum.sigma(rand, sce.BT_S2.pre)))
#cycle = colSums(NormalizedData[cycleA_2, ])
#random =colSums(logNormalizedData[rand, ])
#random =colSums(logNormalizedData[rand, ]),

library(pheatmap)
pm <- pheatmap(t(sumsigma), scale="row", show_colnames = FALSE,
               color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100), 
               cellwidth=1.5, cellheight =60)

