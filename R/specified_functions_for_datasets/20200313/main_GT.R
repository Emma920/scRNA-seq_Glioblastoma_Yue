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

source("/home/hdd/yue/code/R/20200313/biomart.R")
gene.name <- biomart.name()
write.table(gene.name, file="/home/hdd/yue/data/text/20200313/gene_name")


source("yue/code/R/20200313/load_files_GT.R")
countMatrix_GT <- knit_aggregate_GT("yue/data/aligned/20200313", "*ReadsPerGene.out.tab")
write.table(countMatrix_GT, file="data/CM/20200313/countMatrix_GT")


source("/home/hdd/yue/code/R/20200313/load_data_GT.R")
sce <-load.data.GT("/home/hdd/yue/data/CM/20200313/countMatrix_GT", "/home/hdd/yue/data/text/20200313/20200313_meta.csv" )

isSpike(sce, "ERCC") <- grepl("^ERCC-", rownames(sce))
isSpike(sce, "MT") <- rownames(sce) %in% c("MT-TF",   "MT-RNR1", "MT-TV", "MT-RNR2", "MT-TL1",  "MT-ND1", 
                                           "MT-TI",   "MT-TQ",   "MT-TM",  "MT-ND2",  "MT-TW",  "MT-TA",   
                                           "MT-TN",   "MT-TC",  "MT-TY" ,  "MT-CO1",  "MT-TS1",  "MT-TD",
                                           "MT-CO2",  "MT-TK",   "MT-ATP8", "MT-ATP6", "MT-CO3",  "MT-TG",   
                                           "MT-ND3",  "MT-TR",   "MT-ND4L", "MT-ND4",  "MT-TH",   "MT-TS2",
                                           "MT-TL2",  "MT-ND5",  "MT-ND6", "MT-TE",   "MT-CYB",  "MT-TT", "MT-TP")

sce <- calculateQCMetrics(
  sce,
  feature_controls = list(
    ERCC = isSpike(sce, "ERCC"), 
    MT = isSpike(sce, "MT")
  )
)
hist(
  sce$total_counts,
  breaks = 100
)
abline(v = 25000, col = "red")

hist(
  sce$total_features_by_counts,
  breaks = 100
)
abline(v = 7000, col = "red")

plotColData(
  sce,
  x = "total_features_by_counts",
  y = "pct_counts_ERCC"
)
