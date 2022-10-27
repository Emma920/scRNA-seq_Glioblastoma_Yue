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
library(Seurat)


#LOAD DATA#####################################################################################
source("/home/hdd/yue/code/R/Tirosh/load_data_Tirosh.R")

sceT <-load.data.Tirosh("/home/hdd/yue/data/CM/Tirosh/GSE102130_K27Mproject.RSEM.vh20170621.txt", "/home/hdd/yue/data/text/ncbi_Acc_list/Tirosh/meta_Tirosh.csv" )
sceT.total <- logNormCounts(sceT, pseudo_count = 1, exprs_values = "tpm")

sceT.BCH836.total <- subset(sceT.total, , Sample =="BCH836")


source("/home/hdd/yue/code/R/Neftel/edgeR.R")
result836_PLAUR_edgeR <- edgeR(sceT.BCH836.total, "PLAUR")
write.table(result836_PLAUR_edgeR, file="/home/hdd/yue/data/output/result836_PLAUR_edgeR")

