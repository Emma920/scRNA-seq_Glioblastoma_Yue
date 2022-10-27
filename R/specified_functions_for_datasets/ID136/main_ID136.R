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

# source("/home/hdd/yue/code/R/20200313/biomart.R")
# gene.name <- biomart.name()
# write.table(gene.name, file="/home/hdd/yue/data/text/20200313/gene_name")


# source("code/R/ID136/load_files.R")
# countMatrix_GT <- assign_names("data/CM/ID136/ID136.csv")
# write.table(countMatrix_GT, file="data/CM/ID136/ID136_CM")


source("code/R/ID132/load_data.R")
sce_ID136 <-load.data.GT("data/CM/ID136/ID136_CM", "data/text/ID136/ID136_meta.csv" )

isSpike(sce_ID136, "ERCC") <- grepl("^ERCC-", rownames(sce_ID136))
isSpike(sce_ID136, "MT") <- rownames(sce_ID136) %in% c("MT-TF",   "MT-RNR1", "MT-TV", "MT-RNR2", "MT-TL1",  "MT-ND1", 
                                           "MT-TI",   "MT-TQ",   "MT-TM",  "MT-ND2",  "MT-TW",  "MT-TA",   
                                           "MT-TN",   "MT-TC",  "MT-TY" ,  "MT-CO1",  "MT-TS1",  "MT-TD",
                                           "MT-CO2",  "MT-TK",   "MT-ATP8", "MT-ATP6", "MT-CO3",  "MT-TG",   
                                           "MT-ND3",  "MT-TR",   "MT-ND4L", "MT-ND4",  "MT-TH",   "MT-TS2",
                                           "MT-TL2",  "MT-ND5",  "MT-ND6", "MT-TE",   "MT-CYB",  "MT-TT", "MT-TP")

sce_ID136 <- calculateQCMetrics(
  sce_ID136,
  feature_controls = list(
    ERCC = isSpike(sce, "ERCC"), 
    MT = isSpike(sce, "MT")
  )
)
hist(
  sce_ID136$total_counts,
  breaks = 100
)
abline(v = 25000, col = "red")

hist(
  sce_ID136$total_features_by_counts, 
  xlab = "Number of detected genes", labels = FALSE,
  main = "ID136",
  breaks = 100
)
abline(v = 4000, col = "red")

plotColData(
  sce_ID136,
  x = "total_features_by_counts",
  y = "pct_counts_MT"
)
################################################################################
source("code/R/ID132/pre_processing_lib.R")
sce.ID136 <- pre.processing.lib(sce_ID136)
source("code/R/ID132/pre_processing_ERCC.R")
sce.ID136 <- pre.processing.ERCC(sce_ID136)
################################################################################

source("code/R/Darmanis/clusting.R")
sce.ID136@colData@listData[["label"]] <- clusting(sce.ID136)

genes_ID136 <- row.names(logcounts(sce.ID136))
source("code/R/Darmanis/dimension_reduction.R")
sce.ID136<- dimension.reduction(genes_ID136, sce.ID136)
plotReducedDim(sce.ID136, dimred="UMAP", colour_by = "PROM1")

markers.ID136 <- findMarkers(sce.ID136, sce.ID136@colData@listData[["label"]], pval.type="all")



##################################################################################
plot.PLAUR <- function() {
  logcounts_ID132 <- logcounts(sce.ID132["PLAUR",])
  ggplot(data.frame(t(logcounts_ID132)), aes(x=PLAUR)) +
    geom_histogram(bins = 100)+
    scale_y_log10()+
    geom_vline(xintercept=max(logcounts_ID132["PLAUR",])-0.75*(max(logcounts_ID132["PLAUR",])-min(logcounts_ID132["PLAUR",])),
               color="red", linetype="dashed", size=1)
}

source("code/R/Patel/DE_PLAUR.R")
DE_rankByP <- DE.PLAUR(sce.ID132)
rownames_500 <- rownames(DE_rankByP)[1:200]
