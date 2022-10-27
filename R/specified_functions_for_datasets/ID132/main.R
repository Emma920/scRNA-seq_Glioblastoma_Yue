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


#source("code/R/ID132/load_files.R")
#countMatrix_ID132 <- assign_names("data/CM/ID132/ID132.csv")
#write.table(countMatrix_ID132, file="data/CM/ID132/ID132_CM")


source("code/R/ID132/load_data.R")
sce_ID132 <-load.data.GT("data/CM/ID132/ID132_CM", "data/text/ID132/ID132_meta.csv" )


isSpike(sce_ID132, "MT") <- rownames(sce_ID132) %in% c("MT-TF",   "MT-RNR1", "MT-TV", "MT-RNR2", "MT-TL1",  "MT-ND1", 
                                           "MT-TI",   "MT-TQ",   "MT-TM",  "MT-ND2",  "MT-TW",  "MT-TA",   
                                           "MT-TN",   "MT-TC",  "MT-TY" ,  "MT-CO1",  "MT-TS1",  "MT-TD",
                                           "MT-CO2",  "MT-TK",   "MT-ATP8", "MT-ATP6", "MT-CO3",  "MT-TG",   
                                           "MT-ND3",  "MT-TR",   "MT-ND4L", "MT-ND4",  "MT-TH",   "MT-TS2",
                                           "MT-TL2",  "MT-ND5",  "MT-ND6", "MT-TE",   "MT-CYB",  "MT-TT", "MT-TP")
sce_ID132 <- calculateQCMetrics(
  sce_ID132,
  feature_controls = list(
    ERCC = isSpike(sce_ID132, "ERCC"), 
    MT = isSpike(sce_ID132, "MT")
  )
)
hist(
  sce$total_counts,
  breaks = 100
)
abline(v = 25000, col = "red")

total_features_by_counts <- sce_ID132$total_features_by_counts
hist(
  sce_ID132$total_features_by_counts, 
  xlab = "Number of detected genes", labels = FALSE,
  main = "ID132",
  breaks = 100
)
abline(v = 4000, col = "red")


plotColData(
  sce_ID132,
  x = "total_features_by_counts",
  y = "pct_counts_MT",
)

#UMAP################################################################################
source("code/R/ID132/pre_processing_lib.R")
sce.ID132 <- pre.processing.lib(sce_ID132)
#source("code/R/ID132/pre_processing_ERCC.R")
#sce.ID132 <- pre.processing.ERCC(sce)


#####################################################################################

source("code/R/Darmanis/clusting.R")
sce.ID132@colData@listData[["label"]] <- clusting(sce.ID132)

genes_ID132 <- row.names(logcounts(sce.ID132))
source("code/R/Darmanis/dimension_reduction.R")
sce.ID132<- dimension.reduction(genes_ID132, sce.ID132)
plotReducedDim(sce.ID132, dimred="UMAP", colour_by = "PROM1")

markers.ID132 <- findMarkers(sce.ID132, sce.ID132@colData@listData[["label"]], pval.type="all")





plot.PLAUR <- function() {
  logcounts_ID132 <- logcounts(sce.ID132["PLAUR",])
  ggplot(data.frame(t(logcounts_ID132)), aes(x=PLAUR)) +
    geom_histogram(bins = 100)+
    scale_y_log10()+
    geom_vline(xintercept=max(logcounts_ID132["PLAUR",])-0.75*(max(logcounts_ID132["PLAUR",])-min(logcounts_ID132["PLAUR",])),
               color="red", linetype="dashed", size=1)
}

# source("code/R/Patel/DE_PLAUR
# DE_rankByP <- DE.PLAUR(sce.ID132)
# rownames_500 <- rownames(DE_rankByP)[1:200]



### Find cdx33 cluster

umap_x <- sce.ID132@int_colData@listData[["reducedDims"]]@listData[["UMAP"]][,1]
umap_y <- sce.ID132@int_colData@listData[["reducedDims"]]@listData[["UMAP"]][,2]
cdx33_pos <- umap_y > 1 & umap_x < -2

sce.ID132@colData@listData[["label"]]
sce.ID132@colData@listData[["label"]][cdx33_pos] <- 1
sce.ID132@colData@listData[["label"]][!cdx33_pos] <- 2