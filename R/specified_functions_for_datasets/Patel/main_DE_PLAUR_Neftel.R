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
source("/home/hdd/yue/code/R/Darmanis/load_data_Neftel.R")

sceN <-load.data.Neftel("/home/hdd/yue/data/CM/Neftel/GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv", "/home/hdd/yue/data/text/ncbi_Acc_list/Neftel/meta_Neftel.csv" )
sceN.total <- logNormCounts(sceN, pseudo_count = 1, exprs_values = "tpm")

#OVERVIEW########################################################################################
genes_N <- row.names(tpm(sceN.total))
source("/home/hdd/yue/code/R/Darmanis/dimension_reduction.R")
sceN.total <- dimension.reduction(genes_N, sceN.total)
plotReducedDim(sceN.total, dimred="TSNE", colour_by = "tumour.name" )
tumour_names <- unique(sceN.total[["tumour.name"]])

###############################################################################################
#WITH EACH PATIENT#############################################################################
###############################################################################################

#SUBSET SCE WITH PATIENTS ID###################################################################
sceN.MGH101.total <- subset(sceN.total, , tumour.name=="MGH101")
sceN.MGH100.total <- subset(sceN.total, , tumour.name=="MGH100")
sceN.MGH102.total <- subset(sceN.total, , tumour.name=="MGH102")
sceN.MGH104.total <- subset(sceN.total, , tumour.name=="MGH104")
sceN.MGH105.total <- subset(sceN.total, , tumour.name=="MGH105")
sceN.MGH106.total <- subset(sceN.total, , tumour.name=="MGH106")
sceN.MGH110.total <- subset(sceN.total, , tumour.name=="MGH110")
sceN.MGH113.total <- subset(sceN.total, , tumour.name=="MGH113")
sceN.MGH115.total <- subset(sceN.total, , tumour.name=="MGH115")
sceN.MGH121.total <- subset(sceN.total, , tumour.name=="MGH121")
sceN.MGH122.total <- subset(sceN.total, , tumour.name=="MGH122")
sceN.MGH124.total <- subset(sceN.total, , tumour.name=="MGH124")
sceN.MGH125.total <- subset(sceN.total, , tumour.name=="MGH125")
sceN.BT749.total <- subset(sceN.total, , tumour.name=="BT749")
sceN.BT771.total <- subset(sceN.total, , tumour.name=="BT771")
sceN.BT830.total <- subset(sceN.total, , tumour.name=="BT830")
sceN.MGH85.total <- subset(sceN.total, , tumour.name=="MGH85")
sceN.BT1160.total <- subset(sceN.total, , tumour.name=="BT1160")
sceN.BT1187.total <- subset(sceN.total, , tumour.name=="BT1187")
sceN.BT786.total <- subset(sceN.total, , tumour.name=="BT786")
sceN.BT920.total <- subset(sceN.total, , tumour.name=="BT920")
sceN.MGH128.total <- subset(sceN.total, , tumour.name=="MGH128")
sceN.MGH129.total <- subset(sceN.total, , tumour.name=="MGH129")
sceN.MGH136.total <- subset(sceN.total, , tumour.name=="MGH136")
sceN.MGH143.total <- subset(sceN.total, , tumour.name=="MGH143")
sceN.MGH151.total <- subset(sceN.total, , tumour.name=="MGH151")
sceN.MGH152.total <- subset(sceN.total, , tumour.name=="MGH152")
sceN.MGH66.total <- subset(sceN.total, , tumour.name=="MGH66")

genes_N.MGH101 <- row.names(tpm(sceN.MGH101.total))
write.table(genes_N.MGH101, file="/home/hdd/yue/data/output/genes_N.MGH101")

genes_N.MGH100 <- row.names(tpm(sceN.MGH100.total))
genes_P.MGH102 <- row.names(tpm(sceN.MGH102.total))
genes_P.MGH30 <- row.names(counts(sceP.MGH30.total))
genes_P.MGH31 <- row.names(counts(sceP.MGH31.total))


#DIMENSION REDUCTION ON EACH PATIENT CELLS######################################################
source("/home/hdd/yue/code/R/Darmanis/dimension_reduction.R")

sceN.MGH101.total <- dimension.reduction(genes_N.MGH101, sceN.MGH101.total)
sceN.MGH100.total <- dimension.reduction(genes_N.MGH100, sceN.MGH100.total)
sceP.MGH29.total <- dimension.reduction(genes_P.MGH29, sceP.MGH29.total)
sceP.MGH30.total <- dimension.reduction(genes_P.MGH30, sceP.MGH30.total)
sceP.MGH31.total <- dimension.reduction(genes_P.MGH31, sceP.MGH31.total)

plotReducedDim(sceN.MGH101.total, dimred="TSNE", colour_by = "NANOG")
plotReducedDim(sceN.MGH100.total, dimred="TSNE", colour_by = "PLAUR")
plotReducedDim(sceP.MGH29.total, dimred="TSNE", colour_by = "PLAUR")
plotReducedDim(sceP.MGH30.total, dimred="TSNE", colour_by = "PLAUR")
plotReducedDim(sceP.MGH31.total, dimred="TSNE", colour_by = "PLAUR")


plot.PLAUR <- function() {
  logcountsP_PLAUR31 <- logcounts(sceP.MGH31.total["PLAUR",])
  ggplot(data.frame(t(logcountsP_PLAUR31)), aes(x=PLAUR)) +
    geom_histogram(bins = 100)+
    scale_y_log10()+
    geom_vline(xintercept=max(logcountsP_PLAUR31["PLAUR",])-0.75*(max(logcountsP_PLAUR31["PLAUR",])-min(logcountsP_PLAUR31["PLAUR",])),
               color="red", linetype="dashed", size=1)
}

#DIFFERENTIAL EXPRESSION####################################################################################

source("/home/hdd/yue/code/R/Neftel/edgeR.R")
result106_PLAUR_edgeR <- edgeR(sceN.MGH106.total, "PLAUR")
write.table(result106_PLAUR_edgeR, file="/home/hdd/yue/data/output/result106_PLAUR_edgeR")



source("/home/hdd/yue/code/R/Neftel/monocle.R")
result101_PLAUR_monocle <- monocle(sceN.MGH101.total, "PLAUR") 

pVals_rank_monocle <- order(pVals) 
pVals_ranked_monocle <- data.frame(pVals[pVals_rank_monocle])
write.table(pVals_rank_monocle, file="/home/hdd/yue/data/output/pVals_rank_monocle")
#pVals <- p.adjust(pVals, method = "fdr")
#DE_Quality_AUC(pVals)


source("/home/hdd/yue/code/R/Patel/DE_PLAUR.R")
DE_rankByP <- DE.PLAUR(sceP.MGH28.total)
rownames_500 <- rownames(DE_rankByP)[1:200]
write.table(rownames_500, file="/home/hdd/yue/data/output/rownames_500")
write.table(DE_2, file="/home/hdd/yue/data/output/DE_2")

#CENTERING THE CELLS################################################################################
all.genes <- row.names(sce.BT_S2.total)
BT_S2.seurat <- as.Seurat(sce.BT_S2.total, counts = "counts", data = "logcounts")
BT_S2.seurat <- ScaleData(BT_S2.seurat, features = all.genes)
test2 <- BT_S2.seurat@assays[["RNA"]]@scale.data

ggplot(data.frame(PLAUR = test2["PLAUR",]), aes(x=PLAUR)) +
  geom_histogram(bins = 100)+
  scale_y_log10()+
  geom_vline(xintercept=0.25*max(test2["PLAUR",]),
             color="red", linetype="dashed", size=1)

selected_2 <- test2[,test2["PLAUR",] > 0.25*max(test2["PLAUR",])]
omited_2 <- test2[, test2["PLAUR",] <= 0.25*max(test2["PLAUR",])]
factor_test2 <- factor(test2["PLAUR",] > 0.25*max(test2["PLAUR",]))
#####################################################################################################