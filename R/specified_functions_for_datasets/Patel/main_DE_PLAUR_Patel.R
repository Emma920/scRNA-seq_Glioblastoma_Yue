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


source("/home/hdd/yue/code/R/Darmanis/correlation_matrix.R")

#LOAD FILES####################################################################################
source("/home/hdd/yue/code/R/Darmanis/load_files.R")
#countMatrix <- knit_aggregate("/home/hdd/yue/data/aligned/Darmanis/kallisto/", "*.tsv")
#write.table(countMatrix, file="/home/hdd/yue/data/CM/Darmanis/countMatrix")
source("/home/hdd/yue/code/R/Darmanis/load_STARfiles.R")
#countMatrix_S <- knit_STARaggregate("/home/hdd/yue/data/aligned/Patel/STAR/", "*ReadsPerGene.out.tab")
#write.table(countMatrix_S, file="/home/hdd/yue/data/CM/Patel/countMatrix_S")

#LOAD DATA#####################################################################################
source("/home/hdd/yue/code/R/Darmanis/load_data.R")
source("/home/hdd/yue/code/R/Darmanis/pre_processing.R")
source("/home/hdd/yue/code/R/Darmanis/pre_processing_lib.R")
source("/home/hdd/yue/code/R/Darmanis/pre_processing_tpm.R")


sceP <-load.data("/home/hdd/yue/data/CM/Patel/countMatrix_S", "/home/hdd/yue/data/text/ncbi_Acc_list/Patel/Patel_meta.csv" )

#OVERVIEW########################################################################################
sceP.total <- pre.processing(sceP)
source("/home/hdd/yue/code/R/Darmanis/interested_genes.R")
genes_P <- row.names(counts(sceP.total))
source("/home/hdd/yue/code/R/Darmanis/dimension_reduction.R")
sceP.total <- dimension.reduction(genes_P, sceP.total)
plotReducedDim(sceP.total, dimred="TSNE", colour_by = "patient_id" )


###############################################################################################
#WITH EACH PATIENT#############################################################################
###############################################################################################

#SUBSET SCE WITH PATIENTS ID###################################################################
sceP.MGH26.total <- pre.processing.lib(subset(sceP, , patient_id=="MGH26"))
sceP.MGH28.total <- pre.processing.lib(subset(sceP, , patient_id=="MGH28"))
sceP.MGH29.total <- pre.processing.lib(subset(sceP, , patient_id=="MGH29"))
sceP.MGH30.total <- pre.processing.lib(subset(sceP, , patient_id=="MGH30"))
sceP.MGH31.total <- pre.processing.lib(subset(sceP, , patient_id=="MGH31"))

genes_P.MGH26 <- row.names(counts(sceP.MGH26.total))
genes_P.MGH28 <- row.names(counts(sceP.MGH28.total))
genes_P.MGH29 <- row.names(counts(sceP.MGH29.total))
genes_P.MGH30 <- row.names(counts(sceP.MGH30.total))
genes_P.MGH31 <- row.names(counts(sceP.MGH31.total))


#DIMENSION REDUCTION ON EACH PATIENT CELLS######################################################
source("/home/hdd/yue/code/R/Darmanis/dimension_reduction.R")

sceP.MGH26.total <- dimension.reduction(genes_P.MGH26, sceP.MGH26.total)
sceP.MGH28.total <- dimension.reduction(genes_P.MGH28, sceP.MGH28.total)
sceP.MGH29.total <- dimension.reduction(genes_P.MGH29, sceP.MGH29.total)
sceP.MGH30.total <- dimension.reduction(genes_P.MGH30, sceP.MGH30.total)
sceP.MGH31.total <- dimension.reduction(genes_P.MGH31, sceP.MGH31.total)

plotReducedDim(sceP.MGH26.total, dimred="TSNE", colour_by = "RGCC")
plotReducedDim(sceP.MGH28.total, dimred="TSNE", colour_by = "PROM1")
plotReducedDim(sceP.MGH29.total, dimred="TSNE", colour_by = "RGCC")
plotReducedDim(sceP.MGH30.total, dimred="TSNE", colour_by = "RGCC")
plotReducedDim(sceP.MGH31.total, dimred="TSNE", colour_by = "PLAUR")


logcountsP_PLAUR31 <- logcounts(sceP.MGH31.total["PLAUR",])
ggplot(data.frame(t(logcountsP_PLAUR31)), aes(x=PLAUR)) +
  geom_histogram(bins = 100)+
  scale_y_log10()+
  geom_vline(xintercept=max(logcountsP_PLAUR31["PLAUR",])-0.75*(max(logcountsP_PLAUR31["PLAUR",])-min(logcountsP_PLAUR31["PLAUR",])),
             color="red", linetype="dashed", size=1)
factor_test31 <- factor(logcounts(sceP.MGH31.total["PLAUR",]) > max(logcountsP_PLAUR31["PLAUR",])-0.75*(max(logcountsP_PLAUR31["PLAUR",])-min(logcountsP_PLAUR31["PLAUR",])))

source("/home/hdd/yue/code/R/Neftel/edgeR.R")
result31_PLAUR_edgeR <- edgeR(sceP.MGH31.total, "PLAUR")
write.table(result31_PLAUR_edgeR, file="/home/hdd/yue/data/output/result31_PLAUR_edgeR")

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