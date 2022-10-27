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
library(monocle3)







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
source("/home/hdd/yue/code/R/Darmanis/pre_processing_tpm.R")


sce <-load.data("/home/hdd/yue/data/CM/Darmanis/countMatrix_S", "/home/hdd/yue/data/text/ncbi_Acc_list/Darmanis/neo_metadata.csv" )
sce.total <- pre.processing.lib(sce)
matrixtest <- counts(sce[])

countMatrix_aggregate <- read.table("/home/hdd/yue/data/CM/Darmanis/countMatrix_S", header = TRUE) #countMatrix_aggregate was made by collecting est_counts from each cell and transcripts with the same gene names were aggregated. None of the cells or genes were omited
countMatrix_aggregate1 <- countMatrix_aggregate[!duplicated(countMatrix_aggregate[1]), ]
data <- as.matrix(data.frame(countMatrix_aggregate1[,-1], row.names=countMatrix_aggregate1[,1]))
meta_data1 <- read.csv("/home/hdd/yue/data/text/ncbi_Acc_list/Darmanis/neo_metadata.csv")
meta_data <- as.matrix(data.frame(meta_data1[, -1], row.names = meta_data1[,1]))

meta <- meta_data[colnames(data),]
cds <- new_cell_data_set(data, cell_metadata = meta)



#INTERESTED GENES##############################################################################
source("/home/hdd/yue/code/R/Darmanis/interested_genes.R")
genes <- row.names(counts(sce.total))

#DIMENSIONAL REDUCTION ON SELECTED GENE SETS###################################################
source("/home/hdd/yue/code/R/Darmanis/dimension_reduction.R")
sce.total <- dimension.reduction(genes, sce.total)
plotReducedDim(sce.total, dimred="TSNE", colour_by = "CLOCK" )


###############################################################################################
#WITH EACH PATIENT#############################################################################
###############################################################################################

#SUBSET SCE WITH PATIENTS ID###################################################################
sce.BT_S2.total <- pre.processing.lib(subset(sce, , patient_id=="BT_S2"))
sce.BT_S1.total <- pre.processing.lib(subset(sce, , patient_id=="BT_S1"))
sce.BT_S4.total <- pre.processing.lib(subset(sce, , patient_id=="BT_S4"))
sce.BT_S6.total <- pre.processing.lib(subset(sce, , patient_id=="BT_S6"))

genes.BT_S2 <- row.names(counts(sce.BT_S2.total))
genes.BT_S1 <- row.names(counts(sce.BT_S1.total))
genes.BT_S4 <- row.names(counts(sce.BT_S4.total))
genes.BT_S6 <- row.names(counts(sce.BT_S6.total))


#DIMENSION REDUCTION ON EACH PATIENT CELLS######################################################
source("/home/hdd/yue/code/R/Darmanis/dimension_reduction.R")

sce.BT_S2.total <- dimension.reduction(genes.BT_S2, sce.BT_S2.total)
sce.BT_S1.total <- dimension.reduction(genes.BT_S1, sce.BT_S1.total)
sce.BT_S4.total <- dimension.reduction(genes.BT_S4, sce.BT_S4.total)  
sce.BT_S6.total <- dimension.reduction(genes.BT_S6, sce.BT_S6.total)   

plotReducedDim(sce.BT_S2.total, dimred="TSNE", colour_by ="CTNNB1")

TP53 <- c(logcounts(sce.BT_S2.total["TP53",]))
PLAUR <-  c(logcounts(sce.BT_S2.total["PLAUR",])) 
CRY1 <- c(logcounts(sce.BT_S2.total["CRY1",]))
CRY2 <-  c(logcounts(sce.BT_S2.total["CRY2",])) 
POLR2A <- c(logcounts(sce.BT_S2.total["POLR2A",]))
HDAC1 <-  c(logcounts(sce.BT_S2.total["HDAC1",])) 
SOX2 <-  c(logcounts(sce.BT_S2.total["SOX2",])) 
CTNNB1 <-  c(logcounts(sce.BT_S2.total["CTNNB1",])) 


test <- data.frame(SOX2, PLAUR)
ggplot(test, aes(x=SOX2, y=PLAUR)) + 
  geom_point()

logcounts_PLAUR2 <- logcounts(sce.BT_S2.total["PLAUR",])
ggplot(data.frame(t(logcounts_PLAUR2)), aes(x=PLAUR)) +
  geom_histogram(bins = 100)+
  scale_y_log10()+
  geom_vline(xintercept=max(logcounts_PLAUR2["PLAUR",])-0.75*(max(logcounts_PLAUR2["PLAUR",])-min(logcounts_PLAUR2["PLAUR",])),
             color="red", linetype="dashed", size=1)

factor_test2 <- factor(logcounts(sce.BT_S2.total["PLAUR",]) > max(logcounts_PLAUR2["PLAUR",])-0.75*(max(logcounts_PLAUR2["PLAUR",])-min(logcounts_PLAUR2["PLAUR",])))

logcounts_PLAUR1 <- logcounts(sce.BT_S1.total["PLAUR",])
ggplot(data.frame(t(logcounts_PLAUR1)), aes(x=PLAUR)) +
  geom_histogram(bins = 100)+
  scale_y_log10()+
  geom_vline(xintercept=max(logcounts_PLAUR1["PLAUR",])-0.75*(max(logcounts_PLAUR1["PLAUR",])-min(logcounts_PLAUR1["PLAUR",])),
             color="red", linetype="dashed", size=1)


library(MAST)
fData <- data.frame(names = rownames(logcounts(sce.BT_S2.total)))
rownames(fData) <- rownames(logcounts(sce.BT_S2.total))
cData <- data.frame(cond = factor_test2)
rownames(cData) <- colnames(logcounts(sce.BT_S2.total))

obj <- FromMatrix(as.matrix(logcounts(sce.BT_S2.total)), cData, fData)
#colData(obj)$cngeneson <- scale(colSums(assay(obj) > 0))
cond <- factor(colData(obj)$cond)

# Model expression as function of condition & number of detected genes
zlmCond <- zlm(~ cond, obj) 
summaryCond <- summary(zlmCond)
summaryDt <- summaryCond$datatable


source("/home/hdd/yue/code/R/Neftel/edgeR.R")
result6_PLAUR_edgeR <- edgeR(sce.BT_S6.total, "PLAUR")
write.table(result6_PLAUR_edgeR, file="/home/hdd/yue/data/output/result6_PLAUR_edgeR")


#CENTERING THE CELLS##############################################################################
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
#########################################################################################################


