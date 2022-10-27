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


###############################################################################################
#WITH EACH PATIENT#############################################################################
###############################################################################################

#SUBSET SCE WITH PATIENTS ID###################################################################
sce.BT_S2.total <- pre.processing(subset(sce, , patient_id=="BT_S2"))
sce.BT_S1.total <- pre.processing(subset(sce, , patient_id=="BT_S1"))
sce.BT_S4.total <- pre.processing(subset(sce, , patient_id=="BT_S4"))
sce.BT_S6.total <- pre.processing(subset(sce, , patient_id=="BT_S6"))




#DIMENSION REDUCTION ON EACH PATIENT CELLS######################################################
source("/home/hdd/yue/code/R/Darmanis/dimension_reduction.R")

sce.BT_S2.total <- dimension.reduction(genes.BT_S2, sce.BT_S2.total)
sce.BT_S1.total <- dimension.reduction(genes.BT_S1, sce.BT_S1.total)
sce.BT_S4.total <- dimension.reduction(genes.BT_S4, sce.BT_S4.total)
sce.BT_S6.total <- dimension.reduction(genes.BT_S6, sce.BT_S6.total)

plotReducedDim(sce.BT_S2.total, dimred="TSNE", colour_by = "SERPINE1" )




#ALL_GENE LISTS FROM EACH PATIENT###############################################################
genes.BT_S2 <- row.names(counts(sce.BT_S2.total))
genes.BT_S1 <- row.names(counts(sce.BT_S1.total))
genes.BT_S4 <- row.names(counts(sce.BT_S4.total))
genes.BT_S6 <- row.names(counts(sce.BT_S6.total))
#INTERESTED GENE SETS
cscM_2 <- interested.genes("cscM", sce.BT_S2.total)
diffA_2 <- interested.genes("diffA", sce.BT_S2.total)
angioA_2 <- interested.genes("angioA", sce.BT_S2.total)
cycleA_2 <- interested.genes("cycleA", sce.BT_S2.total)
migrHM_2 <- interested.genes("migrHM", sce.BT_S2.total)
#nscA <- interested.genes("angioM", sce.total)
#nscHM <- interested.genes("angioM", sce.total)
srHM_2 <- interested.genes("srHM", sce.BT_S2.total)
cscM_1 <- interested.genes("cscM", sce.BT_S1.total)
diffA_1 <- interested.genes("diffA", sce.BT_S1.total)
angioA_1 <- interested.genes("angioA", sce.BT_S1.total)
migrHM_1 <- interested.genes("migrHM", sce.BT_S1.total)
srHM_1 <- interested.genes("srHM", sce.BT_S1.total)
cscM_4 <- interested.genes("cscM", sce.BT_S4.total)
diffA_4 <- interested.genes("diffA", sce.BT_S4.total)
angioA_4 <- interested.genes("angioA", sce.BT_S4.total)
migrHM_4 <- interested.genes("migrHM", sce.BT_S4.total)
srHM_4 <- interested.genes("srHM", sce.BT_S4.total)
cscM_6 <- interested.genes("cscM", sce.BT_S6.total)
diffA_6 <- interested.genes("diffA", sce.BT_S6.total)
angioA_6 <- interested.genes("angioA", sce.BT_S6.total)
migrHM_6 <- interested.genes("migrHM", sce.BT_S6.total)
srHM_6 <- interested.genes("srHM", sce.BT_S6.total)

#SUM SIGMA TSNE#####################################################################################
source("/home/hdd/yue/code/R/Darmanis/sum_sigma.R")
cscM_sigma <- sum.sigma(cscM, sce.BT_S2.total)
sce.BT_S2.total$cscM_sigma <- runif(ncol(sce.BT_S2.total))
sce.BT_S2.total$cscM_sigma <- cscM_sigma
plotReducedDim(sce.BT_S2.total, dimred="TSNE", colour_by = "cscM_sigma" )

cscM_sigma <- sum.sigma(cscM, sce.BT_S2.total)
sce.BT_S2.total$cscM_sigma <- runif(ncol(sce.BT_S2.total))
sce.BT_S2.total$cscM_sigma <- cscM_sigma
plotReducedDim(sce.BT_S2.total, dimred="TSNE", colour_by = "cscM_sigma" )


#DIMENSION REDUCTION ON SELECTED GENES#########################################################
sce.BT_S2.cscM <- sce.BT_S2.total[cscM, ]
sce.BT_S2.cscM <- dimension.reduction(cscM, sce.BT_S2.cscM)

cscM.BT_S2.sum <- colSums(logcounts(sce.BT_S2.cscM))
sce.BT_S2.cscM$cscM_sum <- runif(ncol(sce.BT_S2.total))
sce.BT_S2.cscM$cscM_sum <- cscM.BT_S2.sum

sce.BT_S2.cscM$cscM_diffA_dcor_sigma <- runif(ncol(sce.BT_S2.cscM))
sce.BT_S2.cscM$cscM_diffA_dcor_sigma <- csc_diffA_sigma

plotReducedDim(sce.BT_S2.cscM, dimred="TSNE", colour_by = "cscM_diffA_dcor_sigma" )

g <- buildSNNGraph(sce.BT_S2.cscM, k=20, use.dimred = 'PCA') #using the top PCs
clust.louvain <- igraph::cluster_louvain(g)$membership
sce.BT_S2.cscM$label <- factor(clust.louvain)
markers.BT_S2 <- findMarkers(sce.BT_S2.cscM, groups=sce.BT_S2.cscM$label)
rm(g, clust.louvain)
plotReducedDim(sce.BT_S2.cscM, dimred="TSNE", colour_by = "label")


#SUM FEATURES PLOTS ON TSNE########################################################################
cscM.BT_S2.sum <- colSums(logcounts(sce.BT_S2.total[cscM, ]))
sce.BT_S2.total$cscM_sum <- runif(ncol(sce.BT_S2.total))
sce.BT_S2.total$cscM_sum <- cscM.BT_S2.sum

plotReducedDim(sce.BT_S2.total, dimred="TSNE", colour_by = "cscM_sum" )

rand <- sample(nrow(sce.BT_S2.total), 100)
sce.rand <- sce.BT_S2.total[rand,]
rand.sum <- colSums(logcounts(sce.BT_S2.total))
sce.BT_S2.total$rand_sum <- runif(ncol(sce.BT_S2.total))
sce.BT_S2.total$rand_sum <- rand.sum
plotReducedDim(sce.BT_S2.total, dimred="TSNE", colour_by = "rand_sum" )