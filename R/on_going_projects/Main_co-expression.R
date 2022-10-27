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
library("CEMiTool")


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



#ALL_GENE LISTS FROM EACH PATIENT###############################################################
genes.BT_S2 <- row.names(counts(sce.BT_S2.total))
genes.BT_S1 <- row.names(counts(sce.BT_S1.total))
genes.BT_S4 <- row.names(counts(sce.BT_S4.total))
genes.BT_S6 <- row.names(counts(sce.BT_S6.total))

#INTERESTED GENE SETS##############################################################################
source("/home/hdd/yue/code/R/Darmanis/interested_genes.R")
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

CscDiffAngio2 <- c(cscM_2, diffA_2, angioA_2, migrHM_2, cycleA_2, srHM_2)
sce.BT_S2.features2 <- data.frame(logcounts(sce.BT_S2.total[CscDiffAngio2]))
cem <- cemitool(sce.BT_S2.features2, set_beta=3)
report <- diagnostic_report(cem)

cem2 <- cemitool(data.frame(logcounts(sce.BT_S2.total)))

find_modules(cem2)
generate_report(cem2)
hubs <- get_hubs(cem2,500)

M1_30 <- rownames(hubs[["M1"]])
common_CSC_M1 <- intersect(M1_30, diffA_2)

generate_report(cem2, directory="/home/hdd/yue/data/co_expression/cem2_report")

# write analysis results into files
write_files(cem2, directory="/home/hdd/yue/data/co_expression/cem2_files")

# save all plots
save_plots(cem2, "all", directory="/home/hdd/yue/data/co_expression/cem2_plots")