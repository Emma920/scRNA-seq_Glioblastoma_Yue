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


#SUBSET SCE WITH PATIENTS ID###################################################################
sce.BT_S2.total <- pre.processing.lib(subset(sce, , patient_id=="BT_S2"))
sce.BT_S1.total <- pre.processing.lib(subset(sce, , patient_id=="BT_S1"))
sce.BT_S4.total <- pre.processing.lib(subset(sce, , patient_id=="BT_S4"))
sce.BT_S6.total <- pre.processing.lib(subset(sce, , patient_id=="BT_S6"))

genes.BT_S2 <- as.data.frame(row.names(counts(sce.BT_S2.total))) 
genes.BT_S1 <- as.data.frame(row.names(counts(sce.BT_S1.total)))
genes.BT_S4 <- as.data.frame(row.names(counts(sce.BT_S4.total)))
genes.BT_S6 <- as.data.frame(row.names(counts(sce.BT_S6.total)))


#DIMENSION REDUCTION ON EACH PATIENT CELLS######################################################
source("/home/hdd/yue/code/R/Darmanis/dimension_reduction.R")

sce.BT_S2.total <- dimension.reduction(genes.BT_S2, sce.BT_S2.total)
sce.BT_S1.total <- dimension.reduction(genes.BT_S1, sce.BT_S1.total)
sce.BT_S4.total <- dimension.reduction(genes.BT_S4, sce.BT_S4.total)
sce.BT_S6.total <- dimension.reduction(genes.BT_S6, sce.BT_S6.total)

plotReducedDim(sce.BT_S2.total, dimred="TSNE", colour_by ="FGF1")

Ramilowski <- read_csv("/home/hdd/yue/data/text/LR/Ramilowski.csv")
ligand <- unique((Ramilowski["Ligand.ApprovedSymbol"]))
ligand_detected <- semi_join(genes.BT_S6, ligand, by = c("row.names(counts(sce.BT_S6.total))"="Ligand.ApprovedSymbol"))
colnames(ligand_detected) <- "ligand.dt"
ligand_dt <- as.character(ligand_detected$ligand.dt)

receptor <- unique((Ramilowski["Receptor.ApprovedSymbol"]))
receptor_detected <- semi_join(genes.BT_S6, receptor, by = c("row.names(counts(sce.BT_S6.total))"="Receptor.ApprovedSymbol"))
colnames(receptor_detected) <- "receptor.dt"
receptor_dt <- as.character(receptor_detected$receptor.dt)

genes_P.MGH26 <- as.data.frame(row.names(counts(sceP.MGH26.total)))
genes_P.MGH28 <- as.data.frame(row.names(counts(sceP.MGH28.total)))
genes_P.MGH29 <- as.data.frame(row.names(counts(sceP.MGH29.total)))
genes_P.MGH30 <- as.data.frame(row.names(counts(sceP.MGH30.total)))
genes_P.MGH31 <- as.data.frame(row.names(counts(sceP.MGH31.total)))

genes_N.MGH100 <- as.data.frame(row.names(tpm(sceN.MGH100.total)))
genes_N.MGH102 <- as.data.frame(row.names(tpm(sceN.MGH102.total)))
genes_N.MGH101 <- as.data.frame(row.names(tpm(sceN.MGH101.total)))
genes_N.MGH105 <- as.data.frame(row.names(tpm(sceN.MGH105.total)))
genes_N.MGH104 <- as.data.frame(row.names(tpm(sceN.MGH104.total)))


ligand <- unique((Ramilowski["Ligand.ApprovedSymbol"]))
ligand_detected <- semi_join(genes_N.MGH104, ligand, by = c("row.names(tpm(sceN.MGH104.total))"="Ligand.ApprovedSymbol"))
colnames(ligand_detected) <- "ligand.dt"
ligand_dt <- as.character(ligand_detected$ligand.dt)

receptor <- unique((Ramilowski["Receptor.ApprovedSymbol"]))
receptor_detected <- semi_join(genes_N.MGH104, receptor, by = c("row.names(tpm(sceN.MGH104.total))"="Receptor.ApprovedSymbol"))
colnames(receptor_detected) <- "receptor.dt"
receptor_dt <- as.character(receptor_detected$receptor.dt)

sum_ligand <- sum(colSums(logcounts(sce.BT_S2.total[ligand_dt, ])))
sum_receptor <- sum(colSums(logcounts(sce.BT_S2.total[receptor_dt, ])))
