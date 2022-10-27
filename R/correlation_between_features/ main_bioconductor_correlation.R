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



#ALL_GENE LISTS FROM EACH PATIENT###############################################################
genes.BT_S2 <- row.names(counts(sce.BT_S2.total))
genes.BT_S1 <- row.names(counts(sce.BT_S1.total))
genes.BT_S4 <- row.names(counts(sce.BT_S4.total))
genes.BT_S6 <- row.names(counts(sce.BT_S6.total))




#DCOR AND BOOTSTRAPPING#########################################################################
#INTERESTED GENE SETS
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
#DCOR
source("/home/hdd/yue/code/R/Darmanis/dcor.R")
cscM_diffA.dcor2 <- distance.correlation(cscM_2, diffA_2, sce.BT_S2.total)
cscM_migrHM.dcor2 <- distance.correlation(cscM_2, migrHM_2, sce.BT_S2.hvg)
cscM_angioA.dcor2 <- distance.correlation(cscM_2, angioA_2, sce.BT_S2.total)
cscM_srHM.dcor2 <- distance.correlation(cscM_2, srHM_2, sce.BT_S2.total)

cscM_diffA.dcor1 <- distance.correlation(cscM_1, diffA_1, sce.BT_S1.total)
cscM_migrHM.dcor1 <- distance.correlation(cscM_1, migrHM_1, sce.BT_S1.total)
cscM_angioA.dcor1 <- distance.correlation(cscM_1, angioA_1, sce.BT_S1.total)
cscM_srHM.dcor1 <- distance.correlation(cscM_1, srHM_1, sce.BT_S1.total)

cscM_diffA.dcor4 <- distance.correlation(cscM_4, diffA_4, sce.BT_S4.total)
cscM_migrHM.dcor4 <- distance.correlation(cscM_4, migrHM_4, sce.BT_S4.total)
cscM_angioA.dcor4 <- distance.correlation(cscM_4, angioA_4, sce.BT_S4.total)
cscM_srHM.dcor4 <- distance.correlation(cscM_4, srHM_4, sce.BT_S4.total)

cscM_diffA.dcor6 <- distance.correlation(cscM_6, diffA_6, sce.BT_S6.total)
cscM_migrHM.dcor6 <- distance.correlation(cscM_6, migrHM_6, sce.BT_S6.total)
cscM_angioA.dcor6 <- distance.correlation(cscM_6, angioA_6, sce.BT_S6.total)
cscM_srHM.dcor6 <- distance.correlation(cscM_6, srHM_6, sce.BT_S6.total)

proteo <- vector(read_tsv("/home/hdd/yue/data/text/Proteoglycans in cancer.tsv"))
proteo <- c("CD44", "SRC", "CTTN", "HCLS1", "ERBB2", "GRB2", "VAV3", "VAV1", "HBEGF", "HRAS", "KRAS", "NRAS", "RRAS", "RRAS2", "MRAS", "BRAF", "RAF1", "ARAF", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "RAC1", "IQGAP1", "CDC42", "ELK1", "ESR1", "CCND1", "ACTG1", "ACTB", "FLNA", "FLNC", "FLNB", "PAK1", "TIAM1", "ARHGEF1", "RHOA", "ROCK1", "ROCK2", "ANK1", "ANK2", "ANK3", "GAB1", "PIK3CA", "PIK3CD", "PIK3CB", "PIK3R1", "PIK3R2", "PIK3R3", "AKT1", "AKT2", "AKT3", "SLC9A1", "PPP1CA", "PPP1CB", "PPP1CC", "PPP1R12A", 
            "PPP1R12B", "PPP1R12C", "ARHGEF12", "PLCE1", "ITPR1", "ITPR2", "ITPR3", "CAMK2A", "CAMK2D", "CAMK2B", "CAMK2G", "NANOG", "NANOGP8", "DDX5", "DROSHA", "STAT3", "MIR21", "TWIST1", "TWIST2", "MIR10A", "MIR10B", "HOXD10", "DCN", "IGF1", "IGF1R", "MTOR", "PDPK1", "RPS6KB1", "RPS6KB2", "EIF4B", "RPS6", "EGFR", "CAV1", "CAV2", "CAV3", "CD63", "CDKN1A", "CASP3", "TGFB1", "TLR2", "TLR4", "PDCD4", "TNF", "IL12B", "ERBB3", "ERBB4", "MYC", "CTNNB1", "HIF1A")
proteo.common <- intersect(genes.BT_S2, proteo)
cscM_proteo.dcor6 <- distance.correlation(cscM_6, proteo.common, sce.BT_S6.total)

#BOOTSTRAPPING
source("/home/hdd/yue/code/R/Darmanis/dcor_bootstrapping.R")

cscM_rand.dcor6 <- dcor.bootstrapping(cscM_6, proteo.common, genes.BT_S6, sce.BT_S6.total)
#(cscM_diffA.dcor2-mean(cscM_rand.dcor2))/sd(cscM_rand.dcor2)
ggplot(data.frame(dcor=cscM_rand.dcor6),
       aes(x=dcor)) +
  geom_histogram(bins = 50)+
  scale_x_log10() +
  geom_vline(xintercept=cscM_proteo.dcor6,
               color="red", linetype="dashed", size=1) +
  geom_vline(xintercept=mean(cscM_rand.dcor6), 
               color="blue", linetype="dashed", size=1) +
 # geom_segment(aes(x = mean, y = 15, xend = cscM_diffA.dcor, yend = 15, color = "grey3"))   



library(ggsignif)
test <- rbind(cscM_rand.dcor, cscM_diffA.dcor)


#KNN dcor#######################################################################################
source("/home/hdd/yue/code/R/Darmanis/KNN_dcor.R")
source("/home/hdd/yue/code/R/Darmanis/dcor.R")

sce.BT_S2.total$cscM_angioA_dcor <- runif(ncol(sce.BT_S2.total))
sce.BT_S2.total$cscM_angioA_dcor <- kNN.dcor(cscM_2, angioA_2, sce.BT_S2.total)

sce.BT_S2.total$cscM_rand_dcor <- runif(ncol(sce.BT_S2.total))
sce.BT_S2.total$cscM_rand_dcor <- kNN.dcor(cscM_2, rand, sce.BT_S2.total)

plotReducedDim(sce.BT_S2.total, dimred="TSNE", colour_by = "cscM_rand_dcor" )
plotReducedDim(sce.BT_S2.total, dimred="TSNE", colour_by = "NOTCH1" )


source("/home/hdd/yue/code/R/Darmanis/KNN_dcor_sigma.R")
test <- kNN.dcor.sigma(cscM_2, angioA_2, genes.BT_S2, sce.BT_S2.total)
csc_diffA_sigma <- test
sce.BT_S2.total$cscM_diffA_dcor_sigma <- runif(ncol(sce.BT_S2.total))
sce.BT_S2.total$cscM_diffA_dcor_sigma <- csc_diffA_sigma
plotReducedDim(sce.BT_S2.total, dimred="TSNE", colour_by = "cscM_diffA_dcor_sigma" )

csc_angioA_sigma <- kNN.dcor.sigma(cscM, angioA, genes.BT_S2, sce.BT_S2.total)
sce.BT_S2.total$cscM_angioA_dcor_sigma <- runif(ncol(sce.BT_S2.total))
sce.BT_S2.total$cscM_angioA_dcor_sigma <- csc_angioA_sigma
plotReducedDim(sce.BT_S2.total, dimred="TSNE", colour_by = "cscM_angioA_dcor_sigma" )

csc_srHM_sigma <- kNN.dcor.sigma(cscM, srHM, genes.BT_S2, sce.BT_S2.total)
sce.BT_S2.total$cscM_srHM_dcor_sigma <- runif(ncol(sce.BT_S2.total))
sce.BT_S2.total$cscM_srHM_dcor_sigma <- csc_srHM_sigma
plotReducedDim(sce.BT_S2.total, dimred="TSNE", colour_by = "cscM_srHM_dcor_sigma" )

csc_hypoA_sigma <- kNN.dcor.sigma(cscM, hypoA, genes.BT_S2, sce.BT_S2.total)
sce.BT_S2.total$cscM_hypoA_dcor_sigma <- runif(ncol(sce.BT_S2.total))
sce.BT_S2.total$cscM_hypoA_dcor_sigma <- csc_hypoA_sigma
plotReducedDim(sce.BT_S2.total, dimred="TSNE", colour_by = "cscM_hypoA_dcor_sigma" )

csc_G1G0A_sigma <- kNN.dcor.sigma(cscM, G1G0A, genes.BT_S2, sce.BT_S2.total)
sce.BT_S2.total$cscM_G1G0A_dcor_sigma <- runif(ncol(sce.BT_S2.total))
sce.BT_S2.total$cscM_G1G0A_dcor_sigma <- csc_G1G0A_sigma
plotReducedDim(sce.BT_S2.total, dimred="TSNE", colour_by = "cscM_G1G0A_dcor_sigma" )

csc_migrHM_sigma <- kNN.dcor.sigma(cscM, migrHM, genes.BT_S2, sce.BT_S2.total)
sce.BT_S2.total$cscM_migrHM_dcor_sigma <- runif(ncol(sce.BT_S2.total))
sce.BT_S2.total$cscM_migrHM_dcor_sigma <- csc_migrHM_sigma
plotReducedDim(sce.BT_S2.total, dimred="TSNE", colour_by = "cscM_migrHM_dcor_sigma" )

sigma_diffA_angioA <- data.frame(diffA = csc_diffA_sigma, angioA = csc_angioA_sigma)

ggplot(sigma_diffA_angioA, aes(x = csc_diffA_sigma, y = csc_angioA_sigma)) +
  geom_point()




#SUM ON PHEATMAP
source("/home/hdd/yue/code/R/Darmanis/sum_sigma.R")

set.seed(100)
rand <- sample(nrow(sce.BT_S2.total), 78)
sumall= colSums(logcounts(sce.BT_S2.total))
sum.all <- as.matrix(data.frame(csc=colSums(z.score[cscM_2, ]),#/sumall, 
                                 diff=colSums(z.score[diffA_2, ]),#/sumall,
                                 migr=colSums(z.score[migrHM_2, ]),#/sumall, 
                                 angio=colSums(z.score[angioA_2, ]),
                                 random =colSums(z.score[rand, ])
                                ))

z.score <- (cpm(sce.BT_S2.total)-colMeans(cpm(sce.BT_S2.total)))/colSds(cpm(sce.BT_S2.total))
#/sumall random =colSums(cpm(sce.BT_S2.hvg[rand, ]))
sum.all <- as.matrix(data.frame(csc=colSums(cpm(sce.BT_S6.total[cscM_6, ])),#/sumall, 
                                diff=colSums(cpm(sce.BT_S6.total[diffA_6, ])),#/sumall,
                                migr=colSums(cpm(sce.BT_S6.total[migrHM_6, ])),#/sumall, 
                                angio=colSums(cpm(sce.BT_S6.total[angioA_6, ]))
))

sumsigma <- as.matrix(data.frame(csc=sum.sigma(cscM_6, sce.BT_S6.total),#/sumall, 
                                  diff=sum.sigma(diffA_6, sce.BT_S6.total),#/sumall,
                                  migr=sum.sigma(migrHM_6, sce.BT_S6.total),#/sumall, 
                                  angio=sum.sigma(angioA_6, sce.BT_S6.total),#/sumall
                                  random =sum.sigma(rand, sce.BT_S6.total)))
library(pheatmap)
pm <- pheatmap(t(sum.all), scale="row", show_colnames = FALSE,
               cellwidth=1.5, cellheight =60)



#CORRELATION MATRIX################################################################################
source("/home/hdd/yue/code/R/Darmanis/correlation_matrix.R")
cscM_diffA.cor <- correlation.matrix(cscM_2, diffA_2, sce.BT_S2.total)
rank_cscM_diffA <- rownames(cscM_diffA.cor[order(rowSums(cscM_diffA.cor),decreasing=T),][1:30,])
cscM_angioA.cor <- correlation.matrix(cscM_2, angioA_2, sce.BT_S2.total)
rank_cscM_angioA <- rownames(cscM_angioA.cor[order(rowSums(cscM_angioA.cor),decreasing=T),][1:30,])
cscM_migrHM.cor <- correlation.matrix(cscM_2, migrHM_2, sce.BT_S2.total) 
rank_cscM_migrHM <- rownames(cscM_migrHM.cor[order(rowSums(cscM_migrHM.cor),decreasing=T),][1:30,])
common2 <- intersect(common, rank_cscM_migrHM)
d3heatmap(cscM_migrHM.cor, scale = "column", dendrogram = "none", color = "Blues")

#HIERACHICAL CLUSTERING FOR EVERY GENE#################################################################

source("/home/hdd/yue/code/R/Darmanis/hvg_standard.R")

library(multiClust)
library(som)
CscDiffAngioMigr2 <- c(cscM_2, diffA_2, angioA_2, migrHM_2)
features <- data.frame(csc = cscM_2, diff = diffA_2, angio = angioA_2, migr = migrHM_2)
CscDiffAngio.hvg4 <- c(common.cscM4, common.diffA4, common.angioA4 )
sce.BT_S2.features2 <- cpm(sce.BT_S2.total[CscDiffAngioMigr2, ])
write.table(sce.BT_S2.features2, file="/home/hdd/yue/data/output/sce.BT_S2.features2")
write.table(cscM_2, file="/home/hdd/yue/data/output/cscM_2")
write.table(diffA_2, file="/home/hdd/yue/data/output/diffA_2")
write.table(angioA_2, file="/home/hdd/yue/data/output/angioA_2")
write.table(migrHM_2, file="/home/hdd/yue/data/output/migrHM_2")




sce.BT_S2.features1 <- logcounts(sce.BT_S2.features)[, 1:2]

sce.BT_S2.rank <- colRanks(sce.BT_S2.features1)

sce.BT_S2.hvg <- hvg.standard(sce.BT_S2.total)
sce.BT_S1.hvg <- hvg.standard(sce.BT_S1.total)
sce.BT_S4.hvg <- hvg.standard(sce.BT_S4.total)
sce.BT_S6.hvg <- hvg.standard(sce.BT_S6.total)

genes.BT_S2.hvg <- row.names(counts(sce.BT_S2.hvg))
genes.BT_S1.hvg <- row.names(counts(sce.BT_S1.hvg))
genes.BT_S4.hvg <- row.names(counts(sce.BT_S4.hvg))
genes.BT_S6.hvg <- row.names(counts(sce.BT_S6.hvg))


common.cscM2 <- intersect(genes.BT_S2.hvg, cscM_2)
common.migrHM2 <- intersect(genes.BT_S2.hvg, migrHM_2)
common.angioA2 <- intersect(genes.BT_S2.hvg, angioA_2)
common.diffA2 <- intersect(genes.BT_S2.hvg, diffA_2)

common.cscM1 <- intersect(genes.BT_S1.hvg, cscM_1)
common.migrHM1 <- intersect(genes.BT_S1.hvg, migrHM_1)
common.angioA1 <- intersect(genes.BT_S1.hvg, angioA_1)
common.diffA1 <- intersect(genes.BT_S1.hvg, diffA_1)

common.cscM4 <- intersect(genes.BT_S4.hvg, cscM_4)
common.migrHM4 <- intersect(genes.BT_S4.hvg, migrHM_4)
common.angioA4 <- intersect(genes.BT_S4.hvg, angioA_4)
common.diffA4 <- intersect(genes.BT_S4.hvg, diffA_4)

common.cscM6 <- intersect(genes.BT_S6.hvg, cscM_6)
common.migrHM6 <- intersect(genes.BT_S6.hvg, migrHM_6)
common.angioA6 <- intersect(genes.BT_S6.hvg, angioA_6)
common.diffA6 <- intersect(genes.BT_S6.hvg, diffA_6)

source("/home/hdd/yue/code/R/Darmanis/dcor.R")
cscM_diffA.dcor2 <- distance.correlation(common.cscM2, common.diffA2, sce.BT_S2.hvg)
cscM_migrHM.dcor2 <- distance.correlation(common.cscM2, common.migrHM2, sce.BT_S2.hvg)
cscM_angioA.dcor2 <- distance.correlation(common.cscM2, common.angioA2, sce.BT_S2.hvg)
cscM_srHM.dcor2 <- distance.correlation(common.cscM2, common.srHM2, sce.BT_S2.hvg)

cscM_diffA.dcor1 <- distance.correlation(common.cscM1, common.diffA1, sce.BT_S1.hvg)
cscM_migrHM.dcor1 <- distance.correlation(common.cscM1, common.migrHM1, sce.BT_S1.hvg)
cscM_angioA.dcor1 <- distance.correlation(common.cscM1, common.angioA1, sce.BT_S1.hvg)
cscM_srHM.dcor1 <- distance.correlation(common.cscM1, common.srHM1, sce.BT_S1.hvg)

cscM_diffA.dcor4 <- distance.correlation(common.cscM4, common.diffA4, sce.BT_S4.hvg)
cscM_migrHM.dcor4 <- distance.correlation(common.cscM4, common.migrHM4, sce.BT_S4.hvg)
cscM_angioA.dcor4 <- distance.correlation(common.cscM4, common.angioA4, sce.BT_S4.hvg)
cscM_srHM.dcor4 <- distance.correlation(common.cscM4, common.srHM4, sce.BT_S4.hvg)

cscM_diffA.dcor6 <- distance.correlation(common.cscM6, common.diffA6, sce.BT_S6.hvg)
cscM_migrHM.dcor6 <- distance.correlation(common.cscM6, common.migrHM6,sce.BT_S6.hvg)
cscM_angioA.dcor6 <- distance.correlation(common.cscM6, common.angioA6, sce.BT_S6.hvg)
cscM_srHM.dcor6 <- distance.correlation(common.cscM6, common.srHM6, sce.BT_S6.hvg)

#BOOTSTRAPPING
source("/home/hdd/yue/code/R/Darmanis/dcor_bootstrapping.R")

cscM_rand.dcor6 <- dcor.bootstrapping(common.cscM6, common.migrHM6, genes.BT_S6.hvg, sce.BT_S6.hvg)
#(cscM_diffA.dcor2-mean(cscM_rand.dcor2))/sd(cscM_rand.dcor2)
ggplot(data.frame(dcor=cscM_rand.dcor2),
       aes(x=dcor)) +
  geom_histogram(bins = 50)+
  scale_x_log10() +
  geom_vline(xintercept=cscM_diffA.dcor2,
             color="red", linetype="dashed", size=1) +
  geom_vline(xintercept=mean(cscM_rand.dcor2), 
             color="blue", linetype="dashed", size=1)

sce.BT_S2.diffA <- sce.BT_S2.total[diffA_2, ]
sce.BT_S2.zscore <- scale(logcounts(sce.BT_S2.features))
hclust_analysis <- cluster_analysis(sel.exp=sce.BT_S2.zscore,
                                    cluster_type="HClust",
                                    distance="euclidean", linkage_type="ward.D2", 
                                    gene_distance="pearson",
                                    num_clusters=3, data_name="sce.BT_S2.total", 
                                    probe_rank="SD_Rank", probe_num_selection="Fixed_Probe_Num",
                                    cluster_num_selection="Fixed_Clust_Num")

pm <- pheatmap(t(logcounts(sce.BT_S2.features4)), scale="column", show_colnames = FALSE,
               show_rownames = FALSE, cluster_rows = TRUE, cluster_cols = FALSE,
               cellwidth=6, cellheight =6, clustering_distance_cols="correlation")
#               color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100), 
proteo_csc <- c(cscM_2, proteo.common)
sce.BT_S2.pc <- sce.BT_S2.total[proteo_csc, ]
d3heatmap(logcounts(sce.BT_S2.pc), scale = "row", dendrogram = "column", color = "RdBu")







