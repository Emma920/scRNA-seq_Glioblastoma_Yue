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




#LOAD FILES####################################################################################
#source("/home/hdd/yue/code/R/Darmanis/load_files.R")
#countMatrix <- knit_aggregate("/home/hdd/yue/data/aligned/Darmanis/kallisto/", "*.tsv")
#write.table(countMatrix, file="/home/hdd/yue/data/CM/Darmanis/countMatrix")
source("/home/hdd/yue/code/R/Darmanis/load_STARfiles.R")
countMatrix_S <- knit_STARaggregate("/home/hdd/yue/data/aligned/Patel/STAR/", "*ReadsPerGene.out.tab")
write.table(countMatrix_S, file="/home/hdd/yue/data/CM/Patel/countMatrix_S")

#LOAD DATA#####################################################################################
source("/home/hdd/yue/code/R/Darmanis/load_data.R")
source("/home/hdd/yue/code/R/Darmanis/pre_processing.R")
source("/home/hdd/yue/code/R/Darmanis/pre_processing_lib.R")


sce.Patel <-load.data("/home/hdd/yue/data/CM/Patel/countMatrix_S", "/home/hdd/yue/data/text/ncbi_Acc_list/Patel/Patel_meta.csv" )
sce.Patel.total <- pre.processing(sce.Patel)



#INTERESTED GENES##############################################################################
source("/home/hdd/yue/code/R/Darmanis/interested_genes.R")
genes <- row.names(counts(sce.total))



#DIMENSIONAL REDUCTION ON SELECTED GENE SETS###################################################
source("/home/hdd/yue/code/R/Darmanis/dimension_reduction.R")
sce.total <- dimension.reduction(genes, sce.total)
plotReducedDim(sce.total, dimred="TSNE", colour_by = "")                 
#saveRDS(sce.total, file = "/home/hdd/yue/data/rds/Darmanis/sce.total.rds")

################################################################################################
#All cells######################################################################################
################################################################################################


#INTERESTED GENE SETS###########################################################################
cscM_total <- interested.genes("cscM", sce.total)
#csc.GC_M <- unique(c(interested.genes("cscGC", sce.hvg), interested.genes("cscM", sce.hvg)))
diffA_total <- interested.genes("diffA", sce.total)
angioA_total <- interested.genes("angioA", sce.total)
angioM_total <- interested.genes("angioM", sce.total)
cycleA_total <- interested.genes("cycleA", sce.total)[1:100,]
hypoA_total <- interested.genes("hypoA", sce.total)
migrHM_total <- interested.genes("migrHM", sce.total)
proHM_total <- interested.genes("proHM", sce.total)
#nscA_total <- interested.genes("nscA", sce.total)
#nscHM_total <- interested.genes("nscHM", sce.total)


#SUM ON TSNE#####################################################################################
cscM.sum <- colSums(logcounts(sce.total[interested.genes("cscM_total", sce.total), ]))
sce.total$cscM_sum <- runif(ncol(sce.total))
sce.total$cscM_sum <- cscM.sum
rm(cscM.sum)

plotReducedDim(sce.total, dimred="TSNE", colour_by = "cscM_sum" )

rand <- sample(nrow(sce.total), 1000)
sce.rand <- sce.total[rand,]
rand.sum <- colSums(logcounts(sce.rand))
sce.total$rand_sum <- runif(ncol(sce.total))
sce.total$rand_sum <- rand.sum
plotReducedDim(sce.total, dimred="TSNE", colour_by = "rand_sum" )


###############################################################################################
#WITH EACH PATIENT#############################################################################
###############################################################################################

#SUBSET SCE WITH PATIENTS ID###################################################################

sce.BT_S2.total <- pre.processing.lib(subset(sce, , patient_id=="BT_S2"))
sce.BT_S1.total <- pre.processing.lib(subset(sce, , patient_id=="BT_S1"))
sce.BT_S4.total <- pre.processing.lib(subset(sce, , patient_id=="BT_S4"))
sce.BT_S6.total <- pre.processing.lib(subset(sce, , patient_id=="BT_S6"))


#ALL_GENE LISTS FROM EACH PATIENT###############################################################
genes.BT_S2 <- row.names(counts(sce.BT_S2.total))
genes.BT_S1 <- row.names(counts(sce.BT_S1.total))
genes.BT_S4 <- row.names(counts(sce.BT_S4.total))
genes.BT_S6 <- row.names(counts(sce.BT_S6.total))


#DIMENSION REDUCTION ON EACH PATIENT CELLS######################################################
sce.BT_S2.total <- dimension.reduction(genes.BT_S2, sce.BT_S2.total)
sce.BT_S1.total <- dimension.reduction(genes.BT_S1, sce.BT_S1.total)
sce.BT_S4.total <- dimension.reduction(genes.BT_S4, sce.BT_S4.total)
sce.BT_S6.total <- dimension.reduction(genes.BT_S6, sce.BT_S6.total)
#plotReducedDim(sce.BT_S2.total, dimred="TSNE", colour_by = "plate_id")   


#DCOR AND BOOTSTRAPPING#########################################################################
#INTERESTED GENE SETS
cscM <- interested.genes("cscM", sce.BT_S2.total)
diffA <- interested.genes("diffA", sce.BT_S2.total)
angioA <- interested.genes("angioA", sce.BT_S2.total)
angioM <- interested.genes("angioM", sce.BT_S2.total)
cycleA <- interested.genes("cycleA", sce.BT_S2.total)
hypoA_0 <- interested.genes("hypoA", sce.BT_S2.total)
migrHM <- interested.genes("migrHM", sce.BT_S2.total)
proHM <- interested.genes("proHM", sce.BT_S2.total)
#nscA <- interested.genes("angioM", sce.total)
#nscHM <- interested.genes("angioM", sce.total)
hypoA <- hypoA_0[!(hypoA_0 %in% angioA)]
srHM <- interested.genes("srHM", sce.BT_S2.total)
G1G0A <- c("SMPD3","CYP27B1","RAB11FIP4","EZH2","RIDA","PHGDH","CAPN3","ZNF503")

cscM_1 <- interested.genes("cscM", sce.BT_S1.total)
diffA_1 <- interested.genes("diffA", sce.BT_S1.total)
angioA_1 <- interested.genes("angioA", sce.BT_S1.total)
migrHM_1 <- interested.genes("migrHM", sce.BT_S1.total)

cscM_4 <- interested.genes("cscM", sce.BT_S4.total)
diffA_4 <- interested.genes("diffA", sce.BT_S4.total)
angioA_4 <- interested.genes("angioA", sce.BT_S4.total)
migrHM_4 <- interested.genes("migrHM", sce.BT_S4.total)

cscM_6 <- interested.genes("cscM", sce.BT_S6.total)
diffA_6 <- interested.genes("diffA", sce.BT_S6.total)
angioA_6 <- interested.genes("angioA", sce.BT_S6.total)
migrHM_6 <- interested.genes("migrHM", sce.BT_S6.total)
#DCOR
source("/home/hdd/yue/code/R/Darmanis/dcor.R")
cscM_diffA.dcor <- distance.correlation(cscM_6, diffA_6, sce.BT_S6.total)
cscM_migrHM.dcor <- distance.correlation(cscM_6, migrHM_6, sce.BT_S6.total)
cscM_angioA.dcor <- distance.correlation(cscM_6, angioA_6, sce.BT_S6.total)
cscM_srHM.dcor <- distance.correlation(cscM, srHM, sce.BT_S1.total)


#BOOTSTRAPPING
source("/home/hdd/yue/code/R/Darmanis/dcor_bootstrapping.R")

cscM_rand.dcor <- data.frame(dcor = dcor.bootstrapping(cscM_6, migrHM_6, genes.BT_S6, sce.BT_S6.total))
ggplot(cscM_rand.dcor,
       aes(x=dcor)) +
  geom_histogram(bins = 50)+
  scale_x_log10() +
  geom_vline(xintercept=cscM_migrHM.dcor,
             color="red", linetype="dashed", size=1) +
  geom_vline(xintercept=mean(cscM_rand.dcor[[1]]), 
             color="blue", linetype="dashed", size=1) +
  # geom_segment(aes(x = mean, y = 15, xend = cscM_diffA.dcor, yend = 15, color = "grey3"))   
  
  
  
  library(ggsignif)
test <- rbind(cscM_rand.dcor, cscM_diffA.dcor)


#KNN dcor#######################################################################################
source("/home/hdd/yue/code/R/Darmanis/KNN_dcor.R")
sce.BT_S2.total$cscM_srHM_dcor <- runif(ncol(sce.BT_S2.total))
sce.BT_S2.total$cscM_srHM_dcor <- kNN.dcor(cscM, srHM, sce.BT_S2.total)

plotReducedDim(sce.BT_S2.total, dimred="TSNE", colour_by = "cscM_diffA_dcor" )

source("/home/hdd/yue/code/R/Darmanis/KNN_dcor_sigma.R")
test <- kNN.dcor.sigma(cscM, diffA, genes.BT_S2, sce.BT_S2.total)
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

#SUM SIGMA#####################################################################################
source("/home/hdd/yue/code/R/Darmanis/sum_sigma.R")
cscM_sigma <- sum.sigma(cscM, sce.BT_S2.total)
sce.BT_S2.total$cscM_sigma <- runif(ncol(sce.BT_S2.total))
sce.BT_S2.total$cscM_sigma <- cscM_sigma
plotReducedDim(sce.BT_S2.total, dimred="TSNE", colour_by = "cscM_sigma" )

cscM_sigma <- sum.sigma(cscM, sce.BT_S2.total)
sce.BT_S2.total$cscM_sigma <- runif(ncol(sce.BT_S2.total))
sce.BT_S2.total$cscM_sigma <- cscM_sigma
plotReducedDim(sce.BT_S2.total, dimred="TSNE", colour_by = "cscM_sigma" )

sum.cscM <- as.matrix(data.frame(csc=colSums(logcounts(sce.BT_S1.total[cscM_1, ])), diff=colSums(logcounts(sce.BT_S6.total[diffA_6, ]))
                                 ,migr=colSums(logcounts(sce.BT_S6.total[migrHM_6, ])), angio=colSums(logcounts(sce.BT_S6.total[angioA_6, ]))))

heatmap.2(sum.cscM, scale="column")
pm <- pheatmap(t(sum.cscM), scale="row", show_colnames = FALSE,
               color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100), cellwidth=1.5, cellheight =60)

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

#CORRELATION MATRIX################################################################################
source("/home/hdd/yue/code/R/Darmanis/correlation_matrix.R")
cscM_diffA.cor <- correlation.matrix(cscM, diffA, sce.BT_S2.total)
rank_cscM_diffA <- rownames(cscM_diffA.cor[order(rowSums(cscM_diffA.cor),decreasing=T),][1:30,])
cscM_angioA.cor <- correlation.matrix(cscM, angioA, sce.BT_S2.total)
rank_cscM_angioA <- rownames(cscM_angioA.cor[order(rowSums(cscM_angioA.cor),decreasing=T),][1:30,])
cscM_migrHM.cor <- correlation.matrix(cscM, migrHM, sce.BT_S2.total) 
rank_cscM_migrHM <- rownames(cscM_migrHM.cor[order(rowSums(cscM_migrHM.cor),decreasing=T),][1:30,])
common2 <- intersect(common, rank_cscM_migrHM)
d3heatmap(cscM_diffA.cor, scale = "column", dendrogram = "none", color = "Blues")


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


#FIND MARKER GENES FOR CLUSTERS#################################################################
g <- buildSNNGraph(sce.BT_S2.total, k=20, use.dimred = 'PCA') #using the top PCs
clust.louvain <- igraph::cluster_louvain(g)$membership
sce.BT_S2.total$label <- factor(clust.louvain)
markers.BT_S2 <- findMarkers(sce.BT_S2.total, groups=sce.BT_S2.total$label)
rm(g, clust.louvain)
plotReducedDim(sce.BT_S2.total, dimred="TSNE", colour_by = "label")


