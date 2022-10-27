library(tidyverse)
library(plyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(msigdbr)


setwd('/home/hdd/alex/ncbi/public/sra/output/Trimmed/aligned')

#setwd('/home/hdd/alex/ncbi/public/sra/output/')

#load('./data/raw*')

#Pre-processing#########################################################################################
library(SingleCellExperiment)
library(scater)
library(scran)
library(uwot)
library(Rtsne)
library(PCAtools)
library("RColorBrewer")
library(d3heatmap)



countMatrix_aggregate <- read.table("countMatrix_aggregate", header = TRUE) #countMatrix_aggregate was made by collecting est_counts from each cell and transcripts with the same gene names were aggregated. None of the cells or genes were omited
data1 <- data.frame(countMatrix_aggregate[,-1], row.names=countMatrix_aggregate[,1])
meta_data1 <- read.csv("neo_metadata.csv")
sce <- SingleCellExperiment(assays = list(counts = as.matrix(data1)), colData = meta_data1 )



#QC matrix
mito <- grepl("^MT-", rownames(sce))
qc <- perCellQCMetrics(sce, subsets=list(Mito=mito))


#fixed criteria
criteria.lib <- qc$sum < 5e4
criteria.nexprs <- qc$detected < 3e3
criteria.mito <- qc$subsets_Mito_percent > 10
discard <- criteria.lib | criteria.nexprs | criteria.mito


#adaptive criteria
criteria.lib2 <- isOutlier(qc$sum, log=TRUE, type="lower")
criteria.nexprs2 <- isOutlier(qc$detected, log=TRUE, type="lower")
criteria.mito2 <- isOutlier(qc$subsets_Mito_percent, type="higher")
discard2 <- criteria.lib2 | criteria.nexprs2 | criteria.mito2


# Summarize the number of cells removed for each reason. https://osca.bioconductor.org/quality-control.html#identifying-low-quality-cells
DataFrame(LibSize=sum(criteria.lib2), NExprs=sum(criteria.nexprs2), MitoProp=sum(criteria.mito2), Total=sum(discard2))


#Filter out cells
filtered <- sce[,!discard2]
sce.filtered <- filtered[!mito, ]


#normalization
set.seed(100)
clust.zeisel <- quickCluster(sce.filtered) # This is needed to calculate size factors
#table(clust.zeisel)

#deconv.sf.zeisel <- calculateSumFactors(sce.final, cluster=clust.zeisel) # This is just to calculate a value called deconv.sf.zeisel, didn't apply on matrix
#summary(deconv.sf.zeisel) # To see it


sce.zeisel <- computeSumFactors(sce.filtered, cluster=clust.zeisel, min.mean=0.1) #Normalize the countMatrix
sce.zeisel <- logNormCounts(sce.zeisel)
assayNames(sce.zeisel)

dec.zeisel <- modelGeneVar(sce.zeisel)
fit.zeisel <- metadata(dec.zeisel)
plot(fit.zeisel$mean, fit.zeisel$var, xlab="Mean", ylab="Variance")
curve(fit.zeisel$trend(x), col="dodgerblue", add=TRUE, lwd=2) # visualizing the fit
dec.zeisel[order(dec.zeisel$bio, decreasing=TRUE),]


hvg <- getTopHVGs(dec.zeisel, var.threshold=0)
sce.hvg <- sce.zeisel[hvg,]



#FIND INTERESTING GENES#############################################################################


#1. CD133+SP+sphere-forming (cancer-stem-cell)
c2.sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP")
head(unique(c2.sets$gs_name))
cscMsigdbr.sets <- c2.sets[grep(pattern ="BEIER_GLIOMA_STEM_CELL_UP|HARRIS_BRAIN_CANCER_PROGENITORS|GUENTHER_GROWTH_SPHERICAL_VS_ADHERENT_UP", c2.sets$gs_name), ]#this gave the desired gene sets


#2. angiogenesis
c5.sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
head(unique(c5.sets$gs_name))
angioMsigdbr.sets <- c5.sets[grep(pattern ="GO_SPROUTING_ANGIOGENESIS", c5.sets$gs_name), ]


#3. drug resistance ABC group
ABC.sets <- read.csv("GeneCards_ABC.csv", header = TRUE) 


#4. hypoxia
hypo.sets <- read.csv("MGI_ResponseToHypoxia.csv", header = TRUE)
hypo.sets[1] <- apply(hypo.sets[1],2,toupper)

#5. csc geneCards
cscGeneCards.sets <- read.csv("GeneCards-GlioblastomaStemCell.csv", header = TRUE) 

#6. csc AmiGo
cscAmigo.sets <- read.delim("NeuralStemCell_AmiGo.csv", header = TRUE)[3]%>% group_by()  


#MATCHING GENES############################################################################################
csc_msigdbr <- semi_join(cscMsigdbr.sets, as.data.frame(hvg), by = c("human_gene_symbol"="hvg"))
csc_literature <- data.frame("PROM1", "SOX2", "")
angio_msigdbr <- semi_join(angioMsigdbr.sets, as.data.frame(hvg), by = c("human_gene_symbol"="hvg"))
ABC <- semi_join(ABC.sets, as.data.frame(hvg), by = c("Gene.Symbol"="hvg"))
hypoMGI <- semi_join(hypo.sets, as.data.frame(hvg), by = c("Symbol"="hvg"))
cscGeneCards <- semi_join(cscGeneCards.sets, as.data.frame(hvg) , by = c("Gene.Symbol"="hvg"))


#PCA ON SELECTED GENES###################################################################################################


csc <- csc_msigdbr[["human_gene_symbol"]]
sce.csc <- sce.zeisel[csc,]
set.seed(100) # See below.
sce.csc <- runPCA(sce.csc, subset_row=csc) 
reducedDim(sce.csc, "PCA_25") <- reducedDim(sce.csc, "PCA")[,1:25]
set.seed(00101001101)
sce.csc <- runTSNE(sce.csc, dimred="PCA")
plotReducedDim(sce.csc, dimred="TSNE", colour_by = "PTPRZ1")


cscGC <- as.character(cscGeneCards[["Gene.Symbol"]])
sce.cscGC <- sce.hvg[cscGC,]
set.seed(100) # See below.
sce.cscGC <- runPCA(sce.cscGC, subset_row=cscGC) 
reducedDim(sce.cscGC, "PCA_25") <- reducedDim(sce.cscGC, "PCA")[,1:25]
set.seed(00101001101)
sce.cscGC <- runTSNE(sce.cscGC, dimred="PCA")
plotReducedDim(sce.cscGC, dimred="TSNE", colour_by = "VIM")


pca <- PCAtools::pca(counts(sce.cscGC))
pca_result <- PCAtools::getLoadings(pca)
sce.csc <- runTSNE(sce.zeisel, dimred="PCA")


set.seed(100) # See below.
sce.zeisel <- runPCA(sce.zeisel, subset_row=csc)
reducedDimNames(sce.zeisel)
ncol(reducedDim(sce.zeisel))
reducedDim(sce.zeisel, "PCA_25") <- reducedDim(sce.zeisel, "PCA")[,1:25]


#PCAtools
pca <- PCAtools::pca(counts(sce.zeisel))
pca_result <- PCAtools::getLoadings(pca)
PCAtools::getComponents(pca) #get pc no.
#tsne
set.seed(00101001101)
sce.zeisel <- runTSNE(sce.zeisel, feature_set = csc[5])
plotReducedDim(sce.zeisel, dimred="TSNE", colour_by = "AC010970.1")


#CORRELATION MATRIX#####################################################################################################

#combine csc genes and angiogenesis genes
cscGC <- as.character(cscGeneCards[["Gene.Symbol"]])
angioM <- angio_msigdbr[["human_gene_symbol"]]
cscGC_angioM <- c(cscGC, angioM)
  
cscGC_angioM.counts <- counts(sce.zeisel[cscGC_angioM,])

tot <- cor(t(counts(sce.hvg)))
res <- cor(t(cscGC_angioM.counts))>0.2 
cscGC_to_angioM <- as.data.frame(res[c(cscGC), c(angioM)])  
#cormat <- round(cor(csc.counts),2)
#round(res, 2)
d3heatmap(cscGC_to_angioM, scale = "column")

CandA <- data.frame(cor_score = matrix(unlist(cscGC_to_angioM), nrow=197830))
hvgCor <- data.frame(cor_score = matrix(unlist(tot)))
ggplot(CandA, aes(cor_score)) + geom_histogram(bins = 100) + scale_y_log10()
ggplot(hvgCor, aes(cor_score)) + geom_histogram(bins = 100) + scale_y_log10()


#############
cscM <- as.character(csc_msigdbr[["human_gene_symbol"]])
angioM <- angio_msigdbr[["human_gene_symbol"]]
cscM_angioM <- c(cscM, angioM)
cscM_angioM.counts <- counts(sce.zeisel[cscM_angioM,])
mean <- rowMeans(cscM_angioM.counts)
res <- cor(t(cscM_angioM.counts) )
cscM_to_angioM <- res[c(cscM), c(angioM)]
#cormat <- round(cor(csc.counts),2)
#round(res, 2)
d3heatmap(cscM_to_angioM, scale = "column")

#dCor#
library("energy")
cscM.counts <- counts(sce.zeisel[cscM,])
angioM.counts <- counts(sce.zeisel[angioM,])

correlation <- dcor(t(cscM.counts), t(angioM.counts))

###########
hypoMG <- as.character(hypoMGI[["Symbol"]]) 
angioM_hypoMG <- c(angioM, hypoMG)

angioM_hypoMG.counts <- counts(sce.zeisel[angioM_hypoMG,])


res <- cor(t(angioM_hypoMG.counts) )
angioM_to_hypoMG <- res[c(angioM), c(hypoMG)]
#cormat <- round(cor(csc.counts),2)
#round(res, 2)
d3heatmap(angioM_to_hypoMG, scale = "column", colors = "Blues")

############################################################################################################
#BIOCONDUCTOR###########################################################################################


#dimension reduction
set.seed(100) # See below.
sce.zeisel <- runPCA(sce.zeisel, subset_row=hvg) 
reducedDimNames(sce.zeisel)
ncol(reducedDim(sce.zeisel))


#choosing amount of PCs
percent.var <- attr(reducedDim(sce.zeisel), "percentVar")
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
chosen.elbow
par(mar = c(5,5,3,1))
plot(percent.var, xlab="PC", ylab="Variance explained (%)")
abline(v=chosen.elbow, col="red")
#Choose top 25 PCs
reducedDim(sce.zeisel, "PCA_25") <- reducedDim(sce.zeisel, "PCA")[,1:25]
reducedDimNames(sce.zeisel)
plotReducedDim(sce.zeisel, dimred="PCA", colour_by = "patient_id")


pca <- PCAtools::pca(counts(sce.hvg))
pca_result <- PCAtools::getLoadings(pca)

set.seed(00101001101)
sce.zeisel <- runTSNE(sce.zeisel, dimred="PCA")
plotReducedDim(sce.zeisel, dimred="TSNE", colour_by = "AC010970.1")

#######################################################################################

