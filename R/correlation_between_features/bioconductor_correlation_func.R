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


load.data <- function(cm.original) {
  countMatrix_aggregate <- read.table(cm.original, header = TRUE) #countMatrix_aggregate was made by collecting est_counts from each cell and transcripts with the same gene names were aggregated. None of the cells or genes were omited
  data1 <- data.frame(countMatrix_aggregate[,-1], row.names=countMatrix_aggregate[,1])
  meta_data1 <- read.csv("neo_metadata.csv")
  sce <- SingleCellExperiment(assays = list(counts = as.matrix(data1)), colData = meta_data1)
  return(sce)
}

sce <-load.data("countMatrix_aggregate")

pre.processing <- function(sce) {
  
  # QC MATRIX
  mito <- grepl("^MT-", rownames(sce))
  qc <- perCellQCMetrics(sce, subsets=list(Mito=mito))
  
  
  # FIXED CRITERIA
  criteria.lib <- qc$sum < 5e4
  criteria.nexprs <- qc$detected < 3e3
  criteria.mito <- qc$subsets_Mito_percent > 10
  discard <- criteria.lib | criteria.nexprs | criteria.mito
  
  
  # ADAPTIVE CRITERIA
  criteria.lib2 <- isOutlier(qc$sum, log=TRUE, type="lower")
  criteria.nexprs2 <- isOutlier(qc$detected, log=TRUE, type="lower")
  criteria.mito2 <- isOutlier(qc$subsets_Mito_percent, type="higher")
  discard2 <- criteria.lib2 | criteria.nexprs2 | criteria.mito2
  # Summarize the number of cells removed for each reason. https://osca.bioconductor.org/quality-control.html#identifying-low-quality-cells
  #DataFrame(LibSize=sum(criteria.lib2), NExprs=sum(criteria.nexprs2), MitoProp=sum(criteria.mito2), Total=sum(discard2))
  
  
  #FILTER OUT CELLS
  filtered <- sce[,!discard2]
  sce.filtered <- filtered[!mito, ]
  
  
  #NORMALIZATION
  set.seed(100)
  clust <- quickCluster(sce.filtered) # This is needed to calculate size factors
  #table(clust.zeisel)
  #deconv.sf.zeisel <- calculateSumFactors(sce.final, cluster=clust.zeisel) # This is just to calculate a value called deconv.sf.zeisel, didn't apply on matrix
  #summary(deconv.sf.zeisel) # To see it
  sce.total <- computeSumFactors(sce.filtered, cluster=clust, min.mean=0.1) #Normalize the countMatrix
  sce.total <- logNormCounts(sce.total)
  
  
  #SELECT HIGHLY VARIABLE GENES
  dec.total <- modelGeneVar(sce.total)
  fit.zeisel <- metadata(dec.total)
  #plot(fit.zeisel$mean, fit.zeisel$var, xlab="Mean", ylab="Variance")
  #curve(fit.zeisel$trend(x), col="dodgerblue", add=TRUE, lwd=2) # visualizing the fit
  #dec.total[order(dec.total$bio, decreasing=TRUE),]
  hvg <- getTopHVGs(dec.total, var.threshold=0)
  sce.hvg <- sce.total[hvg,]
  return(sce.hvg)
}


sce.hvg <- pre.processing(sce)

#FIND INTERESTING GENES AND ORGANIZE########################################################################################
#1. CD133+SP+sphere-forming Msigdbr (cancer-stem-cell)
interested.genes <- function(ch) {
  c2.sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP")
  head(unique(c2.sets$gs_name))
  cscMsigdbr.sets <- c2.sets[grep(pattern ="BEIER_GLIOMA_STEM_CELL_UP|HARRIS_BRAIN_CANCER_PROGENITORS|GUENTHER_GROWTH_SPHERICAL_VS_ADHERENT_UP", c2.sets$gs_name), ]#this gave the desired gene sets
  #5. csc geneCards
  cscGeneCards.sets <- read.csv("GeneCards-GlioblastomaStemCell.csv", header = TRUE) 
  #6. csc AmiGo
  nscAmigo.sets <- read.delim("NeuralStemCell_AmiGo.csv", header = TRUE)[3]%>% group_by()  
  #2. angiogenesis Msigdbr
  c5.sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
  head(unique(c5.sets$gs_name))
  angioMsigdbr.sets <- c5.sets[grep(pattern ="GO_SPROUTING_ANGIOGENESIS", c5.sets$gs_name), ]
  #3. angiogenesis GeneCards
  angioGeneCards.sets <- read.csv("GeneCards_Angiogenesis.csv", header = TRUE)[1:2000, ]
  #4. drug resistance ABC 
  ABCGeneCards.sets <- read.csv("GeneCards_ABC.csv", header = TRUE) 
  #5. hypoxia MGI
  hypoMGI.sets <- read.csv("MGI_ResponseToHypoxia.csv", header = TRUE)
  hypoMGI.sets[1] <- apply(hypo.sets[1],2,toupper)
  
  ###matching with our data sets
  hvg <- row.names(counts(sce.hvg))
  csc.msigdbr <- semi_join(cscMsigdbr.sets, as.data.frame(hvg), by = c("human_gene_symbol"="hvg"))
  csc.GeneCards <- semi_join(cscGeneCards.sets, as.data.frame(hvg) , by = c("Gene.Symbol"="hvg"))
  nsc.Amigo <- semi_join(nscAmigo.sets, as.data.frame(hvg), by = c("c3"="hvg"))
  #csc_literature <- data.frame("PROM1", "SOX2", "")
  angio.msigdbr <- semi_join(angioMsigdbr.sets, as.data.frame(hvg), by = c("human_gene_symbol"="hvg"))
  angio.GeneCards <- semi_join(angioGeneCards.sets, as.data.frame(hvg), by = c("Gene.Symbol"= "hvg"))
  ABC.GeneCards <- semi_join(ABC.sets, as.data.frame(hvg), by = c("Gene.Symbol"="hvg"))
  hypo.MGI <- semi_join(hypo.sets, as.data.frame(hvg), by = c("Symbol"="hvg"))

  gene.list <- list(cscM = as.character(csc.msigdbr[["human_gene_symbol"]]), 
                    cscGC = as.character(csc.GeneCards[["Gene.Symbol"]]), 
                    nscA = as.character(nsc.Amigo[["c3"]]),   
                    angioM = angio_msigdbr[["human_gene_symbol"]],
                    ABCGC = as.character(ABC.GeneCards[["Gene.Symbol"]]),
                    hypoMG = as.character(hypo.MGI[["Symbol"]]) 
                    )
  return(gene.list[[ch]])
  
}

cscM <- interested.genes("cscM")
csc.GC_M <- unique(c(interested.genes("cscGC"), interested.genes("cscM")))
angioM <- interested.genes("angioM")
hvg <- row.names(counts(sce.hvg))


###

#DIMENSIONAL REDUCTION ON SELECTED GENES###################################################################################################


dimension.reduction <- function(selected.set) {
  sce.selected <- sce.hvg[selected.set,]
  set.seed(100) # See below.
  sce.selected <- runPCA(sce.selected, subset_row=selected.set)
  reducedDim(sce.selected, "PCA_25") <- reducedDim(sce.selected, "PCA")[,1:25]
  set.seed(00101001101)
  sce.selected <- runTSNE(sce.selected, dimred="PCA")
  return(sce.selected)
  #
}

sce.cscGC <- dimension.reduction(cscGC)
plotReducedDim(sce.cscGC, dimred="TSNE", colour_by = "PTPRZ1")


#PCAtools
#pca <- PCAtools::pca(counts(sce.GC))
#pca_result <- PCAtools::getLoadings(pca)
#PCAtools::getComponents(pca) #get pc no.


#CORRELATION MATRIX#####################################################################################################

#combine csc genes and angiogenesis genes

correlation.matrix <- function(selected.set1, selected.set2) {
  set1_set2 <- c(selected.set1, selected.set2)
  set1_set2.counts <- counts(sce.hvg[set1_set2,])
  matrix.together <- cor(t(set1_set2.counts), method = "pearson") 
  matrix.set1_set2 <- as.data.frame(matrix.together[selected.set1, selected.set2])  
  #cormat <- round(cor(csc.counts),2)
  #round(res, 2)
  return(matrix.set1_set2)
}


cscA_angioM.cor <- correlation.matrix(cscA, angioM)
cscGC_angioM.cor <- correlation.matrix(interested.genes("cscGC"), angioM)








#dCor##################################################################################
library("energy")
cscM.counts <- counts(sce[cscM,])
angioM.counts <- counts(sce[angioM,])

correlation <- dcor(t(cscM.counts), t(angioM.counts))



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

