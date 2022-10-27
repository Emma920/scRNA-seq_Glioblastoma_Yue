library(tidyverse)
library(plyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(msigdbr)

interested.genes <- function(ch, sce.object=sce.hvg) {
  c2.sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP")
  head(unique(c2.sets$gs_name))
  cscMsigdbr.sets <- c2.sets[grep(pattern ="BEIER_GLIOMA_STEM_CELL_UP|HARRIS_BRAIN_CANCER_PROGENITORS|GUENTHER_GROWTH_SPHERICAL_VS_ADHERENT_UP", c2.sets$gs_name), ]#this gave the desired gene sets
  #5. csc geneCards
  cscGeneCards.sets <- read.csv("/home/hdd/yue/data/text/FeatureList/GeneCards-GlioblastomaStemCell.csv", header = TRUE) 
  #6. csc AmiGo
  nscAmigo.sets <- read.delim("/home/hdd/yue/data/text/FeatureList/NeuralStemCell_AmiGo.csv", header = TRUE)[3]%>% group_by()  
  #2. angiogenesis Msigdbr
  c5.sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
  head(unique(c5.sets$gs_name))
  angioMsigdbr.sets <- c5.sets[grep(pattern ="GO_SPROUTING_ANGIOGENESIS", c5.sets$gs_name), ]
  #3. angiogenesis GeneCards
  angioGeneCards.sets <- read.csv("/home/hdd/yue/data/text/FeatureList/GeneCards_Angiogenesis.csv", header = TRUE)[1:2000, ]
  #4. drug resistance ABC 
  ABCGeneCards.sets <- read.csv("/home/hdd/yue/data/text/FeatureList/GeneCards_DrugResistanceABC.csv", header = TRUE) 
  #5. hypoxia MGI
  hypoMGI.sets <- read.csv("/home/hdd/yue/data/text/FeatureList/MGI_ResponseToHypoxia.csv", header = TRUE)
  hypoMGI.sets[1] <- apply(hypoMGI.sets[1],2,toupper)
  
  ###matching with our data sets
  hvg <- row.names(counts(sce.object))
  csc.msigdbr <- semi_join(cscMsigdbr.sets, as.data.frame(hvg), by = c("human_gene_symbol"="hvg"))
  csc.GeneCards <- semi_join(cscGeneCards.sets, as.data.frame(hvg) , by = c("Gene.Symbol"="hvg"))
  nsc.Amigo <- semi_join(nscAmigo.sets, as.data.frame(hvg), by = c("c3"="hvg"))
  #csc_literature <- data.frame("PROM1", "SOX2", "")
  angio.msigdbr <- semi_join(angioMsigdbr.sets, as.data.frame(hvg), by = c("human_gene_symbol"="hvg"))
  angio.GeneCards <- semi_join(angioGeneCards.sets, as.data.frame(hvg), by = c("Gene.Symbol"= "hvg"))
  ABC.GeneCards <- semi_join(ABCGeneCards.sets, as.data.frame(hvg), by = c("Gene.Symbol"="hvg"))
  hypo.MGI <- semi_join(hypoMGI.sets, as.data.frame(hvg), by = c("Symbol"="hvg"))
  
  gene.list <- list(cscM = as.character(csc.msigdbr[["human_gene_symbol"]]), 
                    cscGC = as.character(csc.GeneCards[["Gene.Symbol"]]), 
                    nscA = as.character(nsc.Amigo[["c3"]]),   
                    angioM = angio.msigdbr[["human_gene_symbol"]],
                    angioGC = as.character(angio.GeneCards[["Gene.Symbol"]]),
                    ABCGC = as.character(ABC.GeneCards[["Gene.Symbol"]]),
                    hypoMG = as.character(hypo.MGI[["Symbol"]]) 
  )
  return(gene.list[[ch]])
  
}