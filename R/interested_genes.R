library(tidyverse)
library(plyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(msigdbr)

interested.genes <- function(ch, sce.object) {
  
  #1. csc Msigd
  c2.sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP")
  head(unique(c2.sets$gs_name))
  cscMsigdbr.sets <- c2.sets[grep(pattern ="BEIER_GLIOMA_STEM_CELL_UP|HARRIS_BRAIN_CANCER_PROGENITORS|GUENTHER_GROWTH_SPHERICAL_VS_ADHERENT_UP", c2.sets$gs_name), ]#this gave the desired gene sets
  #2. nsc AmiGo
  nscAmigo.sets <- read.delim("/home/hdd/yue/data/text/FeatureList/NeuralStemCell_AmiGo.csv", header = TRUE)[3]%>% group_by()  
  #3. glial cell differentiation Amigo
  diffAmigo.sets <- read.csv("/home/hdd/yue/data/text/FeatureList/Amigo_GlialCellDifferentiation.csv", header = TRUE)[3]%>% group_by()  
  #4. positive regulation of angiogenesis Amigo
  angioAmigo.sets <- read.csv("/home/hdd/yue/data/text/FeatureList/Amigo_PositiveRegulationOfAngiogenesis.csv", header = TRUE)[3]%>% group_by()  
  #5. GO_SPROUTING_ANGIOGENESIS Msigdb
  c5.sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
  head(unique(c5.sets$gs_name))
  angioMsigdbr.sets <- c5.sets[grep(pattern ="GO_SPROUTING_ANGIOGENESIS", c5.sets$gs_name), ]
  #5. regulation of cell cycle Amigo
  cycleAmigo.sets <- read.csv("/home/hdd/yue/data/text/FeatureList/Amigo_RegulationOfCellCycle.csv", header = TRUE)[3]%>% group_by()  
  #6. response to hypoxia
  hypoAmigo.sets <- read.csv("/home/hdd/yue/data/text/FeatureList/Amigo_ResponseToHypoxia.csv", header = TRUE)[3]%>% group_by()  
  #7. cancer migration HumanMine
  migrHM.sets <- read_tsv("/home/hdd/yue/data/text/FeatureList/HM_CancerMigration.tsv")[1]%>% group_by()  
  #8. cancer proliferation
  proHM.sets <- read_tsv("/home/hdd/yue/data/text/FeatureList/HM_CancerProliferation.tsv")[1]%>% group_by() 
  #9. neural stem cell
  nscHM.sets <- read_tsv("/home/hdd/yue/data/text/FeatureList/HM_NeuralStemCell.tsv")[1]%>% group_by()
  #10. stem cell self-renewal HumanMine
  srHM.sets <- read_tsv("/home/hdd/yue/data/text/FeatureList/HM_StemCellSelfRenewal.tsv")[1:100, 1]
  #11. G1 to G0 transition
  #G1G0Amigo.sets <- read_tsv("/home/hdd/yue/data/text/FeatureList/Amigo_G1ToG0Transion.tsv")[1]
  

  ###matching with our data sets
  genes <- row.names(counts(sce.object))
  csc.msigdbr <- semi_join(cscMsigdbr.sets, as.data.frame(genes), by = c("human_gene_symbol"="genes"))
  nsc.Amigo <- semi_join(nscAmigo.sets, as.data.frame(genes), by = c("c3"="genes"))
  diff.Amigo <- semi_join(diffAmigo.sets, as.data.frame(genes), by = c("c3"="genes"))
  angio.Amigo <- semi_join(angioAmigo.sets, as.data.frame(genes), by = c("c3"="genes"))
  angio.msigdbr <- semi_join(angioMsigdbr.sets, as.data.frame(genes), by = c("human_gene_symbol"="genes")) 
  cycle.Amigo <- semi_join(cycleAmigo.sets, as.data.frame(genes), by = c("c3"="genes"))
  hypo.Amigo <- semi_join(hypoAmigo.sets, as.data.frame(genes), by = c("c3"="genes"))
  migr.HM <- semi_join(migrHM.sets, as.data.frame(genes), by = c("Gene Symbol"="genes"))
  pro.HM <- semi_join(proHM.sets, as.data.frame(genes), by = c("Gene Symbol"="genes"))
  nsc.HM <- semi_join(nscHM.sets, as.data.frame(genes), by = c("Gene Symbol"="genes"))
  sr.HM <- semi_join(srHM.sets, as.data.frame(genes), by = c("Gene Symbol"="genes"))
  #G1G0.Amigo <- semi_join(G1G0Amigo.sets, as.data.frame(genes),by = c("c3"="genes"))
  
  #csc_literature <- data.frame("PROM1", "SOX2", "")


  cscM <- unique(as.character(csc.msigdbr[["human_gene_symbol"]]))
  nscA <- unique(as.character(nsc.Amigo[["c3"]]))
  nscHM <- unique(as.character(nsc.HM[["Gene Symbol"]]))
  diffA_0 <- unique(as.character(diff.Amigo[["c3"]]))
  diffA <- diffA_0[!(diffA_0 %in% cscM)]
  angioA_0 <- unique(as.character(angio.Amigo[["c3"]]))
  angioA <- angioA_0[!(angioA_0 %in% cscM)] 
  angioM_0 <- unique(as.character(angio.msigdbr[["human_gene_symbol"]]))
  angioM <- angioM_0[!(angioM_0 %in% cscM)]
  cycleA_0 <- unique(as.character(cycle.Amigo[["c3"]]))
  cycleA <- cycleA_0[!(cycleA_0 %in% cscM)]
  hypoA_0 <- unique(as.character(hypo.Amigo[["c3"]]))
  hypoA <- hypoA_0[!(hypoA_0 %in% cscM)]
  migrHM_0 <- unique(as.character(migr.HM[["Gene Symbol"]]))
  migrHM <- migrHM_0[!(migrHM_0 %in% cscM)]
  proHM_0 <- unique(as.character(pro.HM[["Gene Symbol"]]))
  proHM <- proHM_0[!(proHM_0 %in% cscM)]
  srHM_0 <- unique(as.character(sr.HM[["Gene Symbol"]]))
  srHM <- srHM_0[!(srHM_0 %in% cscM)]
  #G1G0A_0 <- unique(as.character(G1G0.Amigo[["c3"]]))
  #G1G0A <- G1G0A_0[!(G1G0A_0 %in% cscM)]
  
  
  gene.list <- list("cscM"=cscM, "nscA"=nscA, "diffA"=diffA, "angioA"=angioA, "angioM"=angioM, "cycleA"=cycleA, "hypoA"=hypoA, "migrHM"=migrHM, "proHM"=proHM, "nscHM"=nscHM, "srHM"=srHM) #, "G1G0A"=G1G0A)
  
  
  return(gene.list[[ch]])
  
}