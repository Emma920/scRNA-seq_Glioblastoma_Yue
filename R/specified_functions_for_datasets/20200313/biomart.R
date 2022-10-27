
biomart.meta <- function() {
  library(biomaRt)
  #listMarts()
  ensembl_hs_mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  
  #listDatasets(ensembl_hs_mart)[1:100,]
  #listAttributes(ensembl_hs_mart)[1:100,]
  #listFilters(ensembl_hs_mart)[1:100,]
  
  gene.meta <- getBM(attributes=c("ensembl_gene_id", "external_gene_name", "name_1006","definition_1006"), mart=ensembl_hs_mart)
  return(gene.meta)
}



biomart.name <- function() {
  library(biomaRt)
  #listMarts()
  ensembl_hs_mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  
  #listDatasets(ensembl_hs_mart)[1:100,]
  #listAttributes(ensembl_hs_mart)[1:100,]
  #listFilters(ensembl_hs_mart)[1:100,]
  
  gene.name <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"), mart=ensembl_hs_mart)
  return(gene.name)
}

