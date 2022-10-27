
library(tidyverse)
library(plyr)
library(dplyr)
library(data.table)
library(Rsubread)

setwd('/home/hdd/alex/ncbi/public/sra/output/Trimmed/STAR_aligned')
dataList <- list.files( pattern="*.sam")
gtf <- list.files(path = "/home/hdd/alex/ncbi/public/sra/output/Trimmed/STAR_aligned", pattern = "*.gtf")
fc <- featureCounts(files=dataList, 
                    annot.ext=gtf, 
                    isGTFAnnotationFile=TRUE,
                    isPairedEnd=TRUE)  #remove [["counts"]] if also want to get info like length, start, end, str, and strand

fc1 <- as.data.frame(fc[["counts"]])
setwd('/home/hdd/alex/ncbi/public/sra/output/Trimmed/STAR_aligned/results')
write.table(fc1, file= sub("Aligned.out.sam","",dataList[1]))
