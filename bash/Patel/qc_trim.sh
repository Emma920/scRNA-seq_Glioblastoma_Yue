#!/usr/bin/env bash
#mkdir /home/hdd/alex/ncbi/public/sra/output/QC_output 
FILES=/home/hdd/yue/data/ 
for file in $(ls $FILES/trimmed/Patel/*P.fastq) 
do
    fastqc -t 10 -o $FILES/qc/Patel/qc_trim $file
    
done
