#!/usr/bin/env bash
mkdir /home/hdd/alex/ncbi/public/sra/output/QC_output 
FILES=/home/hdd/alex/ncbi/public/sra/output 
for file in $(ls $FILES/*.fastq) 
do
    fastqc -t 10 -o $FILES/QC_output $file
    
done
