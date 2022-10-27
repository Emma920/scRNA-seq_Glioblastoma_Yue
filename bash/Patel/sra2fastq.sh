#!/usr/bin/env bash
FILES=/home/hdd/yue/data/ 
for file in $(ls $FILES/ncbi/Patel/*.sra) 
do
    parallel-fastq-dump --sra-id $file --threads 8 --outdir $FILES/fastq/Patel --split-files
    
done

