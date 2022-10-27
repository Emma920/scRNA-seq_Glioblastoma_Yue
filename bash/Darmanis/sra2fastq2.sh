#!/usr/bin/env bash
FILES1=/home/yue/hdd/yue/data/ncbi/Darmanis_normal
FILES2=/home/yue/hdd/yue/data
for file in $(ls $FILES2/ncbi/Darmanis_normal/sra/*.sra) 
do
    parallel-fastq-dump --tmpdir $FILES1/fastq --sra-id $file --threads 10 --outdir $FILES1/fastq --split-files
    
    
done
