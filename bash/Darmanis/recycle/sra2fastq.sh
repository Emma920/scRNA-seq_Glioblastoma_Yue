#!/usr/bin/env bash
FILES=/home/yue/hdd/yue/data/ncbi/Darmanis_normal 
for file in $(ls $FILES/sra/*.sra) 
do
    fastq-dump --split-files $file -O $FILES/fastq
    rm $file
    
done
