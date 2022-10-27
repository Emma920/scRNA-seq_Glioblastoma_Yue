#!/usr/bin/env bash
FILES=/home/hdd/alex/ncbi/public/sra/missed_files_neo 
for file in $(ls $FILES/*.sra) 
do
    parallel-fastq-dump --sra-id $file --threads 8 --outdir $FILES/missed_output --split-files
    
done
