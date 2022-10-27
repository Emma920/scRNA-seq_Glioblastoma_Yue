#!/usr/bin/env bash

FILES=/home/hdd/yue/data/
for file in $(ls $FILES/fastq/20200313/*.fastq.gz) 
do
    fastqc -t 20 -o $FILES/qc/20200313 $file
    
done
