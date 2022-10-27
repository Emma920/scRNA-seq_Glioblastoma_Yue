#!/usr/bin/env bash

FILES=/home/hdd/yue/data/
for file in $(ls $FILES/trimmed/20200313/*P.fastq) 
do
    fastqc -t 20 -o $FILES/qc/20200313/qc_trimmed $file
    
done
