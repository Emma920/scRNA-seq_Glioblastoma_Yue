#!/usr/bin/env bash

FILES=/home/hdd/yue/data/
for file in $(ls $FILES/fastq/20200313/*.fastq.gz)
do
    bname=$(basename $file)
    echo HERE $bname
    echo HERE2 $file
    #/usr/bin/TrimmomaticPE  
    trimmomatic PE -threads 18 ${file%_*}_R1.fastq.gz ${file%_*}_R2.fastq.gz $FILES/trimmed/20200313/${bname%_*}_1P.fastq $FILES/trimmed/20200313/${bname%_*}_1U.fastq $FILES/trimmed/20200313/${bname%_*}_2P.fastq $FILES/trimmed/20200313/${bname%_*}_2U.fastq ILLUMINACLIP:$FILES/text/adapters.fa:2:30:10 LEADING:15 TRAILING:15  MINLEN:36  
done