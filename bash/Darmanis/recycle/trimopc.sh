#!/usr/bin/env bash
FILES=/home/hdd/alex/ncbi/public/sra/OPC 
for file in $(ls $FILES/*.fastq) 
do
    bname=$(basename $file)
    echo HERE $bname
    echo HERE2 $file
    #/home/ema/FastQC/fastqc -t 4 -o $FILES/QC_output $file
    [ ! -d $FILES/Trimmed ] && mkdir $FILES/Trimmed
    trimmomatic PE -threads 18 ${file%_*}_1.fastq ${file%_*}_2.fastq $FILES/Trimmed/${bname%_*}_1P.fastq $FILES/Trimmed/${bname%_*}_1U.fastq $FILES/Trimmed/${bname%_*}_2P.fastq $FILES/Trimmed/${bname%_*}_2U.fastq LEADING:20 SLIDINGWINDOW:4:20 MINLEN:20 
done
