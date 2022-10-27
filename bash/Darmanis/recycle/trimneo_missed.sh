#!/usr/bin/env bash
FILES=/home/hdd/alex/ncbi/public/sra/missed_files_neo/missed_output
for file in $(ls $FILES/*.fastq) 
do
    bname=$(basename $file)
    echo HERE $bname
    echo HERE2 $file
    #/home/ema/FastQC/fastqc -t 4 -o $FILES/QC_output $file
    #[ ! -d $FILES/Trimmed_missed ] && mkdir $FILES/Trimmed_missed
    #/usr/bin/TrimmomaticPE
    trimmomatic PE -threads 18 ${file%_*}_1.fastq ${file%_*}_2.fastq $FILES/Trimmed_missed/${bname%_*}_1P.fastq $FILES/Trimmed_missed/${bname%_*}_1U.fastq $FILES/Trimmed_missed/${bname%_*}_2P.fastq $FILES/Trimmed_missed/${bname%_*}_2U.fastq LEADING:20 SLIDINGWINDOW:4:20 MINLEN:20
done
