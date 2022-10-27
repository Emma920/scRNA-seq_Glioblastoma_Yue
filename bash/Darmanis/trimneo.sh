#!/usr/bin/env bash

FILES=/home/yue/hdd/yue/data
for file in $(ls $FILES/ncbi/Darmanis_normal/fastq/*.fastq)
do
    bname=$(basename $file)
    echo HERE $bname
    echo HERE2 $file
    #/home/ema/FastQC/fastqc -t 4 -o $FILES/QC_output $file
    #[ ! -d $FILES/Trimmed ] && mkdir $FILES/Trimmed
    #/usr/bin/TrimmomaticPE  
    TrimmomaticPE -threads 18 ${file%_*}_1.fastq ${file%_*}_2.fastq $FILES/trimmed/Darmanis_normal_in_neoplastic_paper/${bname%_*}_1P.fastq $FILES/trimmed/Darmanis_normal_in_neoplastic_paper/${bname%_*}_1U.fastq $FILES/trimmed/Darmanis_normal_in_neoplastic_paper/${bname%_*}_2P.fastq $FILES/trimmed/Darmanis_normal_in_neoplastic_paper/${bname%_*}_2U.fastq LEADING:20 SLIDINGWINDOW:4:20 MINLEN:20
done


