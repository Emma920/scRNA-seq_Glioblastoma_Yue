
#!/usr/bin/env bash
FILES=/home/hdd/yue
#mkdir /home/ssd/yue/indices mkdir $FILES/STAR_aligned/results
for file in $(ls $FILES/data/trimmed/Darmanis/*1P.fastq) 
do
    bname=$(basename $file)
    echo HERE $bname
    echo HERE2 $file
    kallisto quant -t 15 -i $FILES/data/text/kallisto_transindex -o $FILES/data/aligned/Darmanis/kallisto  ${file%_*}_1P.fastq ${file%_*}_2P.fastq   
    
    mv $FILES/data/aligned/Darmanis/kallisto/abundance.h5 $FILES/data/aligned/Darmanis/kallisto/${bname%_*}.h5
    mv $FILES/data/aligned/Darmanis/kallisto/abundance.tsv $FILES/data/aligned/Darmanis/kallisto/${bname%_*}.tsv
    mv $FILES/data/aligned/Darmanis/kallisto/run_info.json $FILES/data/aligned/Darmanis/kallisto/${bname%_*}.json
#STAR --runThreadN 20 --runMode genomeGenerate --genomeDir $FILES/STAR_aligned/indices --genomeFastaFiles 
#$FILES/STAR_aligned/GRCh38.primary_assembly.genome.fa --sjdbGTFfile 
#$FILES/STAR_aligned/gencode.v33.primary_assembly.annotation.gtf
done