  
#!/usr/bin/env bash
FILES=/home/hdd/yue
#mkdir /home/ssd/yue/indices mkdir $FILES/STAR_aligned/results
for file in $(ls $FILES/data/trimmed/20200313/*1P.fastq) 
do
    bname=$(basename $file)
    echo HERE $bname
    echo HERE2 $file
    STAR --runThreadN 20 --genomeDir /home/ssd/yue/indices --sjdbGTFfile $FILES/data/text/Homo_sapiens.GRCh38.99.ERCC.gtf --readFilesIn ${file%_*}_1P.fastq ${file%_*}_2P.fastq --quantMode TranscriptomeSAM GeneCounts --outFileNamePrefix $FILES/data/aligned/20200313/${bname%_*}  
   
#STAR --runThreadN 20 --runMode genomeGenerate --genomeDir $FILES/STAR_aligned/indices --genomeFastaFiles 
#$FILES/STAR_aligned/GRCh38.primary_assembly.genome.fa --sjdbGTFfile 
#$FILES/STAR_aligned/gencode.v33.primary_assembly.annotation.gtf
done
