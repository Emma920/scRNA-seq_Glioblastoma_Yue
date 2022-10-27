  
#!/usr/bin/env bash

FILES=/home/hdd/yue/data
#mkdir /home/ssd/yue/indices
#mkdir $FILES/STAR_aligned/results
for file in $(ls $FILES/trimmed/Darmanis/*1P.fastq)
do
    bname=$(basename $file)
    echo HERE $bname
    echo HERE2 $file
    STAR --runThreadN 10 --genomeDir /home/ssd/yue/indices  --sjdbGTFfile $FILES/text/Homo_sapiens.GRCh38.99.gtf   --readFilesIn ${file%_*}_1P.fastq ${file%_*}_2P.fastq --quantMode TranscriptomeSAM GeneCounts --outFileNamePrefix $FILES/aligned/bam/Darmanis/STAR/${bname%_*} --soloFeatures Gene GeneFull SJ Velocyto
    #Rscript $FILES/STAR_aligned/samToCM.R
    #$FILES/aligned/bam/Darmanis/STAR/${bname%_*}Aligned.out.sam
    #$FILES/aligned/bam/Darmanis/STAR/${bname%_*}Aligned.toTranscriptome.out.bam
    #rm -r $FILES/aligned/bam/Darmanis/STAR/${bname%_*}_STARgenome

#STAR --runThreadN 20 --runMode genomeGenerate --genomeDir $FILES/STAR_aligned/indices --genomeFastaFiles $FILES/STAR_aligned/GRCh38.primary_assembly.genome.fa --sjdbGTFfile $FILES/STAR_aligned/gencode.v33.primary_assembly.annotation.gtf 
done

