  
#!/usr/bin/env bash

FILES=/root/Desktop/workspace/hdd/yue/data
#mkdir /home/ssd/yue/indices
#mkdir $FILES/STAR_aligned/results
for file in $(ls $FILES/trimmed/Darmanis_normal_in_neoplastic_paper/*1P.fastq)
do
    bname=$(basename $file)
    echo HERE $bname
    echo HERE2 $file
    STAR --runThreadN 10 --genomeDir /root/Desktop/workspace/ssd/yue/indices  --sjdbGTFfile $FILES/text/Homo_sapiens.GRCh38.99.gtf   --readFilesIn ${file%_*}_1P.fastq ${file%_*}_2P.fastq --quantMode TranscriptomeSAM GeneCounts --outFileNamePrefix $FILES/aligned/Darmanis_normal_in_neoplastic_paper/${bname%_*} 
    #Rscript $FILES/STAR_aligned/samToCM.R
    rm $FILES/aligned/Darmanis_normal_in_neoplastic_paper/${bname%_*}Aligned.out.sam
    rm $FILES/aligned/Darmanis_normal_in_neoplastic_paper/${bname%_*}Aligned.toTranscriptome.out.bam
    rm -r $FILES/aligned/Darmanis_normal_in_neoplastic_paper/${bname%_*}_STARgenome
    rm ${file%_*}_1P.fastq 
    rm ${file%_*}_2P.fastq
    rm ${file%_*}_1U.fastq 
    rm ${file%_*}_2U.fastq

#STAR --runThreadN 20 --runMode genomeGenerate --genomeDir $FILES/STAR_aligned/indices --genomeFastaFiles $FILES/STAR_aligned/GRCh38.primary_assembly.genome.fa --sjdbGTFfile $FILES/STAR_aligned/gencode.v33.primary_assembly.annotation.gtf 
done

















