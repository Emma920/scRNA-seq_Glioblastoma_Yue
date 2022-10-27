#!/usr/bin/env bash
FILES=/home/hdd/yue/ncbi/public/sra/output 
LISTS=/home/hdd/yue/ncbi/public/sra/OPC.txt 
for f in `cat $LISTS` 
do
       echo ${f}
       bname=$(basename $f)
       [ ! -d $FILES/OPC ] && mkdir $FILES/OPC
       mv $FILES/${bname}_* $FILES/OPC
done
