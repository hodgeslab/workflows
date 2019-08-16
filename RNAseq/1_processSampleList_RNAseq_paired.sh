#!/bin/sh

GENOME=hg38
FASTQFOLDER=./fastq

TEMPLATE=slurm_mapAndCount_hisat2_paired

mkdir -p scripts

WORKDIR=`pwd`

while read -r line; do
  SAMPLENAME=`echo "$line" | cut -f 1`
  READ1=`echo "$line" | cut -f 2`
  READ2=`echo "$line" | cut -f 3`

  sed -e s%XXSAMPLENAMEXX%$SAMPLENAME% -e s%XXWORKDIRXX%$WORKDIR% -e s%XXGENOMEXX%$GENOME% -e s%XXREAD1XX%$READ1% -e s%XXREAD2XX%$READ2% $TEMPLATE -e s%XXFASTQFOLDERXX%$FASTQFOLDER% > scripts/mapAndCount_${SAMPLENAME}.slurm
done < sampleList.txt

cut -f 1 sampleList.txt > All_conditions.txt
