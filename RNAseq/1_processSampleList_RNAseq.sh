#!/bin/sh

GENOME=hg38
FASTQFOLDER=./fastq

TEMPLATE=slurm_mapAndCount_hisat2

mkdir -p scripts

WORKDIR=$(pwd)

cat sampleList.txt | sed '$a\' | while read -r line; do
  SAMPLENAME=$(echo "$line" | cut -f 1)
  FASTQFILE=$(echo "$line" | cut -f 2)

  sed -e s%XXSAMPLENAMEXX%$SAMPLENAME% -e s%XXWORKDIRXX%$WORKDIR% -e s%XXGENOMEXX%$GENOME% -e s%XXFASTQFILEXX%$FASTQFILE% $TEMPLATE -e s%XXFASTQFOLDERXX%$FASTQFOLDER% > scripts/mapAndCount_${SAMPLENAME}.slurm
done

cut -f 1 sampleList.txt > All_conditions.txt
