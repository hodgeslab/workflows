#!/bin/sh

GENOME=hg38
WORKDIR=$(pwd)
TEMPLATE=slurm_bowtie2_v2_paired
FASTQDIR=./fastq

mkdir -p scripts

cat sampleList.txt | sed '$a\' | while read -r line; do
  SAMPLENAME=$(echo "$line" | cut -f 1)
  READ1=$(echo "$line" | cut -f 2)
  READ2=$(echo "$line" | cut -f 3)

  sed -e s%XXWORKDIRXX%$WORKDIR% -e s%XXSAMPLENAMEXX%$SAMPLENAME% -e s%XXREAD1XX%$FASTQDIR/$READ1% -e s%XXREAD2XX%$FASTQDIR/$READ2% -e s%XXGENOMEXX%$GENOME% $TEMPLATE > scripts/bowtie2_${SAMPLENAME}.slurm
done
