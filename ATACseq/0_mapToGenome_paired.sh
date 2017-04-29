#!/bin/sh

GENOME=mm9
WORKDIR=`pwd`
TEMPLATE=~/templates/workflows/ATACseq/bowtie2_paired_ATAC

mkdir -p scripts

while read -r line; do
  SAMPLENAME=`echo "$line" | cut -f 1`
  READ1=`echo "$line" | cut -f 2`
  READ2=`echo "$line" | cut -f 3`

  sed -e s%XXWORKDIRXX%$WORKDIR% -e s%XXSAMPLENAMEXX%$SAMPLENAME% -e s%XXREAD1XX%$READ1% -e s%XXREAD2XX%$READ2% -e s%XXGENOMEXX%$GENOME% $TEMPLATE > scripts/bowtie_${SAMPLENAME}
done < sampleList.txt
