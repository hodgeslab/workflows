#!/bin/sh

GENOME=mm9
TEMPLATE=~/templates/workflows/RNAseq/qsub_mapAndCount
COUNTSFILE=htseq_counts.txt

mkdir -p scripts

WORKDIR=`pwd`

while read -r line; do
  SAMPLENAME=`echo "$line" | cut -f 1`
  FASTQFILE=`echo "$line" | cut -f 2`

  sed -e s%XXSAMPLENAMEXX%$SAMPLENAME% -e s%XXWORKDIRXX%$WORKDIR% -e s%XXGENOMEXX%$GENOME% -e s%XXFASTQFILEXX%$FASTQFILE% $TEMPLATE > scripts/qsub_mapAndCount_${SAMPLENAME}
done < sampleList.txt

cut -f 1 sampleList.txt > All_conditions.txt
