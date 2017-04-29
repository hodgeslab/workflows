#!/bin/sh

GENOME=mm9
WORKDIR=`pwd`
TEMPLATE=~/templates/workflows/ATACseq/macs2_paired_ATAC

mkdir -p scripts

while read -r line; do
  SAMPLENAME=`echo "$line" | cut -f 1`
  TREATBED=`echo "${SAMPLENAME}.bed"`
  
  sed -e s%XXWORKDIRXX%$WORKDIR% -e s%XXSAMPLENAMEXX%$SAMPLENAME% -e s%XXTREATBEDXX%$TREATBED% -e s%XXGENOMEXX%$GENOME% $TEMPLATE > scripts/macs_${SAMPLENAME}
done < sampleList.txt
