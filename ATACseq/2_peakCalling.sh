#!/bin/sh

GENOME=hg38
WORKDIR=$(pwd)
TEMPLATE=slurm_macs2

mkdir -p scripts

while read -r line; do
  SAMPLENAME=$(echo "$line" | cut -f 1)
  
  sed -e s%XXWORKDIRXX%$WORKDIR% -e s%XXSAMPLENAMEXX%$SAMPLENAME% -e s%XXGENOMEXX%$GENOME% $TEMPLATE > scripts/macs2_${SAMPLENAME}.slurm
done < sampleList.txt
