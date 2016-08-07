#!/bin/sh

GENOME=$1
WORKDIR=`pwd`
TEMPLATE=~/templates/workflows/ChIPseq/macs_single

mkdir -p scripts

while read -r line; do
  SAMPLENAME=`echo "$line" | cut -f 1`
  TREATBED=`echo "$line" | cut -f 2`
  INPUTBED=`echo "$line" | cut -f 3`
  
  sed -e s%XXWORKDIRXX%$WORKDIR% -e s%XXSAMPLENAMEXX%$SAMPLENAME% -e s%XXTREATBEDXX%$TREATBED% -e s%XXINPUTBEDXX%$INPUTBED% -e s%XXGENOMEXX%$GENOME% $TEMPLATE > scripts/macs_${SAMPLENAME}
done < sampleList.txt
