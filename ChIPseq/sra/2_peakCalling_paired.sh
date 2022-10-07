#!/bin/sh

GENOME=hg38
WORKDIR=$(pwd)
TEMPLATE=slurm_macs2_paired

mkdir -p scripts

while read -r line; do
  SAMPLENAME=$(echo "$line" | cut -f 1)
  TREATFILE=$(echo "$line" | cut -f 2)

  # only process samples that have controls
  # if [ ! -z ${CONTROLNAME} ]; then
  #   sed -e s%XXWORKDIRXX%$WORKDIR% -e s%XXSAMPLENAMEXX%$SAMPLENAME% -e s%XXTREATFILEXX%$TREATFILE% -e s%XXCONTROLNAMEXX%$CONTROLNAME% -e s%XXGENOMEXX%$GENOME% $TEMPLATE > scripts/macs2_${SAMPLENAME}.slurm
  #fi

  sed -e s%XXWORKDIRXX%$WORKDIR% -e s%XXSAMPLENAMEXX%$SAMPLENAME% -e s%XXTREATFILEXX%$TREATFILE% -e s%XXGENOMEXX%$GENOME% $TEMPLATE > scripts/macs2_${SAMPLENAME}.slurm

done < sampleList.txt
