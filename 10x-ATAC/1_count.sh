#!/bin/sh

GENOME=hg38
TEMPLATE=slurm_cellranger_count

mkdir -p scripts

WORKDIR=$(pwd)

cat sampleList.txt | sed '$a\' | while read -r line; do
  SAMPLENAME=$(echo "$line" | cut -f 1)
  SAMPLEID=$(echo "$line" | cut -f 2)
  FASTQDIR=$(echo "$line" | cut -f 3)

  sed -e s%XXSAMPLENAMEXX%$SAMPLENAME% -e s%XXSAMPLEIDXX%$SAMPLEID% -e s%XXWORKDIRXX%$WORKDIR% -e s%XXGENOMEXX%$GENOME% -e s%XXFASTQDIRXX%$FASTQDIR% $TEMPLATE > scripts/count_${SAMPLENAME}.slurm
done

cut -f 1 sampleList.txt > All_conditions.txt
