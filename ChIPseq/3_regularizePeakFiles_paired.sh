#!/bin/sh

WORKDIR=$(pwd)
GENOME=hg38
OUTDIR=fragLibrary
TEMPLATE=slurm_fragLibrary_paired

mkdir -p $OUTDIR
cat sampleList.txt | awk '{if(NF==4)print $1}' > All_conditions.txt

for CONDITION in $(cat All_conditions.txt | sed s/_rep[0-9]// | uniq); do
  sed -e s%XXWORKDIRXX%$WORKDIR% -e s%XXOUTDIRXX%$OUTDIR% -e s%XXCONDITIONXX%$CONDITION% -e s%XXGENOMEXX%$GENOME% $TEMPLATE > scripts/fragLibrary_${CONDITION}.slurm
done

