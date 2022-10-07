#!/bin/sh

WORKDIR=$(pwd)
GENOME=hg38
OUTDIR=fragLibrary
FRAGMENTLENGTH=200
TEMPLATE=slurm_fragLibrary

mkdir -p $OUTDIR
cat sampleList.txt | awk '{if(NF==3)print $1}' > All_conditions.txt

for CONDITION in $(cat All_conditions.txt | sed s/_rep[0-9]// | uniq); do
  sed -e s%XXWORKDIRXX%$WORKDIR% -e s%XXOUTDIRXX%$OUTDIR% -e s%XXCONDITIONXX%$CONDITION% -e s%XXGENOMEXX%$GENOME% $TEMPLATE > scripts/fragLibrary_${CONDITION}.slurm
done

