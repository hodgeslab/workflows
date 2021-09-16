#!/bin/sh

WORKDIR=$(pwd)
GENOME=hg38
OUTDIR=fragLibrary
TEMPLATE=slurm_fragLibrary_paired

mkdir -p $OUTDIR
cut -f 1 sampleList.txt > All_conditions.txt

for CONDITION in $(cat All_conditions.txt | sed s/_rep[0-9]// | uniq); do
  sed -e s%XXWORKDIRXX%$WORKDIR% -e s%XXOUTDIRXX%$OUTDIR% -e s%XXCONDITIONXX%$CONDITION% -e s%XXGENOMEXX%$GENOME% $TEMPLATE > scripts/fragLibrary_${CONDITION}.slurm
done

