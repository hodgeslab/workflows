#!/bin/sh

INFILE=sampleList.txt
INDIR=./scripts

echo "sample	numReads	numMapped	numUnmapped	mapFrequency	processNum"

for NAME in $(cut -f 1 $INFILE); do
  PROCESSES=$(basename -a ${INDIR}/bowtie2_${NAME}.slurm.e* | sed 's/.*\.e//')

  for PROCESS in ${PROCESSES}; do
    FILE=${INDIR}/bowtie2_${NAME}.slurm.e${PROCESS}
    NREADS=$(grep "reads" ${FILE} | cut -f 1 -d " ")
    MAPFREQ=$(grep "overall" ${FILE} | cut -f 1 -d " ")
    NUNMAPPED=$(grep "aligned 0 times" ${FILE} | sed 's/^[ ]*//' | cut -f 1 -d " ")
    NMAPPED=$(expr ${NREADS} - ${NUNMAPPED})
    echo "${NAME}	${NREADS}	${NMAPPED}	${NUNMAPPED}	${MAPFREQ}	${PROCESS}"
  done
done
