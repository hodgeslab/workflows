#!/bin/sh

INFILE=sampleList.txt
INDIR=./scripts

echo "sample	numReads	numMapped	numUnmapped	mapFrequency	numFiltered	processNum"

for NAME in $(cut -f 1 $INFILE); do

  if [ "$(ls ${INDIR}/bowtie2_${NAME}.slurm.e*)" = "" ]; then
    continue;
  fi

  PROCESSES=$(basename -a ${INDIR}/bowtie2_${NAME}.slurm.e* | sed 's/.*\.e//')

  for PROCESS in ${PROCESSES}; do
    FILE=${INDIR}/bowtie2_${NAME}.slurm.e${PROCESS}
    NREADS=$(grep "reads" ${FILE} | cut -f 1 -d " " | head -1)
    MAPFREQ=$(grep "overall" ${FILE} | cut -f 1 -d " " | head -1)
    NUNMAPPED=$(grep "pairs aligned concordantly 0 times" ${FILE} | sed 's/^[ ]*//' | cut -f 1 -d " " | head -1)
    NMAPPED=$(expr ${NREADS} - ${NUNMAPPED})
    NFILTERED=$(wc -l ${NAME}.bed | cut -f 1 -d " ")
    echo "${NAME}	${NREADS}	${NMAPPED}	${NUNMAPPED}	${MAPFREQ}	${NFILTERED}	${PROCESS}"
  done
done
