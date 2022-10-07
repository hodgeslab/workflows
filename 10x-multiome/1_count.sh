#!/bin/sh

GENOME=mm10
TEMPLATE=slurm_cellranger_count
mkdir -p scripts

WORKDIR=$(pwd)
FASTQIN=${WORKDIR}/fastq_orig
FASTQ=${WORKDIR}/fastq
mkdir -p ${FASTQ}
mkdir -p ${FASTQ}/RNA
mkdir -p ${FASTQ}/ATAC

# build libraries.csv file
for SAMPLENAME in $(cat sampleList.txt | sed '$a\' | cut -f 1 | sort | uniq); do
  echo "fastqs,sample,library_type" > ${SAMPLENAME}.csv
done

cat sampleList.txt | sed '$a\' | while read -r line; do
  SAMPLENAME=$(echo "$line" | cut -f 1)
  LIBRARYID=$(echo "$line" | cut -f 2)
  LIBRARYTYPE=$(echo "$line" | cut -f 3)

  # link to new fastq file names
  for FILE in ${FASTQIN}/${LIBRARYID}*.fastq.gz; do
    SUFFIX=$(basename $FILE | sed 's/^.*_\(S[0-9]\{1,2\}_L[0-9]\{3\}.*\)/\1/')
    if [ "$LIBRARYTYPE" == "Gene Expression" ]; then
      ln -s ${FILE} ${FASTQ}/RNA/${SAMPLENAME}_${SUFFIX}
    elif [ "$LIBRARYTYPE" == "Chromatin Accessibility" ]; then
      ln -s ${FILE} ${FASTQ}/ATAC/${SAMPLENAME}_${SUFFIX}
    fi

  done

  if [ "$LIBRARYTYPE" == "Gene Expression" ]; then
    echo "${FASTQ}/RNA,${SAMPLENAME},${LIBRARYTYPE}" >> ${SAMPLENAME}.csv
  elif [ "$LIBRARYTYPE" == "Chromatin Accessibility" ]; then
    echo "${FASTQ}/ATAC,${SAMPLENAME},${LIBRARYTYPE}" >> ${SAMPLENAME}.csv
  fi

  sed -e s%XXSAMPLENAMEXX%$SAMPLENAME% \
      -e s%XXWORKDIRXX%$WORKDIR% \
      -e s%XXGENOMEXX%$GENOME% $TEMPLATE > scripts/count_${SAMPLENAME}.slurm

done

cat sampleList.txt | awk '{if(!seen[$1]) { seen[$1]++; print $1 }}' > All_conditions.txt
