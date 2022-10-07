#!/bin/sh

INFILE=$1
SAMPLES=$(cat $INFILE | cut -f 1)

for SAMPLE in $SAMPLES; do md5sum ${SAMPLE}*; done | awk 'BEGIN{FS=" "; OFS="\t";}{print $2, $1}' | sort -k 1,1 > ${INFILE}.md5
