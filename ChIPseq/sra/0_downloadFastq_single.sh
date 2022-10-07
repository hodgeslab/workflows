#!/bin/sh

FASTQDIR=./fastq
OUTDIR=./fastq
SRACODE=$1

module load sratoolkit

mkdir -p $OUTDIR
cd $OUTDIR

# single-end reads
fastq-dump --gzip -A ${SRACODE} > ${SRACODE}.o 2> ${SRACODE}.e &
