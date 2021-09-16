#!/bin/sh

# for example run as below:
# sh 1_makeTFheatmaps.sh ~/data/ATACseq ATAC_control_condition ATAC_treat_condition treat-increased

WORKDIR=$1
CNAME=$2
TNAME=$3

GENOME=hg38
FRAGMENTLENGTH=200

TFLIST="Jun-AP1 SpiB"

# $4 should be something like treat-decreased, treat-increased, etc.
INDIR=${WORKDIR}/${CNAME}-${TNAME}
BEDFILE=$INDIR/$4.bed

MASKBEDFILE=$INDIR/unchanged.bed
SIZEFACTORS=$INDIR/sizeFactors.txt
TFDIR=/s1/share/homer-tf-${GENOME}


CBIGWIG1=${WORKDIR}/fragLibrary/${CNAME}_rep1.fragments_${FRAGMENTLENGTH}bp.bigWig
CBIGWIG2=${WORKDIR}/fragLibrary/${CNAME}_rep2.fragments_${FRAGMENTLENGTH}bp.bigWig
TBIGWIG1=${WORKDIR}/fragLibrary/${TNAME}_rep1.fragments_${FRAGMENTLENGTH}bp.bigWig
TBIGWIG2=${WORKDIR}/fragLibrary/${TNAME}_rep2.fragments_${FRAGMENTLENGTH}bp.bigWig

module load R
module load bedtools

# bedtools intersect -a $BEDFILE -b $MASKBEDFILE -wa -v > differentialPeaks.bed
# BEDFILE=differentialPeaks.bed

for TF in ${TFLIST}; do
  NAME=$(basename $BEDFILE .bed)
  bedtools intersect -a $BEDFILE -b ${TFDIR}/${TF}.bed -wa -u > ${NAME}_${TF}.bed
  SITESBED=${NAME}_${TF}.bed
  SITESBEDNAME=$(basename $SITESBED .bed)

  sh makeMetaGeneMatrix_customsites.sh $CBIGWIG1 $CBIGWIG2 $CNAME $SITESBED
  sh makeMetaGeneMatrix_customsites.sh $TBIGWIG1 $TBIGWIG2 $TNAME $SITESBED

  Rscript plotDiffHeatmaps_customsites.R ${CNAME} ${TNAME} ${SIZEFACTORS} ${SITESBEDNAME}
  Rscript plotMetaGene_customsites.R ${CNAME} ${TNAME} ${SIZEFACTORS} ${SITESBEDNAME}

done
