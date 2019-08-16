#!/bin/sh

WORKDIR=.
INPUTDIR=fragLibrary
FRAGMENTLENGTH=200
MERGEDIST=1000
SAMPLELIST=$1
EXPERIMENTNAME=$2
SITESBED=$3
OUTPUTDIR=$EXPERIMENTNAME
GENOME=hg38
BLACKLISTBED=~/data/UCSC_Downloads/${GENOME}/${GENOME}-blacklist.bed
MACSTHRESH=2

module load bedtools

cd $WORKDIR

mkdir -p $OUTPUTDIR

# use custom sites
cat $SITESBED | cut -f 1-3 | sort -k 1,1 -k 2,2n -k 3,3n -u > $OUTPUTDIR/mergedPeaks.bed

# calculate background window
# count background reads (use 100 bp window +/- 8 kbp from all TSSs)
DIST=8000
WINDOW=100
BGBED=$OUTPUTDIR/All_background.bed
sort -k 1,1 -k 2,2n -u ~/data/UCSC_Downloads/$GENOME/knownGene_${GENOME}.bed | \
  awk -v dist=$DIST -v win=$WINDOW 'BEGIN{FS=OFS="\t"}{if($2 > dist) print $1,$2-dist,$2-dist+win; print $1,$2+dist,$2+dist+win}' | sort -k 1,1 -k 2,2n -k 3,3n -u > $BGBED
BGSIZE=`awk 'BEGIN{FS="\t"; sum=0}{sum+=$3-$2}END{print sum}' $BGBED`

for SAMPLE in `cat $SAMPLELIST`; do
  echo "Processing $SAMPLE ..."
  FRAGMENTSBED=$INPUTDIR/${SAMPLE}.fragments_${FRAGMENTLENGTH}bp.bed
  MERGEDPEAKSBED=$OUTPUTDIR/mergedPeaks_count-${SAMPLE}.fragments_${FRAGMENTLENGTH}bp.bed
  BACKGROUNDBED=$OUTPUTDIR/mergedPeaks_background-${SAMPLE}.fragments_${FRAGMENTLENGTH}bp.bed
  BGDENSITYFILE=$OUTPUTDIR/mergedPeaks_bgDensity-${SAMPLE}.fragments_${FRAGMENTLENGTH}bp.txt
  COUNTMINUSBGBED=$OUTPUTDIR/mergedPeaks_countMinusBg-${SAMPLE}.fragments_${FRAGMENTLENGTH}bp.bed

  # count overlapping reads (add pseudo-count if zero reads are encountered)
  # bedtools intersect -sorted -c -a $WORKDIR/$OUTPUTDIR/mergedPeaks.bed -b $FRAGMENTSBED | awk 'BEGIN{ OFS="\t" } { pseudo=$NF > 0 ? $NF : 1; $NF=pseudo; print }' > $MERGEDPEAKSBED
  if [ ! -r $MERGEDPEAKSBED ]; then
    bedtools intersect -sorted -c -a $WORKDIR/$OUTPUTDIR/mergedPeaks.bed -b $FRAGMENTSBED > $MERGEDPEAKSBED
  fi

  # count overlapping reads with background
  if [ ! -r $BACKGROUNDBED ] || [ ! -r $BGDENSITYFILE ]; then
    bedtools intersect -sorted -c -a $BGBED -b $FRAGMENTSBED > $BACKGROUNDBED
    awk -v bgsize=$BGSIZE 'BEGIN{FS="\t"; sum=0}{sum+=$NF}END{print sum/bgsize}' $BACKGROUNDBED > $BGDENSITYFILE
  fi

  BGDENSITY=`cat $BGDENSITYFILE`
  if [ ! -r $COUNTMINUSBGBED ]; then
    awk -v density=$BGDENSITY 'BEGIN{FS=OFS="\t"}{$NF=int($NF-($3-$2)*density); print}' $MERGEDPEAKSBED > $COUNTMINUSBGBED
  fi

done

# assemble master bed file with columns: id chr start end label counts ...

echo "Assembling data from all samples ..."

cumFile=$OUTPUTDIR/rawCounts.txt
headerFile=$OUTPUTDIR/rawCounts.header.txt
labelFile=$OUTPUTDIR/rawCounts.labels.txt

echo -en "id\tchr\tstart\tend\tlabel" > $headerFile
for SAMPLE in `cat $SAMPLELIST`; do
  echo -en "\t${SAMPLE}" >> $headerFile
done
echo >> $headerFile

awk 'BEGIN {FS="\t"; OFS="\t" } { print NR,$1,$2,$3,"" }' $OUTPUTDIR/mergedPeaks.bed > $cumFile
for SAMPLE in `cat $SAMPLELIST`; do
  MERGEDPEAKSBED=$OUTPUTDIR/mergedPeaks_count-${SAMPLE}.fragments_${FRAGMENTLENGTH}bp.bed
  COUNTMINUSBGBED=$OUTPUTDIR/mergedPeaks_countMinusBg-${SAMPLE}.fragments_${FRAGMENTLENGTH}bp.bed
  cut -f 4 $COUNTMINUSBGBED | paste $cumFile - > $cumFile.tmp
  mv $cumFile.tmp $cumFile
done

# mv $cumFile $cumFile.tmp

#
# add gene labels
#
GENEBED=~/data/UCSC_Downloads/${GENOME}/knownCanonical_RefSeq_${GENOME}.sort.bed
cut -f 2-4 $cumFile | bedtools closest -sorted -t first -a - -b $GENEBED | cut -f 8 > $labelFile
paste $cumFile $labelFile | \
awk 'BEGIN{FS="\t";OFS="\t"}{$5=$NF; for(i=1;i<=(NF-2);i++) printf $i "\t"; print $(NF-1)}' > $cumFile.tmp
cat $headerFile $cumFile.tmp > $cumFile
# rm $headerFile
# rm $cumFile.tmp

