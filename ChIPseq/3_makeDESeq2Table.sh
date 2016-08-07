#!/bin/sh

PATH=$PATH:~/src/bedtools2/bin
WORKDIR=.
INPUTDIR=regularizedFragLibrary
FRAGMENTLENGTH=200
MERGEDIST=1000
SAMPLELIST=$1
EXPERIMENTNAME=$2
OUTPUTDIR=$EXPERIMENTNAME
GENOME=mm9
BLACKLISTBED=~/data/UCSC_Downloads/${GENOME}/${GENOME}-blacklist.bed
RPMTHRESH=10

cd $WORKDIR

mkdir -p $OUTPUTDIR

for SAMPLE in `cat $SAMPLELIST`; do
  cat $INPUTDIR/$SAMPLE.peaks.bed | awk -v rpmthresh=$RPMTHRESH 'BEGIN{FS="\t";OFS="\t"}{if($NF>rpmthresh) print}'
done | sort -k 1,1 -k 2,2n > $OUTPUTDIR/mergedPeaks.bed.tmp

bedtools merge -d $MERGEDIST -i $OUTPUTDIR/mergedPeaks.bed.tmp | grep -v chrM | bedtools intersect -v -wa -b $BLACKLISTBED -a - | sort -k 1,1 -k 2,2n > $OUTPUTDIR/mergedPeaks.bed
rm $OUTPUTDIR/mergedPeaks.bed.tmp

for SAMPLE in `cat $SAMPLELIST`; do
  echo "Processing $SAMPLE ..."
  FRAGMENTSBED=$INPUTDIR/${SAMPLE}.fragments_${FRAGMENTLENGTH}bp.bed
  MERGEDPEAKSBED=$OUTPUTDIR/mergedPeaks_count-${SAMPLE}.fragments_${FRAGMENTLENGTH}bp.bed

  # count overlapping reads (add pseudo-count if zero reads are encountered)
  # bedtools intersect -sorted -c -a $WORKDIR/$OUTPUTDIR/mergedPeaks.bed -b $FRAGMENTSBED | awk 'BEGIN{ OFS="\t" } { pseudo=$NF > 0 ? $NF : 1; $NF=pseudo; print }' > $MERGEDPEAKSBED
  if [ ! -r $MERGEDPEAKSBED ]; then
    bedtools intersect -sorted -c -a $WORKDIR/$OUTPUTDIR/mergedPeaks.bed -b $FRAGMENTSBED > $MERGEDPEAKSBED
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
  cut -f 4 $MERGEDPEAKSBED | paste $cumFile - > $cumFile.tmp
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

