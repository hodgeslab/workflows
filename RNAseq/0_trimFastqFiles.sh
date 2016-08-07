WORKDIR=XXWORKDIRXX
OUTDIR=.

OFFSET=10
LEN=65

cd $WORKDIR

mkdir -p $OUTDIR

for i in *.fastq.gz; do
  NAME=`basename $i | sed s/.fastq.gz/_trimmed.fastq/`
  zcat $i | awk -v offset=$OFFSET -v len=$LEN 'NR % 2 == 0 { print substr($1, offset+1, len) } NR % 2 == 1' > $OUTDIR/$NAME
done
