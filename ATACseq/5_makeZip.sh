#!/bin/sh

INFILE="All_conditions.txt"

function all_vs_all () {
  iseen=""
  for irep in `cat $INFILE`; do
    iname=`sed s/_rep[0-9]*// <<< $irep`

    if grep -qw "$iname" <<< $iseen; then
      continue;
    fi

    jseen=""
    for jrep in `cat $INFILE`; do
      jname=`sed s/_rep[0-9]*// <<< $jrep`
      if ! (grep -qw "$jname" <<< $iseen || grep -qw "$jname" <<< $iname || grep -qw "$jname" <<< $jseen); then
        exptname="${iname}-${jname}";
        # assemble txt file from i and j
        echo "$exptname"
        grep -e "^${iname}_rep[0-9]*" $INFILE > ${exptname}.txt
        grep -e "^${jname}_rep[0-9]*" $INFILE >> ${exptname}.txt
      fi
    jseen="$jseen $jname"
    done
  iseen="$iseen $iname";
  done
}

# make all vs. all comparisons
EXPERIMENTS="All_conditions `all_vs_all`"

STR="qc"
for EXPERIMENT in $EXPERIMENTS; do
  STR="$STR $EXPERIMENT/*.{txt,pdf,tif} $EXPERIMENT/[tu]*.bed $EXPERIMENT/homer*"
done

echo zip -r archive.zip $STR
