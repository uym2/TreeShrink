#! /bin/bash

intrees=$1
outdir=$2
mode=$3
start=$4
inc=$5
end=$6
name=$7

q=$(seq $start $inc $end)
treeshrink2.py -i $intrees -q "$q" -d $outdir -m $mode

orgMS=`dirname $intrees`/`basename $intrees .trees`.MS
orgOCC=`dirname $intrees`/`basename $intrees .trees`.occ

[ -s $orgOCC ] || nw_stats -fl $intrees | awk '{print $3;}' | numlist -sum > $orgOCC
[ -s $orgMS ] || java -jar $treecmp -i $intrees -o $orgMS -d ms -m -P


for t in $outdir/*trees; do echo $t; java -jar $treecmp -i $t -o $outdir/`basename $t .trees`.MS -d ms -m -P; done
for t in $outdir/*.MS; do paste $orgMS $t | awk '{print $7-$14}' > $outdir/`basename $t .MS`.dMS; done
for t in $outdir/*trees; do nw_stats -fl $t | awk '{print $3;}' | numlist -sum > $outdir/`basename $t .trees`.occ; done
for t in $outdir/*occ; do paste $orgOCC $t | awk '{print $2/$1;}' > $outdir/`basename $t .occ`.occ_norm; done

grep . $outdir/*occ_norm | sed -e "s/.occ_norm:/ /g" > $outdir/occ.dat
grep . $outdir/*dMS | sed -e "s/.dMS:/ /g" > $outdir/dMS.dat

join $outdir/occ.dat $outdir/dMS.dat | sed "s/^.*shrinked_/TreeShrink_$name q/g" > $outdir/occ_dMS.dat
