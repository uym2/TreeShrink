#! /bin/bash

intrees=$1
outdir=$2
start=$3
inc=$4
end=$5


q=$(seq $start $inc $end)
python treeshrink2.py -i $intrees -q "$q" -d $outdir

orgMS=`mktemp`
orgOCC=`mktemp`
nw_stats -fl $intrees | awk '{print $3;}' | numlist -sum > $orgOCC
java -jar $treecmp -i $intrees -o $orgMS -d ms -m -P

cd $outdir

for t in *trees; do echo $t; java -jar $treecmp -i $t -o `basename $t .trees`.MS -d ms -m -P; done
for t in *.MS; do paste $orgMS $t | awk '{print $7-$14}' > `basename $t .MS`.dMS; done
for t in *trees; do nw_stats -fl $t | awk '{print $3;}' | numlist -sum > `basename $t .trees`.occ; done
for t in *occ; do paste $orgOCC $t | awk '{print $2/$1;}' > `basename $t .occ`.occ_norm; done

grep . *occ_norm | sed -e "s/.occ_norm:/ /g" > occ.dat
grep . *dMS | sed -e "s/.dMS:/ /g" > dMS.dat

join occ.dat dMS.dat | sed "s/^.*shrinked_/TreeShrink2 q/g" > occ_dMS.dat

rm $orgMS $orgOCC

