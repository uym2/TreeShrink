#! /bin/bash

kshrinkDir=$1
trees=$2

nw_labels -I $trees | sort | uniq -c | sed "s/^[ ]*//g" > $kshrinkDir/species_occ.txt


#for x in $kshrinkDir/*.ratios; do
#    echo $x
#    base=$kshrinkDir/`basename $x .ratios`
#    paste $x $base.removals | sed "s/k=.*://g" | tr -d "\t" > $base.ratio.RS
#    get_species_max_value.py $base.ratio.RS $base.spMaxRatio
#done 

while read line; do
    arr=($line)
    species=${arr[1]}
    occ=${arr[0]}
    echo $species
    grep $species $kshrinkDir/*.spMaxRatio | awk '{print $2;}' > $kshrinkDir/sp.$species.txt
    count=`wc -l $kshrinkDir/sp.$species.txt | awk '{print $1;}'`
    while [ $count -lt $occ ]; do
        echo 1 >> $kshrinkDir/sp.$species.txt
        count=$((count+1))    
    done
done  < $kshrinkDir/species_occ.txt

cat $kshrinkDir/*spMaxRatio | awk '{print $2;}' > $kshrinkDir/all_maxRatios 
