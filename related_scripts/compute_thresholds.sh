#! /bin/bash

args=("$@")
kshrinkDir=${args[0]}
#q1=${args[1]}


#Rscript find_threshold_kernel.R $kshrinkDir/all_maxRatios $q1 | tail -n1 > $kshrinkDir/geneThreshold.$q1.txt
i=1
while [ $i -lt $# ]; do
    q=${args[$i]}
    for s in  $kshrinkDir/sp.*txt; do
        Rscript find_threshold_kernel.R $s $q | tail -n1 > $kshrinkDir/spThreshold.`basename $s | sed -e "s/^sp\.//g" -e "s/\.txt//g"`.$q.txt
    done
    i=$((i+1))
done    
