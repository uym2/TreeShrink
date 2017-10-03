#! /bin/bash

args=("$@")
trees=${args[0]}
kshrinkDir=${args[1]}
#q1=${args[1]}

i=2
t=1

while [ $i -lt $# ]; do
    q=${args[$i]}
    #geneThres=`Rscript find_threshold_kernel.R $kshrinkDir/all_maxRatios $q | tail -n1`
    #echo $geneThres > $kshrinkDir/geneThreshold.$q.txt
    while read -u 4 tree; do
        featureFile=$kshrinkDir/$t.spMaxRatio
        echo $t
        while read -u 3 line; do
            arr=($line)
            sp=${arr[0]}
            value=${arr[1]}
            for f in $kshrinkDir/spThreshold.$sp*; do
                q1=`basename $f | awk -F "." '{print $4;}'`
                spThres=`cat $f`
                outfile=$kshrinkDir/$t.gThres$q.spThres$q1.removals
                [ -e $outfile ] || touch $outfile
    #            thres=$spThres
    #            c=`echo $value '>' $geneThres '&&' $value '>' $spThres | bc -l`
    #            if [ $c -eq 1 ];then
    #                echo $sp >> $outfile
    #            fi
            done
        done 3< $featureFile
        for r in $kshrinkDir/$t.gThres$q.spThres*.removals ; do
            outtrees=$kshrinkDir/gThres$q.spThres.`basename $r .removals | awk -F "spThres." '{print $2;}'`.trees
            if [ -s $r ]; then
                echo $tree | nw_prune - `cat $r` >> $outtrees
            else
                echo $tree >> $outtrees
            fi
        done
        t=$((t+1))
    done 4< $trees
    i=$((i+1))
done    
