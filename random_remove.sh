#! /bin/bash

# randomly prune the second trees to match the number of leaves in the first trees

fulltrees=$1
inducedtrees=$2
taxafile=`mktemp`


while read tree1 && read -u 3 tree2; do
    echo $tree1 | nw_labels -I - > $taxafile
    n1=`echo $tree1 | nw_stats -fl - | awk '{print $3;}'`
    n2=`echo $tree2 | nw_stats -fl - | awk '{print $3;}'`
    n=$((n1-n2))

    if [ ! $n -eq 0 ]; then
        echo $tree1 | nw_prune - `sampling.py $taxafile $n`
    else
        echo $tree1
    fi
done < $fulltrees 3< $inducedtrees

rm $taxafile      
