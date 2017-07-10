#! /bin/bash

treelist=$1
k=$2

a_tree=`mktemp`

while read line; do
    echo $line > $a_tree
    mean=`nw_distance $a_tree | numlist -avg`
    sd=`nw_distance $a_tree | numlist -std` # terrible solution!
    thres=`echo $mean+$sd*$k | bc`
    python prune_by_threshold.py $a_tree $thres
done < $treelist     
