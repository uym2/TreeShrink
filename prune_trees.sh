#! /bin/bash

trees=$1
RS=$2
outtrees=$3

while read tree && read -u 3 rs; do
    if [[ -z $rs ]]; then
        echo $tree
    else
        echo $tree | nw_prune - `echo $rs`
    fi
done < $trees 3< $RS > $outtrees
