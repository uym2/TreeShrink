#! /bin/bash

treefile=$1
logfile=$2

taxafile=`mktemp`

#python count_removal.py $infofile $logfile

while read tree && read -u 3 log; do
     count=`echo $log | wc -w`
     if [ ! $count -eq 0 ]; then
         echo $tree | nw_labels -I - > $taxafile
         echo $tree | nw_prune - `python sampling.py $taxafile $count`
     else
         echo $tree
     fi
done < $treefile 3< $logfile

rm $taxafile      
