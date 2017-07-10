#! /bin/bash

treefile=$1
logfile=$2


#python count_removal.py $infofile $logfile

while read tree && read -u 3 log; do
     count=`echo $log | wc -w`
     if [ ! $count -eq 0 ]; then
         echo $tree | nw_prune - `echo $log`
     else
         echo $tree
     fi
done < $treefile 3< $logfile

