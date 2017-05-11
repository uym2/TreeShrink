# ! /bin/bash

# two phases filter
# unrooted + rooted (default values)

infile=$1
outfile=$2

temp=`mktemp`

echo Phase1: unrooted filtering
unrooted_filter.py -i $infile -o $temp

echo Phase2: rooted filtering
rooted_filter.py $temp $outfile

rm $temp
