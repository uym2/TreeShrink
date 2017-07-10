#! /bin/bash

d=$1 
r=$2
k=$3

echo pruning
rooted_pruning.sh ../gene_trees/original/$d.gt.trees.rooted.$r $k > ../gene_trees/rooted_pruning/$d.gt.$r.$k.trees 

sed -i -e "s/'//g" -e "s/ /_/g" ../gene_trees/rooted_pruning/$d.gt.$r.$k.trees 

echo astral scoring
java -jar $astral -i ../gene_trees/rooted_pruning/$d.gt.$r.$k.trees -o ../species_trees/$d.astral.$r.$k.score.tre -t 1 -q ../species_trees/$d.astral.tre > ../species_trees/$d.astral.$r.$k.score.log 2<&1

echo extract occupancy
grep "occupancy" ../species_trees/$d.astral.$r.$k.score.log | sed -e "s/Taxon occupancy: {//g" -e "s/}//g" -e "s/=/ /g" | tr "," "\n"  > ../results/$d.$r.$k.occupancy.dat

echo compute MS
java -jar $treecmp -r ../species_trees/$d.astral.tre -i ../gene_trees/rooted_pruning/$d.gt.$r.$k.trees -o ../results/$d.$r.$k.MS.dat -P -d ms

echo compute RF 
compareTrees.missingBranch ../species_trees/$d.astral.tre ../gene_trees/rooted_pruning/$d.gt.$r.$k.trees -simplify > ../results/$d.$r.$k.RF.dat

echo extract quartet support and branch lenghts
compareTrees ../species_trees/$d.astral.score.tre ../species_trees/$d.astral.$r.$k.score.tre > ../results/$d.orig.vs.$r.$k.dat
