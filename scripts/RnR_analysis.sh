d=$1
m=$2


#java -jar $astral -i ../gene_trees/rogueNarok/$d.RnR.$m.trees -o ../species_trees/$d.astral.RnR.$m.score.tre -q ../species_trees/$d.astral.tre -t 1 > ../species_trees/$d.astral.RnR.$m.score.log 2<&1

#grep "occupancy" ../species_trees/$d.astral.RnR.$m.score.log | sed -e "s/Taxon occupancy: {//g" -e "s/}//g" -e "s/=/ /g" | tr "," "\n"  > ../results/$d.RnR.$m.occupancy.dat

java -jar $treecmp -r ../reference_trees/$d.ref.tre -i ../gene_trees/rogueNarok/$d.RnR.$m.trees -o ../results/$d.RnR.$m.MS.dat -P -d ms

compareTrees.missingBranch ../reference_trees/$d.ref.tre ../gene_trees/rogueNarok/$d.RnR.$m.trees -simplify > ../results/$d.RnR.$m.RF.dat

#compareTrees ../species_trees/$d.astral.score.tre ../species_trees/$d.astral.RnR.$m.score.tre > ../results/$d.RnR.orig.vs.$m.dat
