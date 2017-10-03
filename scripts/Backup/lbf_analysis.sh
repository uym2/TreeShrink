d=$1
m=$2
f=$3
q=$4
n=$5

#echo filtering ...
#~/my_gits/LongBranchFiltering/opt_filter_multi.py -i ../gene_trees/original/$d.gt.trees -o ../gene_trees/LBF/$d.gt.LBF.$n.$q.trees -g ../gene_trees/LBF/.gt.LBF.$n.$q.gradient.dat -r ../gene_trees/LBF/$d.gt.LBF.$n.$q.removal.log -m $m -q $q -f $f
#echo done

sed -i -e "s/'//g" -e "s/ /_/g" ../gene_trees/LBF/$d.gt.LBF.$n.$q.trees

echo quartet scoring by astral ...
java -jar $astral -i ../gene_trees/LBF/$d.gt.LBF.$n.$q.trees -o ../species_trees/$d.astral.LBF.$n.$q.score.tre -q ../species_trees/$d.astral.tre -t 1 > ../species_trees/$d.astral.LBF.$n.$q.score.log 2<&1
echo done

echo computing local support by astral ...
java -jar $astral -i ../gene_trees/LBF/$d.gt.LBF.$n.$q.trees -o ../species_trees/$d.astral.LBF.$n.$q.bs.tre -q ../species_trees/$d.astral.tre > ../species_trees/$d.astral.LBF.$n.$q.bs.log 2<&1
echo done

echo extracting occupancy ...
grep "occupancy" ../species_trees/$d.astral.LBF.$n.$q.score.log | sed -e "s/Taxon occupancy: {//g" -e "s/}//g" -e "s/=/ /g" | tr "," "\n"  > ../results/$d.$n.$q.occupancy.dat
echo done

echo computing MS ...
java -jar $treecmp -r ../species_trees/$d.astral.tre -i ../gene_trees/LBF/$d.gt.LBF.$n.$q.trees -o ../results/$d.$n.$q.MS.dat -P -d ms
echo done

echo computing RF ...
compareTrees.missingBranch ../species_trees/$d.astral.tre ../gene_trees/LBF/$d.gt.LBF.$n.$q.trees -simplify > ../results/$d.$n.$q.RF.dat
echo done

echo comparing quartet score to original ...
compareTrees ../species_trees/$d.astral.score.tre ../species_trees/$d.astral.LBF.$n.$q.score.tre > ../results/$d.orig.vs.$n.$q.dat
echo done

echo comparing local support to original ...
compareTrees ../species_trees/$d.astral.score.tre ../species_trees/$d.astral.LBF.$n.$q.score.tre > ../results/$d.orig.vs.$n.$q.dat
echo done

echo All done!
