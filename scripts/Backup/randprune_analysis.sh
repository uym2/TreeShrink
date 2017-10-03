d=$1
m=$2
q=$3

echo pruning ...
./random_remove.py ../gene_trees/original/$d.gt.trees ../gene_trees/LBF/$d.gt.LBF.$m.$q.removal.log > ../gene_trees/randprune/$d.gt.randprune.$m.$q.trees
echo done

echo quartet scoring by astral ...
java -jar $astral -i ../gene_trees/randprune/$d.gt.randprune.$m.$q.trees -o ../species_trees/$d.astral.randprune.$m.$q.score.tre -q ../species_trees/$d.astral.tre -t 1 > ../species_trees/$d.astral.randprune.$m.$q.score.log 2<&1
echo done

echo computing local support by astral ...
java -jar $astral -i ../gene_trees/randprune/$d.gt.randprune.$m.$q.trees -o ../species_trees/$d.astral.randprune.$m.$q.bs.tre -q ../species_trees/$d.astral.tre  > ../species_trees/$d.astral.randprune.$m.$q.bs.log 2<&1
echo done

echo extracting occupancy ...
grep "occupancy" ../species_trees/$d.astral.randprune.$m.$q.score.log | sed -e "s/Taxon occupancy: {//g" -e "s/}//g" -e "s/=/ /g" | tr "," "\n"  > ../results/$d.randprune.$m.$q.occupancy.dat
echo done

echo computing MS ...
java -jar $treecmp -r ../species_trees/$d.astral.tre -i ../gene_trees/randprune/$d.gt.randprune.$m.$q.trees -o ../results/$d.randprune.$m.$q.MS.dat -P -d ms
echo done

echo computing RF ...
compareTrees.missingBranch ../species_trees/$d.astral.tre ../gene_trees/randprune/$d.gt.randprune.$m.$q.trees -simplify > ../results/$d.randprune.$m.$q.RF.dat
echo done

echo comparing quartet score to original ...
compareTrees ../species_trees/$d.astral.score.tre ../species_trees/$d.astral.randprune.$m.$q.score.tre > ../results/$d.randpruneorig.vs.$m.$q.dat
echo done

echo comparing local support to original ...
compareTrees ../species_trees/$d.astral.score.tre ../species_trees/$d.astral.randprune.$m.$q.score.tre > ../results/$d.randprune.orig.vs.$m.$q.dat
echo done

echo All done!
