d=$1
q=$2

echo pruning ...
./random_remove.py ../gene_trees/original/$d.gt.trees ../gene_trees/LBF/$d.gt.LBF.$q.removal.log > ../gene_trees/randprune/$d.gt.randprune.$q.trees
echo done

#echo quartet scoring by astral ...
#java -jar $astral -i ../gene_trees/randprune/$d.gt.randprune.$q.trees -o ../results/$d.randprune.$q.score.tre -q ../reference_trees/$d.ref.tre -t 1 > ../results/$d.randprune.$q.score.log 2<&1
#echo done

#echo computing local support by astral ...
#java -jar $astral -i ../gene_trees/randprune/$d.gt.randprune.$q.trees -o ../results/$d.randprune.$q.bs.tre -q ../reference_trees/$d.ref.tre  > ../results/$d.randprune.$q.bs.log 2<&1
#echo done

#echo extracting occupancy ...
#grep "occupancy" ../reults/$d.randprune.$q.score.log | sed -e "s/Taxon occupancy: {//g" -e "s/}//g" -e "s/=/ /g" | tr "," "\n"  > ../results/$d.randprune.$q.occupancy.dat
#echo done

echo computing MS ...
java -jar $treecmp -r ../reference_trees/$d.ref.tre -i ../gene_trees/randprune/$d.gt.randprune.$q.trees -o ../results/$d.randprune.$q.MS.dat -P -d ms
echo done

echo computing RF ...
compareTrees.missingBranch ../reference_trees/$d.ref.tre ../gene_trees/randprune/$d.gt.randprune.$q.trees -simplify > ../results/$d.randprune.$q.RF.dat
echo done

#echo comparing quartet score to original ...
#compareTrees ../reference_trees/$d.score.tre ../reference_trees/$d.randprune.$q.score.tre > ../results/$d.randpruneorig.vs.$q.dat
#echo done

#echo comparing local support to original ...
#compareTrees ../reference_trees/$d.score.tre ../reference_trees/$d.randprune.$q.score.tre > ../results/$d.randprune.orig.vs.$q.dat
#echo done

echo All done!
