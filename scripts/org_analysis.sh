d=$1

#echo pruning ...
#./random_remove.py ../gene_trees/original/$d.gt.trees ../gene_trees/LBF/$d.gt.LBF.$q.removal.log > ../gene_trees/original/$d.gt.trees
#echo done

#echo quartet scoring by astral ...
#java -jar $astral -i ../gene_trees/original/$d.gt.trees -o ../results/$d.score.tre -q ../reference_trees/$d.ref.tre -t 1 > ../results/$d.score.log 2<&1
#echo done

#echo computing local support by astral ...
#java -jar $astral -i ../gene_trees/original/$d.gt.trees -o ../results/$d.bs.tre -q ../reference_trees/$d.ref.tre  > ../results/$d.bs.log 2<&1
#echo done

#echo extracting occupancy ...
#grep "occupancy" ../reults/$d.score.log | sed -e "s/Taxon occupancy: {//g" -e "s/}//g" -e "s/=/ /g" | tr "," "\n"  > ../results/$d.occupancy.dat
#echo done

echo computing MS ...
java -jar $treecmp -r ../reference_trees/$d.ref.tre -i ../gene_trees/original/$d.gt.trees -o ../results/$d.MS.org.dat -P -d ms
echo done

echo computing RF ...
compareTrees.missingBranch ../reference_trees/$d.ref.tre ../gene_trees/original/$d.gt.trees -simplify > ../results/$d.RF.org.dat
echo done

#echo comparing quartet score to original ...
#compareTrees ../reference_trees/$d.score.tre ../reference_trees/$d.score.tre > ../results/$d.originalorig.vs.$q.dat
#echo done

#echo comparing local support to original ...
#compareTrees ../reference_trees/$d.score.tre ../reference_trees/$d.score.tre > ../results/$d.original.orig.vs.$q.dat
#echo done

echo All done!
