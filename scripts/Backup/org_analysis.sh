d=$1

echo computing Astral species tree
java -jar $astral -i ../gene_trees/original/$d.gt.trees -o ../species_trees/$d.astral.tre > ../species_trees/$d.astral.log 2<&1
echo done

echo computing quartet score by astral ...
java -jar $astral -i ../gene_trees/original/$d.gt.trees -o ../species_trees/$d.astral.score.tre -q ../species_trees/$d.astral.tre -t 1 > ../species_trees/$d.astral.score.log 2<&1
echo done

echo computing local support by astral ...
java -jar $astral -i ../gene_trees/original/$d.gt.trees -o ../species_trees/$d.astral.bs.tre -q ../species_trees/$d.astral.tre > ../species_trees/$d.astral.bs.log 2<&1
echo done

echo extracting occupancy ...
grep "occupancy" ../species_trees/$d.astral.score.log | sed -e "s/Taxon occupancy: {//g" -e "s/}//g" -e "s/=/ /g" | tr "," "\n"  > ../results/$d.org.occupancy.dat
echo done

echo computing MS ...
java -jar $treecmp -r ../species_trees/$d.astral.tre -i ../gene_trees/original/$d.gt.trees -o ../results/$d.org.MS.dat -P -d ms
echo done

echo computing RF ...
compareTrees.missingBranch ../species_trees/$d.astral.tre ../gene_trees/original/$d.gt.trees -simplify > ../results/$d.org.RF.dat
echo done

echo All done!
