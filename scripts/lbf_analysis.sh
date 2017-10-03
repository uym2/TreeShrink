d=$1
q=$2

ftree=$d.gt.LBF.$q.trees

echo filtering ...
~/my_gits/LongBranchFiltering/opt_filter_multi.py -i ../gene_trees/original/$d.gt.trees -o ../gene_trees/LBF/$d.gt.LBF.$q.trees -g ../gene_trees/LBF/$d.gt.LBF.$q.gradient.dat -r ../gene_trees/LBF/$d.gt.LBF.$q.removal.log -q $q 
echo done

sed -i -e "s/'//g" -e "s/ /_/g" ../gene_trees/LBF/$d.gt.LBF.$q.trees

#echo quartet scoring by astral ...
#java -jar $astral -i ../gene_trees/LBF/$ftree -o ../results/$d.LBF.$q.score.tre -q ../reference_trees/$d.ref.tre -t 1 > ../results/$d.LBF.$q.score.log 2<&1
#echo done

#echo computing local support by astral ...
#java -jar $astral -i ../gene_trees/LBF/$ftree -o ../results/$d.LBF.$q.bs.tre -q ../reference_trees/$d.ref.tre > ../results/$d.LBF.$q.bs.log 2<&1
#echo done

#echo extracting occupancy ...
#grep "occupancy" ../results/$d.LBF.$q.score.log | sed -e "s/Taxon occupancy: {//g" -e "s/}//g" -e "s/=/ /g" | tr "," "\n"  > ../results/$d.$q.occupancy.dat
#echo done

echo computing MS ...
java -jar $treecmp -r ../reference_trees/$d.ref.tre -i ../gene_trees/LBF/$ftree -o ../results/$d.$q.MS.dat -P -d ms
echo done

echo computing RF ...
compareTrees.missingBranch ../reference_trees/$d.ref.tre ../gene_trees/LBF/$ftree -simplify > ../results/$d.$q.RF.dat
echo done

#echo comparing quartet score to original ...
#compareTrees ../reference_trees/$d.astral.score.tre ../reference_trees/$d.LBF.$q.score.tre > ../results/$d.orig.vs.$q.dat
#echo done

#echo comparing local support to original ...
#compareTrees ../reference_trees/$d.astral.score.tre ../reference_trees/$d.LBF.$q.score.tre > ../results/$d.orig.vs.$q.dat
#echo done

echo All done!
