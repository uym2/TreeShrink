#! /bin/bash



#q=0.01; n=med; m=med; f=lkernel; d=xen.can; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
#q=0.025; n=med; m=med; f=lkernel; d=xen.can; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
#q=0.05; n=med; m=med; f=lkernel; d=xen.can; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
#q=0.075; n=med; m=med; f=lkernel; d=xen.can; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
#q=0.1; n=med; m=med; f=lkernel; d=xen.can; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &

m=default;d=xen.can; ./RnR_analysis.sh $d $m > Logs/$d.RnR.$m.log 2<&1 &
m=Lone;d=xen.can; ./RnR_analysis.sh $d $m > Logs/$d.RnR.$m.log 2<&1 &
m=Lhalf;d=xen.can; ./RnR_analysis.sh $d $m > Logs/$d.RnR.$m.log 2<&1 &
m=Lquart;d=xen.can; ./RnR_analysis.sh $d $m > Logs/$d.RnR.$m.log 2<&1 &
m=L3quart;d=xen.can; ./RnR_analysis.sh $d $m > Logs/$d.RnR.$m.log 2<&1 &

#k=2.0;r=trees;d=xen.can; ./rooted_prune_analysis.sh $d $r $k > Logs/$d.$r.$k.log 2<&1 &
#k=2.5;r=trees;d=xen.can; ./rooted_prune_analysis.sh $d $r $k > Logs/$d.$r.$k.log 2<&1 &
#k=3.0;r=trees;d=xen.can; ./rooted_prune_analysis.sh $d $r $k > Logs/$d.$r.$k.log 2<&1 &
#k=3.5;r=trees;d=xen.can; ./rooted_prune_analysis.sh $d $r $k > Logs/$d.$r.$k.log 2<&1 &
#k=4.0;r=trees;d=xen.can; ./rooted_prune_analysis.sh $d $r $k > Logs/$d.$r.$k.log 2<&1 &

#./randprune_analysis.sh xen.can med 0.01 &
#./randprune_analysis.sh xen.can med 0.025 &
#./randprune_analysis.sh xen.can med 0.05 &
#./randprune_analysis.sh xen.can med 0.075 &
#./randprune_analysis.sh xen.can med 0.1 &
