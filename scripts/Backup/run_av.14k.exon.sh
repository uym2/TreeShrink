#! /bin/bash

#./org_analysis.sh av.14k.exon > Logs/av.14k.exon.org.log 2<&1 & 

q=0.01; n=med; m=med; f=lkernel; d=av.14k.exon; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
q=0.025; n=med; m=med; f=lkernel; d=av.14k.exon; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
q=0.05; n=med; m=med; f=lkernel; d=av.14k.exon; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
q=0.075; n=med; m=med; f=lkernel; d=av.14k.exon; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
q=0.1; n=med; m=med; f=lkernel; d=av.14k.exon; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &

m=default;d=av.14k.exon; ./RnR_analysis.sh $d $m > Logs/$d.RnR.$m.log 2<&1 &
m=Lone;d=av.14k.exon; ./RnR_analysis.sh $d $m > Logs/$d.RnR.$m.log 2<&1 &
m=Lhalf;d=av.14k.exon; ./RnR_analysis.sh $d $m > Logs/$d.RnR.$m.log 2<&1 &
m=Lquart;d=av.14k.exon; ./RnR_analysis.sh $d $m > Logs/$d.RnR.$m.log 2<&1 &
m=L3quart;d=av.14k.exon; ./RnR_analysis.sh $d $m > Logs/$d.RnR.$m.log 2<&1 &

k=2.0;r=mixed;d=av.14k.exon; ./rooted_prune_analysis.sh $d $r $k > Logs/$d.$r.$k.log 2<&1 &
k=2.5;r=mixed;d=av.14k.exon; ./rooted_prune_analysis.sh $d $r $k > Logs/$d.$r.$k.log 2<&1 &
k=3.0;r=mixed;d=av.14k.exon; ./rooted_prune_analysis.sh $d $r $k > Logs/$d.$r.$k.log 2<&1 &
k=3.5;r=mixed;d=av.14k.exon; ./rooted_prune_analysis.sh $d $r $k > Logs/$d.$r.$k.log 2<&1 &
k=4.0;r=mixed;d=av.14k.exon; ./rooted_prune_analysis.sh $d $r $k > Logs/$d.$r.$k.log 2<&1 &

d=av.14k.exon;m=med;q=0.01; ./randprune_analysis.sh $d $m $q > Logs/$d.rand.$m.$q 2<&1 &
d=av.14k.exon;m=med;q=0.025; ./randprune_analysis.sh $d $m $q > Logs/$d.rand.$m.$q 2<&1 &
d=av.14k.exon;m=med;q=0.05; ./randprune_analysis.sh $d $m $q > Logs/$d.rand.$m.$q 2<&1 &
d=av.14k.exon;m=med;q=0.075; ./randprune_analysis.sh $d $m $q > Logs/$d.rand.$m.$q 2<&1 &
d=av.14k.exon;m=med;q=0.1; ./randprune_analysis.sh $d $m $q > Logs/$d.rand.$m.$q 2<&1 &
