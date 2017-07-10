#! /bin/bash

d=mm

#for q in $(seq 0.005 0.005 0.15); do
#  #./lbf_analysis.sh $d $q > Logs/$d.lbf.$q.log 2<&1 &
#  ./randprune_analysis.sh $d $q > Logs/$d.rand.$q.log 2<&1 &
#done

for f in $(seq 0 5 101); do
  ./RnR_analysis.sh $d $f > Logs/$d.rnr.$f.log 2<&1 &
done

#for k in $(seq 2.0 0.06 4.0); do
#  ./rooted_prune_analysis.sh $d GAL $k > Logs/$d.rp.$k.log 2<&1 &
#done  


#q=0.01; n=medl; m=med; f=lnorm; d=kp.c1c2; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
#q=0.025; n=medl; m=med; f=lnorm; d=kp.c1c2; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
#q=0.05; n=medl; m=med; f=lnorm; d=kp.c1c2; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
#q=0.075; n=medl; m=med; f=lnorm; d=kp.c1c2; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
#q=0.1; n=medl; m=med; f=lnorm; d=kp.c1c2; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &


#q=0.01; n=med; m=med; f=kernel; d=kp.c1c2; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
#q=0.025; n=med; m=med; f=kernel; d=kp.c1c2; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
#q=0.05; n=med; m=med; f=kernel; d=kp.c1c2; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
#q=0.075; n=med; m=med; f=kernel; d=kp.c1c2; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
#q=0.1; n=med; m=med; f=kernel; d=kp.c1c2; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &


#q=0.01; n=stsd; m=sts; f=kernel; d=kp.c1c2; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
#q=0.025; n=stsd; m=sts; f=kernel; d=kp.c1c2; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
#q=0.05; n=stsd; m=sts; f=kernel; d=kp.c1c2; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
#q=0.075; n=stsd; m=sts; f=kernel; d=kp.c1c2; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
#q=0.1; n=stsd; m=sts; f=kernel; d=kp.c1c2; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &


#q=0.01; n=ind; m=ind; f=lnorm; d=kp.c1c2; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
#q=0.025; n=ind; m=ind; f=lnormm; d=kp.c1c2; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
#q=0.05; n=ind; m=ind; f=lnorm; d=kp.c1c2; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
#q=0.075; n=ind; m=ind; f=lnorm; d=kp.c1c2; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
#q=0.1; n=ind; m=ind; f=lnorm; d=kp.c1c2; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &


#k=2.0;r=mixed;d=kp.c1c2; ./rooted_prune_analysis.sh $d $r $k > Logs/$d.$r.$k.log 2<&1 &
#k=2.5;r=mixed;d=kp.c1c2; ./rooted_prune_analysis.sh $d $r $k > Logs/$d.$r.$k.log 2<&1 &
#k=3.0;r=mixed;d=kp.c1c2; ./rooted_prune_analysis.sh $d $r $k > Logs/$d.$r.$k.log 2<&1 &
#k=3.5;r=mixed;d=kp.c1c2; ./rooted_prune_analysis.sh $d $r $k > Logs/$d.$r.$k.log 2<&1 &
#k=4.0;r=mixed;d=kp.c1c2; ./rooted_prune_analysis.sh $d $r $k > Logs/$d.$r.$k.log 2<&1 &
