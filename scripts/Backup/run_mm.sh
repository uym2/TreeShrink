#! /bin/bash

q=0.01; n=medl; m=med; f=lnorm; d=mm; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
q=0.025; n=medl; m=med; f=lnorm; d=mm; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
q=0.05; n=medl; m=med; f=lnorm; d=mm; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
q=0.075; n=medl; m=med; f=lnorm; d=mm; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
q=0.1; n=medl; m=med; f=lnorm; d=mm; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &


q=0.01; n=medd; m=med; f=kernel; d=mm; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
q=0.025; n=medd; m=med; f=kernel; d=mm; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
q=0.05; n=medd; m=med; f=kernel; d=mm; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
q=0.075; n=medd; m=med; f=kernel; d=mm; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
q=0.1; n=medd; m=med; f=kernel; d=mm; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &


q=0.01; n=stsd; m=sts; f=kernel; d=mm; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
q=0.025; n=stsd; m=sts; f=kernel; d=mm; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
q=0.05; n=stsd; m=sts; f=kernel; d=mm; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
q=0.075; n=stsd; m=sts; f=kernel; d=mm; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &
q=0.1; n=stsd; m=sts; f=kernel; d=mm; ./lbf_analysis.sh $d $m $f $q $n > Logs/$d.$n.$q.log 2<&1 &


