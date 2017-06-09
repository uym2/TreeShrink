#! /usr/bin/env python

from optimal_filter_lib import TreeFilter
from sys import argv
from math import sqrt
from subprocess import check_output,call
import argparse
from sys import stdout
from dendopy import Tree

parser = argparse.ArgumentParser()

parser.add_argument("-i","--input",required=True,help="input trees")
parser.add_argument("-m","--method",required=True,help="method: ind,med,sts")
parser.add_argument("-o","--output",required=False,help="output trees")
parser.add_argument("-r","--removal",required=False,help="list of the removals")
parser.add_argument("-t","--threshold",required=False,help="the cut-off threshold of the gradient")

args = vars(parser.parse_args())

intree = args["input"]
outtree = open(args["output"],'a')
method = args["method"]
thres = args["threshold"]

myfilters = []
mydata = []
datafile = check_output(["mktemp"]).rstrip()

with open(intree,"r") as f:
    for line in f:
        a_tree = Tree.get(data=line,schema="newick")
        a_filter = TreeFilter(ddpTree=a_tree)
        a_filter.optFilter()

        branch_list = []

        for br in myfilter.ddpTree.preorder_edge_iter():
            if br.tail_node is not None:
                branch_list.append(br.length)

        branch_list.sort()
        if len(branch_list)%2:
            med_br = branch_list[len(branch_list)/2]
        else: 
            med_br = (branch_list[len(branch_list)/2] + branch_list[len(branch_list)/2-1])/2

        data = [(a_filter.min_diams[i-1]-a_filter.min_diams[i])/med_br for i in range(1,len(a_filter.min_diams)]

        if method == "sts":
            factor = sorted(data[-3:])[1]
            data = [d/factor for d in data]

        fout = open(datafile,"w" if method == "ind" else "a")     
        
        for x in data:
           fout.write(str(d) + "\n")
        fout.close()

        if method == "ind":
            opt_k = int(check_output(["Rscript","/Users/uym2/my_gits/LongBranchFiltering/find_d.R",datafile,thres])[4:].rstrip())

            fTree = myfilter.filterOut(d=opt_k, fout=open(args["removal"],"a") if args["removal"] else stdout)
            outtree.write(fTree.as_string("newick") + "\n")
        else:
            myfilters.append(a_filter)
            mydata.append(data)        

if method != "ind":
    opt_t = check_output(["Rscript","/Users/uym2/my_gits/LongBranchFiltering/find_threshold.R",datafile,thres])
    for i in range(len(mydata)):
        k = len(mydata[i])
        while k > 0:
            if mydata[i][k-1] > opt_t:
                break
            k = k-1
    if k > 0:
        fTree = myfilters[i].filterOut(d=k,fout=open(args["removal"],"a") if args["removal"] else stdout)
    else:
        fTree = myfilters[i].ddpTree
    outtree.write(fTree.as_string("newick") + "\n")

outtree.close()
call(["rm",datafile])
