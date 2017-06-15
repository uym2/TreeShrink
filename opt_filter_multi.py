#! /usr/bin/env python

from optimal_filter_lib import TreeFilter
from sys import argv
from math import sqrt
from subprocess import check_output,call
import argparse
from sys import stdout
from dendropy import Tree


def find_k(gradient_list,threshold):
    k = len(gradient_list)
    while k > 0:
        if gradient_list[k-1] > threshold:
            break
        k = k-1
    return k

parser = argparse.ArgumentParser()

parser.add_argument("-i","--input",required=True,help="input trees")
parser.add_argument("-m","--method",required=True,help="method: ind,med,sts")
parser.add_argument("-f","--function",required=True,help="a function to fit to data: lnorm, kernel,lkernel")
parser.add_argument("-o","--output",required=True,help="output trees")
parser.add_argument("-r","--removal",required=False,help="list of the removals")
parser.add_argument("-q","--quantile",required=False,help="the cut-off quantile of the gradient to be used as threshold")
parser.add_argument("-g","--gradient",required=False,help="list of the gradient of the diameter by level")

args = vars(parser.parse_args())

intree = args["input"]
outtree = open(args["output"],'a')
method = args["method"]
thres = args["quantile"] if args["quantile"] else "0.05"

pwd="/Users/uym2/my_gits/LongBranchFiltering/" # hack it for now!
Rfunction = pwd + "find_threshold_" + args["function"] + ".R"
print(Rfunction)

myfilters = []
mydata = []
datafile = check_output(["mktemp"]).rstrip()

with open(intree,"r") as f:
    i = 1
    for line in f:
        a_tree = Tree.get(data=line,schema="newick")
        a_filter = TreeFilter(ddpTree=a_tree)
        a_filter.optFilter()

        branch_list = []

        for br in a_filter.ddpTree.preorder_edge_iter():
            if br.tail_node is not None:
                branch_list.append(br.length)

        branch_list.sort()
        if len(branch_list)%2:
            med_br = branch_list[len(branch_list)/2]
        else: 
            med_br = (branch_list[len(branch_list)/2] + branch_list[len(branch_list)/2-1])/2

        data = [(a_filter.min_diams[i-1]-a_filter.min_diams[i])/med_br for i in range(1,len(a_filter.min_diams))]

        if method == "sts":
#            factor = sorted(data[-5:])[1]
            factor = min([x for x in data if x > 0])
            data = [d/factor for d in data]

        fout = open(datafile,"w" if method == "ind" else "a")     
        
        for x in data:
           fout.write(str(x) + "\n")
        fout.close()

        if method == "ind":
            fg = open(args["gradient"],'a') if args["gradient"] else None
            fr = open(args["removal"],'a') if args["removal"] else None
#            f=open(args["removal"],"a") if args["removal"]
#            f.write("Tree " + str(i) + "\n")
            if fr:
                fr.write("Tree " + str(i) + "\n")
            if fg:
                fg.write("Tree " + str(i) + "\n")
                for d in data:
                    fg.write(str(d) + "\t")
                fg.write("\n")
            i = i + 1
#            opt_k = int(check_output(["Rscript","/Users/uym2/my_gits/LongBranchFiltering/find_d.R",datafile,thres])[4:].rstrip())
#            opt_t = float(check_output(["Rscript","/Users/uym2/my_gits/LongBranchFiltering/find_d.R",datafile,thres])[4:].rstrip())
            opt_t=float(check_output(["Rscript",Rfunction,datafile,thres]).lstrip().rstrip()[5:])
            opt_k=find_k(data,opt_t)
            fTree = a_filter.filterOut(d=opt_k, fout=fr)
            outtree.write(fTree.as_string("newick"))
        else:
            myfilters.append(a_filter)
            mydata.append(data)        

if method != "ind":
    fr = open(args["removal"],'a') if args["removal"] else None
    fg = open(args["gradient"],'a') if args["gradient"] else None
    
    opt_t=float(check_output(["Rscript",Rfunction,datafile,thres]).lstrip().rstrip()[5:])
    for i in range(len(mydata)):
        if fr:
            fr.write("Tree " + str(i+1) + "\n")
        if fg:
            fg.write("Tree " + str(i+1) + "\n")
            for d in mydata[i]:
                fg.write(str(d) + "\t")
            fg.write("\n")
        
        opt_k = find_k(mydata[i],opt_t)
        '''
        k = len(mydata[i])
        while k > 0:
            if mydata[i][k-1] > opt_t:
                break
            k = k-1
            
        if k > 0:
            fTree = myfilters[i].filterOut(d=k,fout=fr)
        else:
            fTree = myfilters[i].ddpTree
        '''            
        fTree = myfilters[i].filterOut(d=opt_k,fout=fr)

        outtree.write(fTree.as_string("newick"))
    if fr:
        fr.close()
    if fg:
        fg.close()
outtree.close()
call(["rm",datafile])
