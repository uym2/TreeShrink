#! /usr/bin/env python

from optimal_filter_lib import TreeFilter
from sys import argv
from math import sqrt
from subprocess import check_output,call
import argparse
from sys import stdout
from dendropy import Tree
from os import remove,getcwd,path

def find_k(gradient_list,threshold):
    k = len(gradient_list)
    while k > 0:
        if gradient_list[k-1] > threshold:
            break
        k = k-1
    return k

parser = argparse.ArgumentParser()

parser.add_argument("-i","--input",required=True,help="input trees")
parser.add_argument("-m","--method",required=False,help="method: ind,med,sts. Default: med")
parser.add_argument("-f","--function",required=False,help="a function to fit to data: lnorm, kernel,lkernel. Default: lkernel")
parser.add_argument("-o","--output",required=False,help="output trees")
parser.add_argument("-r","--removal",required=False,help="the removing set")
parser.add_argument("-q","--quantile",required=False,help="the cut-off quantile of the gradient to be used as threshold")
parser.add_argument("-d","--diameter",required=False,help="list of the optimal diameter by level")
parser.add_argument("-g","--gradient",required=False,help="list of the gradient of the diameter by level")
parser.add_argument("-c","--centroid",required=False,action='store_true',help="do centroid reroot in preprocessing")

args = vars(parser.parse_args())

intree = args["input"]
outtree = open(args["output"],'a') if args["output"] else None
method = args["method"] if args["method"] else "med"
thres = args["quantile"] if args["quantile"] else "0.05"
function = args["function"] if args["function"] else "lkernel"

path = path.dirname(path.realpath(__file__))
Rfunction = path + "/find_threshold_" + function + ".R"

myfilters = []
mydata = []
datafile = check_output(["mktemp"]).rstrip()


fr = open(args["removal"],'w') if args["removal"] else stdout
fg = open(args["gradient"],'w') if args["gradient"] else None
fd = open(args["diameter"],'w') if args["diameter"] else None


with open(intree,"r") as f:
    i = 1
    for line in f:
        a_tree = Tree.get(data=line,schema="newick",preserve_underscores=True)
        a_filter = TreeFilter(ddpTree=a_tree,centroid_reroot=args["centroid"])
        a_filter.optFilter()


        if a_filter.ddpTree.seed_node.num_child_nodes() > 2:
            branch_list = [ch.edge_length for ch in a_filter.ddpTree.seed_node.child_node_iter()]
        else:
            branch_list = [sum([ch.edge_length for ch in a_filter.ddpTree.seed_node.child_node_iter()])]

        for br in a_filter.ddpTree.preorder_edge_iter():
            if br.tail_node is not None and br.tail_node is not a_filter.ddpTree.seed_node:
                branch_list.append(br.length)

        branch_list.sort()
        if len(branch_list)%2:
            med_br = branch_list[len(branch_list)//2]
        else: 
            med_br = (branch_list[len(branch_list)//2] + branch_list[len(branch_list)//2-1])/2

#        print(med_br)

        data = [(a_filter.min_diams[i-1]-a_filter.min_diams[i])/med_br for i in range(1,len(a_filter.min_diams))]

        if method == "sts":
#            factor = sorted(data[-5:])[1]
            factor = min([x for x in data if x > 0])
            data = [d/factor for d in data]

        fout = open(datafile,"w" if method == "ind" else "a")     
        
        for x in data:
           fout.write(str(x) + "\n")
        fout.close()

        if fg:
            for d in data:
                fg.write(str(d) + "\t")
            fg.write("\n")
        if fd:
            for d in a_filter.min_diams:
                fd.write(str(d) + "\t")
            fd.write("\n")

        if method == "ind":
#            fg = open(args["gradient"],'a') if args["gradient"] else None
#            fr = open(args["removal"],'a') if args["removal"] else stdout
#            f=open(args["removal"],"a") if args["removal"]
#            f.write("Tree " + str(i) + "\n")
#            if fr:
#                fr.write("Tree " + str(i) + "\n")
            i = i + 1
            opt_t=float(check_output(["Rscript",Rfunction,datafile,thres]).lstrip().rstrip()[5:])
            opt_k=find_k(data,opt_t)
            fTree = a_filter.filterOut(d=opt_k, fout=fr)
            fr.write("\n")
            if outtree:
                outtree.write(fTree.as_string("newick"))
        else:
            myfilters.append(a_filter)
            mydata.append(data)        

if method != "ind":
#    fr = open(args["removal"],'w') if args["removal"] else stdout
#    fg = open(args["gradient"],'w') if args["gradient"] else None
    
    opt_t=float(check_output(["Rscript",Rfunction,datafile,thres]).lstrip().rstrip()[5:])
#    print(opt_t)
    for i in range(len(mydata)):
        opt_k = find_k(mydata[i],opt_t)
        fTree = myfilters[i].filterOut(d=opt_k,fout=fr)
        fr.write("\n")
        if outtree:
            outtree.write(fTree.as_string("newick"))
if fr is not stdout:
    fr.close()
if fg:
    fg.close()
if fd:
    fd.close()
if outtree:
    outtree.close()
remove(datafile)    
#call(["rm",datafile])
