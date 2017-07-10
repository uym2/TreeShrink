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
parser.add_argument("-r","--removal",required=False,help="the optimal removing sets by level")
parser.add_argument("-q","--quantile",required=False,help="the cut-off quantile of the gradient to be used as threshold")
parser.add_argument("-d","--diameter",required=False,help="list of the optimal diameters by level")
parser.add_argument("-g","--gradient",required=False,help="list of the gradients of the diameter by level")
parser.add_argument("-a","--ratio",required=False,help="list of the ratios of the diameter by level")
parser.add_argument("-c","--centroid",required=False,action='store_true',help="do centroid reroot in preprocessing")
parser.add_argument("-k","--k",required=False,help="the maximum number of leaves that can be removed")

args = vars(parser.parse_args())

intree = args["input"]

k = int(args["k"]) if args["k"] else None

fr = open(args["removal"],'w') if args["removal"] else stdout
fg = open(args["gradient"],'w') if args["gradient"] else None
fd = open(args["diameter"],'w') if args["diameter"] else None
fa = open(args["ratio"],'w') if args["ratio"] else None

with open(intree,"r") as f:
    for line in f:
        a_tree = Tree.get(data=line,schema="newick",preserve_underscores=True)
        a_filter = TreeFilter(ddpTree=a_tree,centroid_reroot=args["centroid"])
        a_filter.optFilter(d=k)


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


        if fa: 
            ratios = [(a_filter.min_diams[i-1]/a_filter.min_diams[i]) for i in range(1,len(a_filter.min_diams))]
            for r in ratios:
                fa.write(str(r)+"\t")
            fa.write("\n")
        if fg:
            gradients = [(a_filter.min_diams[i-1]-a_filter.min_diams[i])/med_br for i in range(1,len(a_filter.min_diams))]

            for g in gradients:
                fg.write(str(g) + "\t")
            fg.write("\n")

        if fd:
            for d in a_filter.min_diams:
                fd.write(str(d/med_br) + "\t")
            fd.write("\n")
        
        if fr:
            for i in range(1,len(a_filter.min_diams)):
                fr.write("k="+str(i) + ": ")
                a_filter.list_removals(d=i,fout=fr)
                fr.write("\n")

if fr is not stdout:
    fr.close()
if fg:
    fg.close()
if fd:
    fd.close()
if fa:
    fa.close()    
