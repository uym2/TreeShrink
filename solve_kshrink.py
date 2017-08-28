#! /usr/bin/env python

from optimal_filter_lib import TreeFilter
from sys import argv
from math import sqrt
from subprocess import check_output,call
import argparse
from sys import stdout
from dendropy import Tree
from os.path import basename, dirname, splitext
from os import mkdir
parser = argparse.ArgumentParser()

parser.add_argument("-i","--input",required=True,help="input trees")
parser.add_argument("-o","--outdir",required=False,help="output directory")
#parser.add_argument("-r","--removal",required=False,help="the optimal removing sets by level")
#parser.add_argument("-d","--diameter",required=False,help="list of the optimal diameters by level")
#parser.add_argument("-a","--ratio",required=False,help="list of the ratios of the diameter by level")
parser.add_argument("-c","--centroid",required=False,action='store_true',help="do centroid reroot in preprocessing")
parser.add_argument("-k","--k",required=False,help="the maximum number of leaves that can be removed")

args = vars(parser.parse_args())

intrees = args["input"]

k = int(args["k"]) if args["k"] else None

outdir = args["outdir"] if args["outdir"] else splitext(intrees)[0] + "_kshrink"
mkdir(outdir)

with open(intrees,"r") as f:
    t = 1
    for line in f:
        a_tree = Tree.get(data=line,schema="newick",preserve_underscores=True)
        a_filter = TreeFilter(ddpTree=a_tree,centroid_reroot=args["centroid"])
        a_filter.optFilter(d=k)

        with open(outdir+"/"+str(t)+".ratios","w") as fa:
            ratios = [(a_filter.min_diams[i-1]/a_filter.min_diams[i]) for i in range(1,len(a_filter.min_diams))]
            for r in ratios:
                fa.write(str(r)+"\n")

        with open(outdir+"/"+str(t)+".diams","w") as fd:
            for d in a_filter.min_diams:
                fd.write(str(d) + "\n")
        
        with open(outdir+"/"+str(t)+".removals","w") as fr:
            for i in range(1,len(a_filter.min_diams)):
                fr.write("k="+str(i) + ": ")
                a_filter.list_removals(d=i,fout=fr)
                fr.write("\n")
        t = t+1

