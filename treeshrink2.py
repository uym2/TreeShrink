#! /usr/bin/env python

from optimal_filter_lib import TreeFilter
from sys import argv
from math import sqrt
from subprocess import check_output,call
import argparse
from sys import stdout
from dendropy import Tree, TreeList
from os.path import basename, dirname, splitext
from os import mkdir
from copy import deepcopy

parser = argparse.ArgumentParser()

parser.add_argument("-i","--input",required=True,help="input trees")
parser.add_argument("-d","--outdir",required=False,help="output directory")
parser.add_argument("-o","--output",required=False,help="output trees")
#parser.add_argument("-r","--removal",required=False,help="the optimal removing sets by level")
#parser.add_argument("-d","--diameter",required=False,help="list of the optimal diameters by level")
#parser.add_argument("-a","--ratio",required=False,help="list of the ratios of the diameter by level")
parser.add_argument("-c","--centroid",required=False,action='store_true',help="do centroid reroot in preprocessing")
parser.add_argument("-k","--k",required=False,help="the maximum number of leaves that can be removed")
parser.add_argument("-q","--quantiles",required=False,help="the quantile(s) to set threshold")


args = vars(parser.parse_args())

quantiles = [ q for q in args["quantiles"].split()]
print(quantiles)

intrees = args["input"]
treeName,treeExt = splitext(basename(intrees))
outtrees = args["output"] if args["output"] else treeName + "_shrinked" + treeExt


k = int(args["k"]) if args["k"] else None

outdir = args["outdir"] if args["outdir"] else splitext(intrees)[0] + "_kshrink"
mkdir(outdir)

trees = TreeList.get_from_path(intrees,'newick')
gene_list = [[] for i in range(len(trees))]
species_map = {}
occ = {}
for t,a_tree in enumerate(trees):
    for n in a_tree.leaf_node_iter():
        s = n.taxon.label
        occ[s] = 1 if not s in occ else occ[s]+1

    a_filter = TreeFilter(ddpTree=a_tree,centroid_reroot=args["centroid"])
    a_filter.optFilter(d=k)

    mapping = {}
    for i in range(1,len(a_filter.min_diams)):
        r = a_filter.min_diams[i-1]/a_filter.min_diams[i]
        removals = a_filter.list_removals(d=i)
        for s in removals:
            mapping[s] = r if s not in mapping else max(mapping[s],r)

    for s in mapping:
        species_map[s] = [mapping[s]] if s not in species_map else species_map[s]+[mapping[s]]
        gene_list[t].append((s,mapping[s]))

for s in species_map:
    l = len(species_map[s])
    for i in range(occ[s]-l):
        species_map[s].append(1)
    filename = outdir+"/" + s + ".dat"
    with open(filename,'w') as f:
        for v in species_map[s]:
            f.write(str(v))
            f.write("\n")
    thresholds = [ 0 for i in range(len(quantiles)) ]        
    for i,q in enumerate(quantiles): 
        thresholds[i] = float(check_output(["Rscript","/Users/uym2/my_gits/TreeShrink/find_threshold_lkernel.R",filename,q]).lstrip().rstrip()[5:])
    species_map[s] = (species_map[s],thresholds)
#print(occ)
#print(gene_list)
#print(species_map)
removing_sets = [ [ [ ] for i in range(len(trees)) ] for j in range(len(quantiles)) ]

#with open(outdir + "/removing_sets.txt",'w') as f:
for t,gene in enumerate(gene_list):
    #f.write(str(t) + ": ")
    for s,r in gene:
        for i,threshold in enumerate(species_map[s][1]):
            if r > threshold:
                removing_sets[i][t].append(s)
                #f.write(s + " ")        
    #filt = lambda node: False if (node.taxon is not None and node.taxon.label in removing_sets[t]) else True
    #trees[t].filter_leaf_nodes(filt)
    #f.write("\n")

#trees.write_to_path(outtrees,'newick')

treeName,treeExt = splitext(outtrees)
for i,RS in enumerate(removing_sets):
    trees_shrinked = deepcopy(trees)
    for t,tree in enumerate(trees_shrinked):
        filt = lambda node: False if (node.taxon is not None and node.taxon.label in RS[t]) else True 
        tree.filter_leaf_nodes(filt)
    trees_shrinked.write_to_path(outdir + "/" + treeName + "_" + quantiles[i] + treeExt,'newick')   
#print(removing_sets)
