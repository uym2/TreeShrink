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

parser.add_argument("-i","--input",required=True,help="Input trees")
parser.add_argument("-d","--outdir",required=False,help="Output directory")
parser.add_argument("-o","--output",required=False,help="Output trees")
parser.add_argument("-c","--centroid",required=False,action='store_true',help="Do centroid reroot in preprocessing")
parser.add_argument("-k","--k",required=False,help="The maximum number of leaves that can be removed")
parser.add_argument("-q","--quantiles",required=False,help="The quantile(s) to set threshold")
parser.add_argument("-m","--mode",required=False,help="Filtering mode: 'per-species', 'per-gene', 'all-genes'. Default: 'per-species'")

wdir = "/Users/uym2/my_gits/TreeShrink" 


args = vars(parser.parse_args())

quantiles = [ q for q in args["quantiles"].split()]
#print(quantiles)

intrees = args["input"]
treeName,treeExt = splitext(basename(intrees))
outtrees = args["output"] if args["output"] else treeName + "_shrinked" + treeExt

mode = args["mode"] if args["mode"] else 'per-species'

print(mode)

k = int(args["k"]) if args["k"] else None

outdir = args["outdir"] if args["outdir"] else splitext(intrees)[0] + "_kshrink"
mkdir(outdir)

trees = TreeList.get_from_path(intrees,'newick',preserve_underscores=True)
gene_list = [[] for i in range(len(trees))]
species_map = {}
occ = {}
removing_sets = [ [ [ ] for i in range(len(trees)) ] for j in range(len(quantiles)) ]

for t,a_tree in enumerate(trees):
    # solve k-shrink
    a_filter = TreeFilter(ddpTree=a_tree,centroid_reroot=args["centroid"])
    a_filter.optFilter(d=k)

    # compute species feature (i.e. the max ratio associated with each species for this gene tree)
    mapping = {}
    for i in range(1,len(a_filter.min_diams)):
        r = a_filter.min_diams[i-1]/a_filter.min_diams[i]
        removals = a_filter.list_removals(d=i)
        for s in removals:
            mapping[s] = r if s not in mapping else max(mapping[s],r)
    
    # gather per-species distributions and per-gene species features
    for s in mapping:
        if mode == 'per-species':
            species_map[s] = [mapping[s]] if s not in species_map else species_map[s]+[mapping[s]]
        if mode == 'per-species' or mode == 'all-genes':
            gene_list[t].append((s,mapping[s]))
    
    # fit kernel density to this gene's species features (per-gene mode)
    if mode == 'per-gene':
    	filename = outdir + "/" + "gene_" + str(t) + ".dat"
        with open(filename,'w') as f:
            for s in mapping:
                f.write(str(mapping[s]))
                f.write("\n")
            #n_missing = len(list(a_tree.leaf_node_iter())) - len(mapping)
            #for i in range(n_missing):
            #    f.write("1.0")
            #    f.write("\n")
        if len(mapping) > 1:
            for i,q in enumerate(quantiles):
                threshold = float(check_output(["Rscript",wdir + "/find_threshold_IQR.R",filename,q]).lstrip().rstrip()[4:]) 
                print("Threshold: ", threshold)
                for s in mapping:
                    if mapping[s] > threshold: 
                        removing_sets[i][t].append(s)
    # update taxon occupancy (only for per-species mode)
    if mode == 'per-species':
        for n in a_tree.leaf_node_iter():
            s = n.taxon.label
            occ[s] = 1 if not s in occ else occ[s]+1


# fit kernel density to the per-species distributions and compute per-species threshold (per-species mode)
if mode == 'per-species':
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
            thresholds[i] = float(check_output(["Rscript",wdir + "/find_threshold_lkernel.R",filename,q]).lstrip().rstrip()[5:])
        species_map[s] = (species_map[s],thresholds)

    for t,gene in enumerate(gene_list):
        for s,r in gene:
            for i,threshold in enumerate(species_map[s][1]):
                if r > threshold:
                    removing_sets[i][t].append(s)

# fit kernel density to all the species features across all genes and compute the global threshold (all-gene mode) 
if mode == 'all-genes':
    filename = outdir + "/" + "all_genes" + ".dat"
    with open(filename,'w') as f:
        for gene in gene_list:
            for s,r in gene:
                f.write(str(r))
                f.write("\n")
    for i,q in enumerate(quantiles):
        threshold = float(check_output(["Rscript",wdir + "/find_threshold_lkernel.R",filename,q]).lstrip().rstrip()[5:])
        for t,gene in enumerate(gene_list):
            for s,r in gene:
                if r > threshold:
                    removing_sets[i][t].append(s)


# the below code is locked now because Dendropy's filter_leaf_nodes() seems to have problem
# i.e. it produces the trees that the treecmp tool cannot compute the MS distance (need further exploration)
'''
treeName,treeExt = splitext(outtrees)
for i,RS in enumerate(removing_sets):
    trees_shrinked = deepcopy(trees)
    for t,tree in enumerate(trees_shrinked):
        filt = lambda node: False if (node.taxon is not None and node.taxon.label in RS[t]) else True 
        tree.filter_leaf_nodes(filt,update_bipartitions=True)
    trees_shrinked.write_to_path(outdir + "/" + treeName + "_" + quantiles[i] + treeExt,'newick')   
#print(removing_sets)
''' 

# prune trees according to the removing sets. 
# calling nw_prune here as a (terrible) temporary sollution, because Dendropy's filter_leaf_nodes() seems to have problems
fName,ext = splitext(outtrees)
for i,RS in enumerate(removing_sets):
    outfile = outdir + "/" + fName + "_RS_" + quantiles[i]
    with open(outfile,'w') as f:
        for item in RS:
            for s in item:
                f.write(s + "\t")
            f.write("\n")
    trees_shrinked =  outdir + "/" + fName + "_" + quantiles[i] + ext
    call(["prune_trees.sh",intrees,outfile,trees_shrinked]) 
