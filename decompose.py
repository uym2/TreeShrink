#! /usr/bin/env python

from treeshrink.decompose_lib import decompose
from treeshrink.alignment import Alignment
from treeswift import *
import argparse
from os.path import basename, dirname, splitext,realpath,join,normpath,isdir,isfile,exists
from os import mkdir,getcwd,rmdir,listdir
from brlen_filter import find_thresholds

parser = argparse.ArgumentParser()
parser.add_argument("-i","--indir",required=False,help="The parent input directory where the trees (and alignments) can be found")
parser.add_argument("-t","--tree",required=False,default="input.tree",help="The name of the input tree/trees. Each subdirectory under must contain a tree with this name. Default: input.tree")
parser.add_argument("-a","--alignment",required=False,help="The name of the input alignment. Each subdirectory must contain an alignment with this name. Default: None")
parser.add_argument("--minSize",required=False,type=int,default=20,help="The minimum number of taxa in each subtree. Default: 20")
parser.add_argument("--minBranch",required=False,type=float,default=1.0,help="The minimum branch length that could be cut. Default: 1.0")
parser.add_argument("-o","--outdir",required=False,help="Output directory.")
    
args = vars(parser.parse_args())

treename = splitext(args["tree"])[0]
min_nleaf = args["minSize"]
min_brlen = args["minBranch"]

if args["outdir"]:
    outdir = args["outdir"]
elif args["indir"]:
    outdir = args["indir"] + "_decomposed"  
else:    
    outdir = treename + "_decomposed"
mkdir(outdir)
mkdir(normpath(join(outdir,'annotated_trees')))

if args["indir"]:
    subdirs = [d for d in listdir(args["indir"]) if exists(normpath(join(args["indir"],d,args["tree"])))]
    tree_strs = []
    gene_names = []
    aln_files = []
    for d in listdir(args["indir"]):
        if exists(normpath(join(args["indir"],d,args["tree"]))):
            gene_names.append(d)
            tree_strs.append(open(normpath(join(args["indir"],d,args["tree"])),'r').read())
            if args["alignment"]:
                aln_file = normpath(join(args["indir"],d,args["alignment"]))
                aln_files.append(aln_file if exists(aln_file) else None)
else:
    tree_strs = open(args["tree"],'r').readlines()    
    gene_names = ["gene_" + str(i+1).rjust(4,'0') for i in range(len(tree_strs))]
    aln_files = None

ntrees = len(tree_strs)


thresholds = find_thresholds(tree_strs,p=0.01)

for i in range(ntrees):    
    tree = read_tree_newick(tree_strs[i])
    gene = gene_names[i]
    thres = thresholds[i+1] if (i+1) in thresholds else 1000
    decomposed_trees, ann_tree = decompose(tree,min_nleaf=min_nleaf,min_brlen=thres)
    
    # produce annotated trees
    if len(decomposed_trees) > 1:
        with open(normpath(join(outdir,'annotated_trees',gene + ".tre")),'w') as fout_ann:
            fout_ann.write("#NEXUS\n")
            fout_ann.write("begin trees;\n")
            fout_ann.write("tree 1 = ")
            fout_ann.write(ann_tree)
            fout_ann.write("\n")
            fout_ann.write("end;")
    
    alnfile = aln_files[i] if aln_files else None
    if alnfile:
        aln = Alignment()
        aln.read_filepath(alnfile)
    for t,T in enumerate(decomposed_trees):
        od = normpath(join(outdir,gene + "_decomposed_" + str(t+1)))
        mkdir(od)
        T.write_tree_newick(normpath(join(od,"tree.tre")))
        if alnfile:
            A = aln.sub_alignment([node.label for node in T.traverse_leaves()])
            A.write_filepath(normpath(join(od,"aln.fasta")))
