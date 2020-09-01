#! /usr/bin/env python

from treeshrink.decompose_lib import decompose
from treeshrink.alignment import Alignment
from treeswift import *
import argparse
from os.path import basename, dirname, splitext,realpath,join,normpath,isdir,isfile,exists
from os import mkdir,getcwd,rmdir,listdir

parser = argparse.ArgumentParser()
parser.add_argument("-i","--indir",required=False,help="The parent input directory where the trees (and alignments) can be found")
parser.add_argument("-t","--tree",required=False,default="input.tree",help="The name of the input tree/trees. Each subdirectory under must contain a tree with this name. Default: input.tree")
parser.add_argument("-a","--alignment",required=False,help="The name of the input alignment. Each subdirectory must contain an alignment with this name. Default: input.fasta")
parser.add_argument("-o","--outdir",required=True,help="Output directory.")
    
args = vars(parser.parse_args())


treename = splitext(args["tree"])[0]
subdirs = [d for d in listdir(args["indir"]) if exists(normpath(join(args["indir"],d,args["tree"])))] #if args["tree"] else "input.tre")))]
outdir = args["outdir"]
mkdir(outdir)

for d in subdirs:
    treefile = normpath(join(args["indir"],d,args["tree"]))
    alnfile = normpath(join(args["indir"],d,args["alignment"]))
    aln = Alignment()
    aln.read_filepath(alnfile)
    tree = read_tree_newick(treefile)
    decomposed_trees = decompose(tree,min_nleaf=20,min_brlen=1.5)
    for t,T in enumerate(decomposed_trees):
        od = normpath(join(args["outdir"],d + "_decomposed_" + str(t+1)))
        mkdir(od)
        T.write_tree_newick(normpath(join(od,"tree.nwk")))
        A = aln.sub_alignment([node.label for node in T.traverse_leaves()])
        A.write_filepath(normpath(join(od,"aln.fasta")))
