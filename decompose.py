#! /usr/bin/env python

from treeshrink.decompose_lib import decompose
from treeshrink.alignment import Alignment
from treeshrink.sequence_lib import sample_from_list
from treeswift import *
import argparse
from os.path import basename, dirname, splitext,realpath,join,normpath,isdir,isfile,exists
from os import mkdir,getcwd,rmdir,listdir
from shutil import copyfile
#from brlen_filter import find_thresholds

parser = argparse.ArgumentParser()
parser.add_argument("-i","--indir",required=False,help="The parent input directory where the trees (and alignments) can be found")
parser.add_argument("--range",required=False,help="Only decompose genes in this range. Default: decompose all")
parser.add_argument("--genes",required=False,help="Only decompose genes in this list. Default: decompose all")
parser.add_argument("-t","--tree",required=False,default="input.tree",help="The name of the input tree/trees. Each subdirectory under must contain a tree with this name. Default: input.tree")
parser.add_argument("-a","--alignment",required=False,help="The name of the input alignment. Each subdirectory must contain an alignment with this name. Default: None")
parser.add_argument("--minSize",required=False,default=None,help="The minimum number of taxa in each subtree. Default: square root of tree size")
parser.add_argument("--minBranch",required=False,type=float,default=1.0,help="The minimum branch length that could be cut. Default: 1.0")
parser.add_argument("-o","--outdir",required=False,help="Output directory.")
parser.add_argument("-T","--outTree",required=False,help="The name of the output tree/trees. Each subdirectory under outdir will contain a tree with this name. Default: the same name as input.")
parser.add_argument("-A","--outAlignment",required=False,help="The name of the output alignment(s). Each subdirectory will contain this/these alignment(s). If more than one alignment is supplied per gene, please put in quotes. Default: the same name(s) as input.")
    
args = vars(parser.parse_args())

treename = splitext(args["tree"])[0]
min_nleaf = int(args["minSize"]) if args["minSize"] else None
min_brlen = args["minBranch"]
idx_range = tuple(float(x) for x in args["range"].split()) if args["range"] else (0,float("inf"))
print(idx_range)

if args["outdir"]:
    outdir = args["outdir"]
elif args["indir"]:
    outdir = args["indir"] + "_decomposed"  
else:    
    outdir = treename + "_decomposed"
if not exists(outdir):
    mkdir(outdir)
if not exists(normpath(join(outdir,'annotated_trees'))):
    mkdir(normpath(join(outdir,'annotated_trees')))

print("Reading input trees ...")
if args["indir"]:
    if not args["genes"]:
        subdirs = listdir(args["indir"]) 
    else:
        subdirs = args["genes"].strip().split()
    tree_strs = []
    gene_names = []
    aln_files = []
    count = 1
    #for d in listdir(args["indir"]):
    for d in subdirs:
        count += 1
        if exists(normpath(join(args["indir"],d,args["tree"]))):
            gene_names.append(d)
            tree_strs.append(open(normpath(join(args["indir"],d,args["tree"])),'r').read())
            if args["alignment"]:
                A = []
                for a in args["alignment"].strip().split():
                    aln_file = normpath(join(args["indir"],d,a))
                    if exists(aln_file):
                        A.append(aln_file)
                aln_files.append(A)
else:
    tree_strs = open(args["tree"],'r').readlines()    
    gene_names = ["gene_" + str(i+1).rjust(4,'0') for i in range(len(tree_strs))]
    aln_files = None

ntrees = len(tree_strs)
outTree = args["outTree"] if args["outTree"] else args["tree"]
outAlns = args["outAlignment"] if args["outAlignment"] else args["alignment"]


#thresholds = find_thresholds(tree_strs,p=0.01)
print("Decomposing ...")
for i in range(ntrees):  
    if i+1 < idx_range[0]:
        continue
    if i+1 > idx_range[1]:
        break
    print(i+1,gene_names[i])    
    tree = read_tree_newick(tree_strs[i])
    gene = gene_names[i]
    #thres = thresholds[i+1] if (i+1) in thresholds else 1000 
    #if exists(normpath(join(outdir,gene + "_decomposed_1"))):
    #    continue
    decomposed_trees, ann_tree = decompose(tree,min_nleaf=min_nleaf,min_brlen=min_brlen)
    
    # produce annotated trees
    if len(decomposed_trees) > 1:
        with open(normpath(join(outdir,'annotated_trees',gene + ".tre")),'w') as fout_ann:
            fout_ann.write("#NEXUS\n")
            fout_ann.write("begin trees;\n")
            fout_ann.write("tree 1 = ")
            fout_ann.write(ann_tree)
            fout_ann.write("\n")
            fout_ann.write("end;")

    # output trees
    for t,treestr in enumerate(decomposed_trees):
        od = normpath(join(outdir,gene + "_decomposed_" + str(t+1)))
        mkdir(od)
        with open(normpath(join(od,outTree)),'w') as f:
            f.write(treestr)
        
    # output alignments
    if len(decomposed_trees) == 1:
        od = normpath(join(outdir,gene + "_decomposed_1"))
        for alnfile,alnName in zip(aln_files[i],outAlns.strip().split()):
            alnfile_out = (normpath(join(od,alnName)))
            copyfile(alnfile,alnfile_out)
    else:        
        for t,treestr in enumerate(decomposed_trees):
            od = normpath(join(outdir,gene + "_decomposed_" + str(t+1)))
            T = read_tree_newick(treestr)    
            for alnfile,alnName in zip(aln_files[i],outAlns.strip().split()):
                alnfile_out = (normpath(join(od,alnName)))
                sample_from_list(alnfile,[node.label for node in T.traverse_leaves()],alnfile_out,store_index_file=True,renew_index_file=False)
                #aln = Alignment()
                #aln.read_filepath(alnfile)
                #A = aln.sub_alignment([node.label for node in T.traverse_leaves()])
                #A.write_filepath(normpath(join(od,alnName))) 
