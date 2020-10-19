#! /usr/bin/env python

import argparse
from os import mkdir,listdir
from shutil import copyfile
from os.path import exists,splitext,normpath,join

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--indir",required=True,help="Input folder which contains all trees and alignments")
    parser.add_argument("-o","--outdir",required=True,help="Output folder to be created, must be non-exist")
    parser.add_argument("-t","--treeExt",required=False,default="tree",help="Tree extension. Default: .tre")
    parser.add_argument("-a","--alnExt",required=False,default="fasta",help="Alignment extension. Default: .fasta")
    
    args = vars(parser.parse_args())
    
    if exists(args["outdir"]):
        print("Output folder already exists. Abort now!")
        exit(0)

    mkdir(args["outdir"])

    for f in listdir(args["indir"]):
        name,ext = splitext(f)
        if ext != args["treeExt"] and ext != args["alnExt"]:
            continue    
        if not exists(normpath(join(args["outdir"],name))):
            mkdir(normpath(join(args["outdir"],name)))
        if ext == args["treeExt"]:
            copyfile(normpath(join(args["indir"],f)),normpath(join(args["outdir"],name,"input.tree")))
        else:
            copyfile(normpath(join(args["indir"],f)),normpath(join(args["outdir"],name,"input.fasta")))        

if __name__ == "__main__":
    main()
