#! /usr/bin/env python

from filter_lib import filter_branch
from sys import argv
from dendropy import Tree

from os.path import splitext
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-i','--input',required=True,help="input file")
parser.add_argument('-o','--outfile',required=True,help="output file")
parser.add_argument('-p','--percentile',required=False,help="default is 0.5. Fill in -1 to use the mean instead of percentile")
parser.add_argument('-factor','--factor',required=False,help="default is 1")

args = vars(parser.parse_args())

infile = args['input']
outfile = args['outfile']
a_tree = Tree.get_from_path(infile,"newick")
p = float(args['percentile']) if args['percentile'] else 0.5
f = float(args['factor']) if args['factor'] else 1

filter_branch(a_tree,percentile=p,factor=f)

a_tree.write_to_path(outfile,"newick")
