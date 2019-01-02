#! /usr/bin/env python

from treeshrink.filter_lib import filter_branch
from sys import argv
from dendropy import Tree

from os.path import splitext
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-i','--input',required=True,help="input file")
parser.add_argument('-o','--outfile',required=True,help="output file")
parser.add_argument('-u','--unit',required=False,help="unit-length for unit-based filter")
parser.add_argument('-l','--lowthres',required=False,help="low threshold")
parser.add_argument('-g','--highthres',required=False,help="high threshold")
parser.add_argument('-f','--factor',required=False,help="factor")

args = vars(parser.parse_args())

infile = args['input']
outfile = args['outfile']
a_tree = Tree.get_from_path(infile,"newick",preserve_underscores=True)
unit=args['unit'] if args['unit'] else None
low = float(args['lowthres']) if args['lowthres'] else 0
high = float(args['highthres']) if args['highthres'] else 1
factor = float(args['factor']) if args['factor'] else 1

filter_branch(a_tree,unit_length=args['unit'],low_percentile=low,high_percentile=high,factor=factor)

a_tree.write_to_path(outfile,"newick")
