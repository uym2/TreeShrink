#! /usr/bin/env python

from filter_lib import filter_branch
from sys import argv
from dendropy import Tree

infile = argv[1]
outfile = argv[2]
a_tree = Tree.get_from_path(infile,"newick")

filter_branch(a_tree)

a_tree.write_to_path(outfile,"newick")
