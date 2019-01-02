#! /usr/bin/env python

# usage: python filter_branch.py <tree_file>

from treeshrink.Tree_extend import MVR_Tree
from sys import argv
from os.path import splitext

tree_file = argv[1]
outfile = argv[2]

base_name,ext = splitext(tree_file)

a_tree = MVR_Tree(tree_file=tree_file)
a_tree.filter_branch()

a_tree.tree_as_newick(outfile=outfile)


