# ! /usr/bin/env python

from optimal_filter_lib import TreeFilter
from sys import argv
from Tree_extend import MV00_Tree
from math import sqrt

intree = argv[1]
#outtree = argv[2]

#MVTree = MV00_Tree(tree_file=tree)
#MVTree.find_root()

#std = sqrt(MVTree.minVAR)

myfilter = TreeFilter(tree_file=intree)
myfilter.optFilter(d=30)

#myfilter.list_removals(d=2)

#fTree = myfilter.filterOut(d=172)

#fTree.write_to_path(outtree,"newick")

#print(myfilter.best_entries[1].removed)
for i in range(1,len(myfilter.min_diams)):
   print(myfilter.min_diams[i-1]/ myfilter.min_diams[i])
#   print(myfilter.min_diams[i-1] - myfilter.min_diams[i])
   #print(myfilter.min_diams[i-1] - myfilter.min_diams[i])/myfilter.min_diams[i]
   #print(myfilter.min_diams[i-1] - myfilter.min_diams[i])/myfilter.min_diams[i]/std
#   print(myfilter.min_diams[i-1])


