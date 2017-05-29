# ! /usr/bin/env python

from optimal_filter_lib import TreeFilter
from sys import argv
from Tree_extend import MV00_Tree
from math import sqrt

tree = argv[1]


#MVTree = MV00_Tree(tree_file=tree)
#MVTree.find_root()

#std = sqrt(MVTree.minVAR)

myfilter = TreeFilter(tree_file=tree)
myfilter.optFilter(d=60)

#print(myfilter.best_entries[1].removed)
for i in range(1,len(myfilter.min_diams)):
   #print(myfilter.min_diams[i-1]/ myfilter.min_diams[i])
   #print(myfilter.min_diams[i-1] - myfilter.min_diams[i])
   print(myfilter.min_diams[i-1] - myfilter.min_diams[i])/myfilter.min_diams[i]
   #print(myfilter.min_diams[i-1] - myfilter.min_diams[i])/myfilter.min_diams[i]/std
   #print(myfilter.min_diams[i-1])


