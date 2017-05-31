# ! /usr/bin/env python

from optimal_filter_lib import TreeFilter
from sys import argv
from math import sqrt
from subprocess import check_output,call

intree = argv[1]
outtree = argv[2]

myfilter = TreeFilter(tree_file=intree)
myfilter.optFilter()

datafile = check_output(["mktemp"]).rstrip()
fout = open(datafile,"w")

for i in range(1,len(myfilter.min_diams)):
   fout.write(str(myfilter.min_diams[i-1]-myfilter.min_diams[i]) + "\n")

fout.close()

opt_k = int(check_output(["Rscript","find_d.R",datafile,"0.05"])[4:].rstrip())

fTree = myfilter.filterOut(d=opt_k)
fTree.write_to_path(outtree,"newick")

call(["rm",datafile])
