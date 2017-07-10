#! /usr/bin/env python

from sys import argv

infile = argv[1]
outfile = argv[2]

fout = open(outfile,"w")

with open(infile,"r") as f:
    flag = False
    for line in f:
        if line.split()[0] == "Tree":
            if flag:
                fout.write(str(count)+"\n")
            else:
                flag = True
            count = 0
        else:
            count = count+1
    # write the last time
    fout.write(str(count)+"\n") 
