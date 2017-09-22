#! /usr/bin/env python

from sys import argv

filein = argv[1]
fileout = argv[2]


d = {}
with open(filein,'r') as fin:
    for line in fin:
        method, condition, occ, dMS = line.rstrip().split()
        key = method + "#" + condition
        occ = float(occ)
        dMS = float(dMS)
        if key not in d:
            d[key] = (occ,dMS,1)
        else:
            prev_occ, prev_dMS, prev_count = d[key]
            d[key] = (occ+prev_occ, dMS+prev_dMS, prev_count+1)
            
with open(fileout,'w') as fout:
    for key in d:
        occ,dMS,count = d[key]
        method,condition = key.split("#")
        occ = occ/count
        dMS = dMS/count
        fout.write(method + " " + condition + " " + str(occ) + " " + str(dMS) + "\n")
                            
