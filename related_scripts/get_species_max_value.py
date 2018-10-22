#! /usr/bin/env python

from sys import argv

filein = argv[1]
fileout = argv[2]

with open(filein,'r') as f:
	d = {}
	for line in f:
		L = line.split()
		r = float(L[0])
		for s in L[1:]:
			if not s in d:
				d[s] = r
			else:
				d[s] = max(r,d[s])

with open(fileout,'w') as f:
	for s in d:
		f.write(s + " " + str(d[s]) + "\n")
		
