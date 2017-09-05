#! /usr/bin/env python

import random
from sys import argv

infile=argv[1]
k=int(argv[2])
L = []


with open(infile,"r") as f:
    for line in f:
        L.append(line.rstrip())
    L1 = random.sample(L,k)
    for item in L1:
        print(item)

