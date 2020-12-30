#! /usr/bin/env python
from numpy import median, mean, std, quantile
from scipy.stats import norm
from math import log10
from sys import argv
from treeswift import *

def collect_features(tree,minLen=0.001):
    # collect branch lengths of the input tree
    # put them into a list of tuples (idx,brlen)
    brlens = [node.edge_length for node in tree.traverse_preorder() if not node.is_root() and not node.is_leaf() and node.edge_length > minLen]
    brlens_log = [log10(x) for x in brlens]
    q3 = 10**find_cutoff(brlens_log,0.25)
    q1 = 10**find_cutoff(brlens_log,0.75)
    m = mean([x for x in brlens if x > q1 and x < q3])
    normBrlens = [(br/m-1) for br in brlens]
    return brlens,m,normBrlens

def find_cutoff(numlist,p=0.01):
    # find the outliers of a list of numbers
    # fit here a Gaussian model and use p-value to determine outliers
    mu = mean(numlist)
    sigma = std(numlist)
    N = norm(loc=mu,scale=sigma)
    
    return N.ppf(1-p)

def find_thresholds(treeStrs,p=0.01):
    # collect features of all trees
    # identify by tuples (treeIdx,brlens,normBrlens)
    features = []
    for i,s in enumerate(treeStrs):
        tree = read_tree_newick(s)
        brlens, m, normBrlens = collect_features(tree,minLen=0.001)
        features += [(i+1,br,normbr) for (br,normbr) in zip(brlens,normBrlens) if normbr>0] # only include the features of branches above median
        
    q = find_cutoff([log10(x) for (_,_,x) in features],p=p)
    cutoff = 10**q
    thresholds = {}
    
    for i,br,normbr in features:
        if normbr > cutoff:
            if i not in thresholds:
                thresholds[i] = br
            elif br < thresholds[i]:
                thresholds[i] = br

    return thresholds            

def main():
    inFile = argv[1]
    outFile = argv[2]

    with open(inFile,'r') as fin:
        treeStrs = fin.read().split()

    thresholds = find_thresholds(treeStrs,p=0.01)
    with open(outFile,'w') as fout:
        for i in thresholds:
            fout.write(str(i) + " " + str(thresholds[i]) + "\n")

if __name__ == "__main__":
    main()
