from treeswift import *
from heapq import *
from treeshrink.alignment import Alignment
from copy import deepcopy

def bisect(tree,min_nleaf=50,min_brlen=1.5):
    # assume there is no unifurcation
    B = []
    R = tree.root.child_nodes()
    if len(R) == 2:
        el = R[0].get_edge_length() + R[1].get_edge_length()
        B.append((el,R[0]))
        B.append((el,R[1]))
    else:
        B += [(node.get_edge_length(),node) for node in R]

    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.nleaf = 1
        else:
            node.nleaf = sum(c.nleaf for c in node.child_nodes())    
        if not node.is_root(): # and not node in R:
            B.append((node.get_edge_length(),node))
    
    # sort by edge length
    B.sort()
    # the total number of leaf
    n_total = tree.root.nleaf
    # find a branch to bisect
    output = {'T1':None,'T2':None,'cut_branch':None,'success':False}
    while B:
        br,node = B.pop()
        nIn,nOut = node.nleaf,n_total-node.nleaf
        if br < min_brlen:
            break
        if nIn >= min_nleaf and nOut >= min_nleaf:
            # bisect on this branch
            p = node.get_parent()
            p.remove_child(node)
            new_tree = Tree()
            new_tree.root = node
            tree.suppress_unifurcations()
            new_tree.suppress_unifurcations()
            output['T1'] = tree
            output['T2'] = new_tree
            output['cut_branch'] = br
            output['success'] = True
            output['cut_pos'] = node.id # mark that we cut on the branch above this node
            return output
    return output        

def decompose(tree,min_nleaf=20,min_brlen=1):
    # give each node an ID
    ID = 0
    for node in tree.traverse_preorder():
        node.id = ID
        ID += 1
    backup_tree = deepcopy(tree)    
    
    stack = [tree]
    tree_list = []
    cut_pos = []
    
    while stack:
        T = stack.pop()
        output = bisect(T,min_nleaf=min_nleaf,min_brlen=min_brlen)
        if output['success']:
            stack += [output['T1'],output['T2']]
            cut_pos.append(output['cut_pos'])
        else:
            tree_list.append(T)    
    
    # annotate
    for node in backup_tree.traverse_preorder():
        if node.id in cut_pos:
            lb = node.label if node.label is not None else ''
            node.label = lb + '[&!color=#0000ff]' 
                        
    return tree_list, backup_tree.newick()

def main():                
    tree = read_tree_newick("test.tre")
    algn = Alignment()
    algn.read_filepath("test.fasta")

    tree_list = decompose(tree)

    with open("decomposed.trees",'w') as fout:
        for t,T in enumerate(tree_list):
            fout.write(T.newick() + "\n")
            L = [node.label for node in T.traverse_leaves()]
            sub_aln = algn.sub_alignment(L)
            sub_aln.write_filepath("decomposed_" + str(t+1) + ".fasta")

if __name__ == "__main__":
    main()        
