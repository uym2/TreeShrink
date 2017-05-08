import numpy as np
import operator
from dendropy import Tree
from copy import deepcopy

def filter_branch(a_tree,percentile=0.5,factor=1):
    a_tree.deroot()
    branch_list = list_branch(a_tree,sort=(percentile>=0))
    d = compute_diameter(a_tree,branch_list,percentile=percentile)
    thres = d*factor
    #print(thres)
    count_leaves(a_tree)
    for br in branch_list:
        #print(br.length)
        if br.length > thres:
            remove_branch(a_tree,br)

def count_leaves(a_tree):
    for node in a_tree.postorder_node_iter():
        if node.is_leaf():
            node.nleaf = 0
        else:
            node.nleaf = sum([ch.nleaf for ch in node.child_node_iter()])


def remove_branch(a_tree,br):
            p = br.tail_node
            c = br.head_node
            if c.nleaf > a_tree.seed_node.nleaf/2:
                p.remove_child(c)
                a_tree.seed_node = c
            elif p is not None:
                p.remove_child(c,suppress_unifurcations=True)
                while p and p.num_child_nodes() == 0:
                    c = p
                    p = p.parent_node
                    if p is not None:
                        p.remove_child(c,supress_unifurcations=True)


def list_branch(a_tree,sort=True):
    branch_list = []
    for br in a_tree.preorder_edge_iter():
        if br.tail_node is not None:
            branch_list.append(br)
    if sort:
        branch_list.sort(key=lambda x: x.length,reverse=True)
    
    return branch_list


def compute_diameter(a_tree,branch_list,percentile=0.5):
    # percentile < 0: take mean instead of percentile. Default is the median
    if percentile >= 0:
        unit_length = branch_list[int(len(branch_list)*percentile)].length
    else:
        unit_length = 0
        for br in branch_list:
            unit_length += br.length
            unit_length /= len(branch_list)

    max_br_distance = 0

    for node in a_tree.postorder_node_iter():
        if node.is_leaf():
            node.max_br_below = 0
        else:
            children = node.child_nodes()
            max1 = max(children[0].max_br_below+1,children[1].max_br_below+1)
            max2 = min(children[0].max_br_below+1,children[1].max_br_below+1)
            i = 2
            while i < len(children):
                if children[i].max_br_below+1 > max1:
                    max1 = children[i].max_br_below+1
                elif children[i].max_br_below+1 > max2:
                    max2 = children[i].max_br_below+1
                i += 1

            max_br_distance = max(max_br_distance,max1 + max2)
            node.max_br_below = max1
    d =  max_br_distance*unit_length

    #print("Estimated diameter: ", d)

    return d

