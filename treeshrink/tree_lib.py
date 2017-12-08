import sys
from dendropy import Tree


def prune_node(T,node):
    if node is not T.seed_node:
        p = node.parent_node
        p.remove_child(node)
        if p.num_child_nodes() == 1:
            v = p.child_nodes()[0]
            p.remove_child(v)
            if p is T.seed_node:
                T.seed_node = v
            #    p.remove_child(v)
            else:
                u = p.parent_node
                l = p.edge_length + v.edge_length
                u.remove_child(p)
                u.add_child(v)
                v.edge_length = l



def prune_tree(T,RS):
# prune the taxa in the removing set RS from tree T
    L = list(T.leaf_node_iter())
    for leaf in L:
        if leaf.taxon.label in RS:
            prune_node(T,leaf)

def get_taxa(tree_file,scheme='newick'):
	a_tree = Tree.get_from_path(tree_file,scheme,preserve_underscores=True)
	return [leaf.taxon.label for leaf in a_tree.leaf_nodes()]

def report_taxa(tree_file,scheme='newick',listing=True,counting=True):
	a_tree = Tree()
	a_tree.read_from_path(tree_file,scheme)
	if listing:
		for leaf in a_tree.leaf_nodes():
			print(leaf.taxon.label)
	if counting:
		print('Taxa #: ' + str(len(a_tree.leaf_nodes())))

def tree_as_newick(a_tree,outfile=None,append=False):
# dendropy's method to write newick seems to have problem ...
	if outfile:
		outstream = open(outfile,'a') if append else open(outfile,'w')
	else:
		outstream = sys.stdout

	__write_newick(a_tree.seed_node,outstream)

	outstream.write(";\n")
	if outfile:
		outstream.close()	

def __write_newick(node,outstream):
	if node.is_leaf():
			if node.taxon:
				outstream.write(node.taxon.label)
			else:
				outstream.write(str(node.label))
	else:
		outstream.write('(')
		is_first_child = True
		for child in node.child_node_iter():
			if is_first_child:
				is_first_child = False
			else:
				outstream.write(',')
			__write_newick(child,outstream)
		outstream.write(')')
	if not node.is_leaf() and node.label is not None:
			outstream.write(str(node.label))

	if not node.edge_length is None:
		outstream.write(":"+str(node.edge_length))
