from dendropy import Tree
from math import sqrt
from Queue import Queue
from copy import deepcopy

class TreeInduced:
    def __init__(self,bestLCA=None):
        self.bestLCA = None
        self.records = {}

class Entry:
    def __init__(self,bestLCA=None):
        self.level = 0
        self.info = TreeInduced(bestLCA=bestLCA)
        self.backtrack = None

class TreeFilter:       
    def __updateNode__(self,node,records=None):
        if node.is_leaf():
            records[node] = (node,node,0,0)
            return

        if records is None:
            records = self.records

        max1 = -1
        max2 = -1
        node1 = None
        node2 = None

        for ch in node.child_node_iter():
               if records[ch] is None:
                    # this child node was removed!
                    continue
               l = records[ch][2] + ch.edge_length
               if l > max1:
                   max1 = l
                   node1 = ch
               elif l > max2:
                   max2 = l
                   node2 = ch

        MAX = max(max1,max1+max2)
        records[node] = (node1,node2,max1,MAX) if node1 else None

   def __substitute_anchor__(self,entry,anchor,inherit_info=True):
        anchor1 = self.__get_anchor1__(entry)
        anchor2 = self.__get_anchor2__(entry)

        if anchor is anchor1:
            retained_anchor = anchor2
        elif anchor is anchor2:
            retained_anchor = anchor1
        else:
            print("ERROR : TreeFilter.__substitute_anchor__: The request anchor to be removed is not one of the two anchors!")
            return None

        new_entry = Entry()
        if inherit_info:
            new_entry.info = entry.info
        else:
            new_entry.info = deepcopy(entry.info)
        
        # update node records from anchor to root
        new_entry.info.records[anchor] = None 
        curr_node = anchor.parent_node
        new_best_LCA = curr_node
        new_best_record = self.__lookup__(new_entry,new_best_LCA)

        while curr_node:
            self.__updateNode__(curr_node,records=new_entry.info.records)
            node_record = self.__lookup__(new_entry,curr_node)
            if node_record is not None and (new_best_record is None or node_record[3] > new_best_record[3]): 
                new_best_LCA = curr_node
                new_best_record = self.__lookup__(new_entry,new_best_LCA)
            curr_node = curr_node.parent_node

        # lookup for the new bestLCA: it must be on the path from the retained anchor to the current bestLCA of the entry
        curr_node = retained_anchor.parent_node
        node_record = self.__lookup__(new_entry,curr_node)

        while curr_node is not entry.info.bestLCA:
            if node_record is not None and (new_best_record is None or node_record[3] > new_best_record[3]):
                new_best_LCA = curr_node
                new_best_record = self.__lookup__(new_entry,new_best_LCA)
            curr_node = curr_node.parent_node        
            node_record = self.__lookup__(new_entry,curr_node)
         
        if node_record is not None and (new_best_record is None or node_record[3] > new_best_record[3]):
            new_entry.info.bestLCA = curr_node
        else:
            new_entry.info.bestLCA = new_best_LCA

        # backtrack
        new_entry.backtrack.entry = entry
        new_entry.backtrack.removed = anchor
        new_entry.backtrack.retained = retained_anchor

        return new_entry

   def __substitute_anchor1__(self,entry,inherit_info=False):
        return self.__substitute_anchor__(self,entry,self.__get_anchor1__(entry),inherit_info=inherit_info)

   def __substitute_anchor2__(self,entry,inherit_info=True):
        return self.__substitute_anchor__(self,entry,self.__get_anchor2__(entry),inherit_info=inherit_info)
   
   def __get_diam__(self,entry):
        try:
            d = entry.info.records[entry.info.bestLCA][3]
        except:
            d = self.records[entry.info.bestLCA][3]
        return d

    def __get_anchor1__(self,entry):
        try:
            node = entry.info.records[entry.info.bestLCA][0]
        except:
            node = self.records[entry.info.bestLCA][0]
        return node


    def __get_anchor2__(self,entry):
        try:
            node = entry.info.records[entry.info.bestLCA][1]
        except:
            node = self.records[entry.info.bestLCA][1]
        return node
    
    def __lookup__(self,entry,node):
        try:
            node_record = entry.info.records[node]
        except:
            node_record = self.records[node]

        return node_record
    
    def __init__(self, ddpTree = None, tree_file = None, schema = "newick"):
        if tree_file:
            self.ddpTree = Tree.get_from_path(tree_file,schema)
        else:  
            self.ddpTree = ddpTree

        self.bestLCA = None
        self.nleaf = 0
        self.records = {}

        diam = -1

        for node in self.ddpTree.postorder_node_iter():
               if node.is_leaf():
                   self.nleaf += 1
                   continue
               self.__updateNode__(self,node)
               if self.records[node][3] > diam:
                   diam = self.records[node][3]
                   self.bestLCA = node

        first_entry = Entry(bestLCA=self.bestLCA)
        self.myQueue = [first_entry]
        self.best_entries = []
        self.min_diams = []

   def optFilter(self,d=None):
        if d is None:
            d = 2*int(sqrt(self.nleaf))

        curr_level = -1
        min_diam = self.__get_diam__(first_entry)
        best_entry = first_entry

        while 1:
            if self.myQueue.empty():
                break
            curr_entry = myQueue.get()
            diam = self.__get_diam(curr_entry)

            if curr_entry.level != curr_level:
                # reached a next level
                curr_level = curr_entry.level
                self.best_entries.append(best_entry)
                self.min_diams.append(min_diam)
                min_diam = diam
                best_entry = curr_entry
                if curr_level < d:
                    # add 2 entries
                    self.myQueue.put(self.__substitute_anchor1__(curr_entry))
                    self.myQueue.put(self.__substitute_anchor2__(curr_entry))
            else:
                if diam < min_diam:
                    min_diam  = diam
                    best_entry = curr_entry
                if curr_level < d:
                    # add 1 entry
                    anchor1 = self.__get_anchor1__(curr_entry)
                    anchor2 = self.__get_anchor2__(curr_entry)
                    if anchor1 is curr_entry.backtrack.retained:
                        anchor = anchor1
                    elif anchor2 is curr_entry.backtrack.retained:
                        anchor = anchor2
                    else:
                        print("FilterTree.optFilter: anchors inconsistent between child and parent!")
                        return False
                    self.myQueue.put(self.__substitute_anchor__(curr_entry),anchor)

            curr_entry.info = None
