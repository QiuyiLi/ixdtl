import skbio
from collections import defaultdict
from statistics import mean
from .tree_table import *


class SpeciesTree:
    def __init__(self):
        self.__skbioTree = None
        self.__treeTable = None
        self.__parameter = None

    @property
    def skbioTree(self):
        return self.__skbioTree

    @property
    def treeTable(self):
        return self.__treeTable

    @property
    def parameter(self):
        return self.__parameter
        
    @parameter.setter
    def parameter(self, parameter):
        self.__parameter = parameter

    def getNodeById(self, id):
        return self.__treeTable.getEntryById(id)

    def getNodes(self):
        return self.__treeTable.table

    def getNodeByName(self, name):
        return self.__treeTable.getEntryByName(name)

    def getRoot(self):
        return self.__treeTable.getRoot()
    
    def readFromNewickFile(self, path):
        self.__treeTable = TreeTable()
        self.__skbioTree = self.__treeTable.createFromNewickFile(path)

        print(self.__treeTable)
        print(self.getRoot())

    def coalescent(self, distance_above_root):
        """
        the main multi-species coalecent function
        """
        nodes = self.getNodes()
        root = self.getRoot()
        coalescent_process = defaultdict(list)
        clade_set_into_root = None

        old_leaves = [node.id for node in nodes if not node.children]      # leaves of the given species tree
        new_leaves = []     # leaves set will be updated in the loop
        clade_set = {}      # set of extant species that an ancestral gene will finally be fixed in
        labelled = {}       # avoid doing repeated coalescence
        for node in nodes:
            labelled[node.id] = False
            clade_set[node.id] = [str(node.id) + '*'] if not node.children else []

        while True:
            for leaf in old_leaves:
                if (leaf == root.id):
                    clade_set_into_root = self.__coalescentRecurse(id=root.id, 
                                                                   distance=distance_above_root, 
                                                                   clade_set=clade_set,
                                                                   coalescent_process=coalescent_process)
                    break
                else:
                    parent = self.getNodeById(leaf).parent
                    children = self.getNodeById(parent).children
                    if (labelled[leaf]):
                        continue
                    labelled[leaf] = True
                    if (len(clade_set[children[0]]) != 0 
                        and len(clade_set[children[1]]) != 0):
                        self.__coalescentRecurse(id=children[0], 
                                                 distance=self.getNodeById(children[0]).distanceToParent,
                                                 clade_set=clade_set,
                                                 coalescent_process=coalescent_process)
                        self.__coalescentRecurse(id=children[1], 
                                                 distance=self.getNodeById(children[1]).distanceToParent,
                                                 clade_set=clade_set,
                                                 coalescent_process=coalescent_process)
                        # the clade set of the parent before coalescence is the union of the clade set of its children after coalescence                        
                        clade_set[parent] = list(set().union(clade_set[children[0]], clade_set[children[1]]))    
                        if (len(new_leaves) > 0):
                            new_leaves = [e for e in new_leaves if e != children[0] and e != children[1]]
                        new_leaves.append(parent)
                    else:
                        # updating leaves set
                        new_leaves.append(leaf)

            if (leaf == root.id):
                break

            temp_new_leaves = []
            for new_leaf in new_leaves:
                if (new_leaf not in temp_new_leaves):
                    temp_new_leaves.append(new_leaf)
            old_leaves = temp_new_leaves.copy()
            new_leaves = []
            labelled = {}
            for node in nodes:
                labelled[node.id] = False

        return coalescent_process, clade_set_into_root

    def __coalescentRecurse(self, id, distance, clade_set, coalescent_process):
        """
        This is the recursive part of the multi-species coalescent process:
            Given a set of n genes gathering into a branch in the species tree from the bottom,
            whenever we come across a point of coalescence, we randomly merge 2 elements in the gene sets,
            and record the set before the new coalescence, named "from_set", and the set after the coalescence,
            named "to_set", and the distance from the last coalescent event or the bottom of the branch.
        """
        if (len(clade_set[id]) <= 1):
            return clade_set[id]
        else:
            # rate of coalescence
            lambda_c = len(clade_set[id]) * self.__getLambdaCoal(clade_set[id])    
            distance_fake = np.random.exponential(scale=1.0/lambda_c)

            # no coalescent event anymore in this branch
            if (distance < distance_fake):      
                return clade_set[id]
            else:
                # when coalescent, randomly merge 2 elements in the gene sets
                if (len(clade_set[id]) >= 2):   
                    temp_set = sorted(clade_set[id])
                    couple = np.random.choice(clade_set[id], size=2, replace=False)
                    clade_set[id] = [''.join(self.__starSorted(couple))] + [e for e in clade_set[id] if e not in couple]

                    # print process
                    # Debug.log(header="initial node " + str(id) + ": " + str(temp_set) + '\n')
                    # Debug.log(header="coalescent at node " + str(id) + ": " + str(clade_set[id]) + ", " + "distance = " + str(distance_fake) + '\n')

                    # save process
                    coalescent_process[str(id)].append({
                        'from_set': temp_set, 
                        'to_set': clade_set[id].copy(),
                        'distance': distance_fake
                    })
                else:
                    # stop when gene set only has one single element
                    return clade_set[id]

                distance = distance - distance_fake

                # use recursion to simulate the case when there is more than one coalescent events in the branch
                self.__coalescentRecurse(id=id, 
                                         distance=distance, 
                                         clade_set=clade_set,
                                         coalescent_process=coalescent_process)     
        return clade_set[id]

    def __getLambdaCoal(self, clade_set):
        indices = []
        for clade in clade_set:
            splited = clade.split('*')[:-1]
            for index in splited:
                indices.append(int(index))
        return mean(self.__parameter['coalescent'][indices])
    
    def __starInSet(self, target, clade):
        """
        checking whether a given clade is in the target set
        modified for the "*" representation
        """
        if (len(target) <= len(clade)):
            splited_target = target.split('*')[:-1]
            splited_clade = clade.split('*')[:-1]
            return set(splited_target).issubset(set(splited_clade))
        else:
            return False

    def __starSorted(self, couple):
        string = ''
        for e in couple:
            string += e
        splited = string.split('*')[:-1]
        splited = sorted([int(e) for e in splited])
        return [str(e) + '*' for e in splited]