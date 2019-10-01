import skbio
import numpy as np
from collections import defaultdict
from statistics import mean
from .tree_table import *


class SpeciesTree:
    """
    Nodes are represented in Tree Table which is introduced in tree_table.py
    """

    def __init__(self, randomState):
        self.__randomState = randomState

        self.__treeTable = None
        self.__lambdaCoalescent = None

    @property
    def treeTable(self):
        return self.__treeTable

    @property
    def lambdaCoalescent(self):
        return self.__lambdaCoalescent

    def setLambdaCoalescent(self, parameter):
        self.__lambdaCoalescent = self.__randomState.gamma(
            shape=parameter['shape'], scale=parameter['scale'],
            size=len(self.getLeaves()))

    def getSkbioTree(self):
        return self.__treeTable.skbioTree

    def getNodeById(self, id):
        return self.__treeTable.getEntryById(id)

    def getNodes(self):
        return self.__treeTable.table

    def getNodeByName(self, name):
        return self.__treeTable.getEntryByName(name)

    def getRoot(self):
        return self.__treeTable.root

    def getLeaves(self):
        return self.__treeTable.leaves

    def readFromNewickFile(self, path):
        self.__treeTable = TreeTable()
        self.__treeTable.createFromNewickFile(path)

        print(self.__treeTable)

    def coalescent(self, distanceAboveRoot):
        """
        the main multi-species coalecent function
        """
        nodes = self.getNodes()
        root = self.getRoot()
        coalescentProcess = defaultdict(list)
        cladeSetIntoRoot = None

        # leaves of the given species tree
        oldLeaves = [node.id for node in nodes if not node.children]

        # leaves set will be updated in the loop
        newLeaves = []

        # set of extant species that an ancestral gene will finally be fixed in
        cladeSet = {}

        # avoid doing repeated coalescence
        labelled = {}

        # initialize labelled set and clade set
        for node in nodes:
            labelled[node.id] = False
            cladeSet[node.id] = \
                [str(node.id) + '*'] if not node.children else []

        while True:
            for leaf in oldLeaves:
                if leaf == root.id:
                    cladeSetIntoRoot = self.__coalescentRecurse(
                        id=root.id, distance=distanceAboveRoot,
                        cladeSet=cladeSet, coalescentProcess=coalescentProcess)
                    break
                else:
                    parent = self.getNodeById(leaf).parent
                    children = self.getNodeById(parent).children
                    if labelled[leaf]:
                        continue
                    labelled[leaf] = True
                    if (len(cladeSet[children[0]]) != 0 
                        and len(cladeSet[children[1]]) != 0):
                        self.__coalescentRecurse(
                            id=children[0], 
                            distance=self.getNodeById(
                                children[0]).distanceToParent,
                            cladeSet=cladeSet, 
                            coalescentProcess=coalescentProcess)
                        self.__coalescentRecurse(
                            id=children[1], 
                            distance=self.getNodeById(
                                children[1]).distanceToParent,
                            cladeSet=cladeSet, 
                            coalescentProcess=coalescentProcess)
                            
                        # the clade set of the parent before coalescence
                        # is the union of the clade set of its children 
                        # after coalescence
                        cladeSet[parent] = list(set().union(
                            cladeSet[children[0]], cladeSet[children[1]]))
                            
                        if len(newLeaves) > 0:
                            newLeaves = [e for e in newLeaves if e !=
                                         children[0] and e != children[1]]
                        newLeaves.append(parent)
                    else:
                        # updating leaves set
                        newLeaves.append(leaf)

            if leaf == root.id:
                break

            tempNewLeaves = []
            for newLeaf in newLeaves:
                if newLeaf not in tempNewLeaves:
                    tempNewLeaves.append(newLeaf)
            oldLeaves = tempNewLeaves.copy()

            newLeaves = []
            labelled = {}
            for node in nodes:
                labelled[node.id] = False

        return coalescentProcess, cladeSetIntoRoot

    def __coalescentRecurse(self, id, distance, cladeSet, coalescentProcess):
        """
        This is the recursive part of the multi-species coalescent process:
        Given a set of n genes gathering into a branch in the species tree 
        from the bottom, whenever we come across a point of coalescence, 
        we randomly merge 2 elements in the gene sets, and record the set 
        before the new coalescence, named "from_set", and the set after 
        the coalescence, named "to_set", and the distance from the last 
        coalescent event or the bottom of the branch.
        """
        if len(cladeSet[id]) <= 1:
            return cladeSet[id]
        else:
            # rate of coalescence
            lambdaC = len(cladeSet[id]) \
                      * self.__getLambdaCoalescentByCladeSet(cladeSet[id])
            fakeDistance = self.__randomState.exponential(scale=1.0 / lambdaC)

            # no coalescent event anymore in this branch
            if distance < fakeDistance:
                return cladeSet[id]
            else:
                # when coalescent, randomly merge 2 elements in the gene sets
                if len(cladeSet[id]) >= 2:
                    temp_set = sorted(cladeSet[id])
                    couple = self.__randomState.choice(
                        cladeSet[id], 
                        size=2, 
                        replace=False)
                    cladeSet[id] = [''.join(self.__starSorted(couple))] \
                        + [e for e in cladeSet[id] if e not in couple]

                    # save process
                    coalescentProcess[id].append({
                        'fromSet': temp_set,
                        'toSet': cladeSet[id].copy(),
                        'distance': fakeDistance
                    })
                else:
                    # stop when gene set only has one single element
                    return cladeSet[id]

                distance = distance - fakeDistance

                # use recursion to simulate the case when there is
                # more than one coalescent events in the branch
                self.__coalescentRecurse(
                    id=id, 
                    distance=distance,
                    cladeSet=cladeSet, 
                    coalescentProcess=coalescentProcess)

        return cladeSet[id]

    def __getLambdaCoalescentByCladeSet(self, cladeSet):
        indices = []
        for clade in cladeSet:
            splited = clade.split('*')[:-1]
            for index in splited:
                indices.append(int(index))
        return mean(self.__lambdaCoalescent[indices])

    def __starInSet(self, target, clade):
        """
        checking whether a given clade is in the target set
        modified for the "*" representation
        """
        if len(target) <= len(clade):
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

    def __find_ancestors(self, leaf_name, coalescent_process):
        """
        find the ancestors of the given leaf in reverse time order
        """
        sequence = []
        for speciesNodeId, v in coalescent_process.items():
            branch_distance = 0.0
            for elem in v:
                branch_distance += elem['distance']
                if (leaf_name in elem['fromSet'] 
                    and leaf_name not in elem['toSet']):
                    for e in elem['toSet']:
                        if (len(leaf_name) < len(e) 
                            and self.star_in_set(leaf_name, e)):
                            coal_height = super().distance_to_leaf(
                                node_id=speciesNodeId, branch_distance=branch_distance)
                            pair = (e, coal_height)
                            sequence.append(pair)
                            sequence += self.find_ancestors(
                                leaf_name=e, 
                                coalescent_process=coalescent_process)
        return sequence

    def time_sequences(self, coalescent_process):
        """
        backward-in-time coalescent process modified data structure 
        for constructing the coalescent tree in newick format
        """
        time_sequences = {}
        for leaf in self.leaves:
            time_sequences[str(leaf)] = self.find_ancestors(
                leaf_name=str(leaf) + '*', 
                coalescent_process=coalescent_process)
        return time_sequences