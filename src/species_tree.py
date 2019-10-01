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
        self.__coalescentRate = None

    def __repr__(self):
        return str(self.__treeTable)

    def __str__(self):
        return str(self.__treeTable)

    @property
    def randomState(self):
        return self.__randomState

    @property
    def treeTable(self):
        return self.__treeTable

    @property
    def coalescentRate(self):
        return self.__coalescentRate

    def setCoalescentRate(self, parameter):
        self.__coalescentRate = self.randomState.gamma(
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

    def getTreeHeight(self):
        return self.__treeTable.treeHeight

    def getDistanceToLeaf(self, nodeId, branchDistance):
        return self.__treeTable.distanceToLeaf(nodeId, branchDistance)

    def readFromNewickFile(self, path):
        self.__treeTable = TreeTable()
        self.__treeTable.createFromNewickFile(path)

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
            coalescentRate = len(cladeSet[id]) \
                      * self.__getCoalescentRateInAncestralBranch(cladeSet[id])
            fakeDistance = self.randomState.exponential(
                scale=1.0 / coalescentRate)

            # no coalescent event anymore in this branch
            if distance < fakeDistance:
                return cladeSet[id]
            else:
                # when coalescent, randomly merge 2 elements in the gene sets
                if len(cladeSet[id]) >= 2:
                    temp_set = sorted(cladeSet[id])
                    couple = self.randomState.choice(
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

    def __getCoalescentRateInAncestralBranch(self, cladeSet):
        indices = []
        for clade in cladeSet:
            splited = clade.split('*')[:-1]
            for index in splited:
                indices.append(int(index))
        return mean(self.coalescentRate[indices])

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

    def getTimeSequences(self, coalescentProcess):
        """
        backward-in-time coalescent process modified data structure 
        for constructing the coalescent tree in newick format
        """
        timeSequences = {}
        for leaf in self.getLeaves():
            sequence = self.__findAncestors(
                leafName=str(leaf.id) + '*', 
                coalescentProcess=coalescentProcess)
            if sequence:
                timeSequences[leaf.id] = sequence
        return timeSequences

    def __findAncestors(self, leafName, coalescentProcess):
        """
        find the ancestors of the given leaf in reverse time order
        """
        sequence = []
        for speciesNodeId, mergingSets in coalescentProcess.items():
            branchDistance = 0.0
            for mergingSet in mergingSets:
                branchDistance += mergingSet['distance']
                if (leafName in mergingSet['fromSet'] 
                    and leafName not in mergingSet['toSet']):
                    for element in mergingSet['toSet']:
                        if (len(leafName) < len(element) 
                            and self.__starInSet(leafName, element)):
                            coalescentHeight = self.getDistanceToLeaf(
                                nodeId=speciesNodeId, 
                                branchDistance=branchDistance)
                            pair = (element, coalescentHeight)
                            sequence.append(pair)
                            sequence += self.__findAncestors(
                                leafName=element, 
                                coalescentProcess=coalescentProcess)
        return sequence
