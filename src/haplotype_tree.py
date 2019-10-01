import skbio
import numpy as np
from collections import defaultdict
from statistics import mean
from .tree_table import *


class HaplotypeTree:
    """
    Nodes are represented in Tree Table which is introduced in tree_table.py
    """

    def __init__(self, randomState):
        self.__randomState = randomState

        self.__treeTable = None
        self.__eventRates = {}
        self.__recombination = None
        self.__hemiplasy = None

    def __repr__(self):
        return str(self.__treeTable)

    def __str__(self):
        return str(self.__treeTable)

    def initialize(self, locusTree):
        coalescentProcess, cladeSetIntoRoot = locusTree.coalescent(
            distanceAboveRoot=float('inf'))
        print('coalescent process:')
        print(coalescentProcess)
        print()

        timeSequences = locusTree.getTimeSequences(
            coalescentProcess=coalescentProcess)
        print('time sequences:')
        print(timeSequences)
        print()

        skbioTree = self.createSkbioTree(timeSequences)
        self.readFromSkbioTree(skbioTree)

    def setEventRates(self, duplicationRate, transferRate, lossRate):
        self.__eventRates['d'] = duplicationRate
        self.__eventRates['t'] = transferRate
        self.__eventRates['l'] = lossRate

    def setRecombination(self, recombination):
        self.__recombination = recombination

    def setHemiplasy(self, hemiplasy):
        self.__hemiplasy = hemiplasy

    def createSkbioTree(self, timeSequences):
        skbioTree = skbio.tree.TreeNode()   # root node
        if len(timeSequences) > 0:
            skbioTree.name = next(iter(timeSequences.values()))[-1][0]
            self.__createSkbioTreeRecurse(
                skbioTree=skbioTree, timeSequences=timeSequences)
            skbioTree.length = None
        else:
            skbioTree.name = str(next(iter(timeSequences))) + '*'
            skbioTree.length = None
        return skbioTree

    def readFromSkbioTree(self, skbioTree):
        self.__treeTable = TreeTable()
        self.__treeTable.createFromSkbioTree(skbioTree)

    def __createSkbioTreeRecurse(self, skbioTree, timeSequences):
        # one node (leaf)
        if (skbioTree.name.count('*') == 1):
            skbioTree.length = self.__distanceToParent(
                nodeName=skbioTree.name, parentName=skbioTree.parent.name, 
                timeSequences=timeSequences)
            return
        # two nodes
        elif (len(skbioTree.name) == 4):
            childOneName = skbioTree.name[:2]
            childTwoName = skbioTree.name[2:]
            childOne = skbio.tree.TreeNode(
                name=childOneName, 
                length=self.__distanceToParent(
                    nodeName=childOneName, parentName=skbioTree.name, 
                    timeSequences=timeSequences), 
                parent=skbioTree)
            childTwo = skbio.tree.TreeNode(
                name=childTwoName, 
                length=self.__distanceToParent(
                    nodeName=childTwoName, parentName=skbioTree.name, 
                    timeSequences=timeSequences),
                parent=skbioTree)
            skbioTree.children = [childOne, childTwo]
            return
        # otherwise
        else:
            isFound = False
            for _, sequence in timeSequences.items():
                prevPair = None
                for pair in sequence:
                    if (prevPair != None and skbioTree.name == pair[0]):
                        childOneName = prevPair[0]
                        childTwoName = self.__starReplace(
                            skbioTree.name, prevPair[0])
                        childOne = skbio.tree.TreeNode(
                            name=childOneName, 
                            length=self.__distanceToParent(
                                nodeName=childOneName, 
                                parentName=skbioTree.name, 
                                timeSequences=timeSequences), 
                            parent=skbioTree)
                        childTwo = skbio.tree.TreeNode(
                            name=childTwoName, 
                            length=self.__distanceToParent(
                                nodeName=childTwoName, 
                                parentName=skbioTree.name, 
                                timeSequences=timeSequences),
                            parent=skbioTree)
                        self.__createSkbioTreeRecurse(
                            skbioTree=childOne, timeSequences=timeSequences)
                        self.__createSkbioTreeRecurse(
                            skbioTree=childTwo, timeSequences=timeSequences)
                        skbioTree.children = [childOne, childTwo]
                        isFound = True
                        break
                    prevPair = pair
                if isFound: break

    def __distanceToParent(self, nodeName, parentName, timeSequences):
        for leaf, sequence in timeSequences.items():
            if (nodeName.count('*') == 1 and nodeName[0] == str(leaf)):
                for pair in sequence:
                    if pair[0] == parentName:
                        return pair[1]
            else:
                prevPair = None
                for pair in sequence:
                    if (prevPair != None 
                        and prevPair[0] == nodeName 
                        and pair[0] == parentName):
                        return pair[1] - prevPair[1]
                    prevPair = pair
        return None

    def __starReplace(self, string, substring):
        a = string.split('*')[:-1]
        b = substring.split('*')[:-1]
        diff = set(a).difference(set(b))
        return ''.join([e + '*' for e in sorted(list(diff))])