import skbio
import numpy as np
from collections import defaultdict
from statistics import mean
from .tree_table import *
from .locus_tree import *


class HaplotypeTree:
    """
    Nodes are represented in Tree Table which is introduced in tree_table.py
    """

    def __init__(self, randomState, speciesTree):
        self.__randomState = randomState
        self.__speciesTree = speciesTree

        self.__coalescentProcess = None
        self.__treeTable = None
        self.__eventRates = {}
        self.__recombination = None
        self.__hemiplasy = None

    def __repr__(self):
        return str(self.__treeTable)

    def __str__(self):
        return str(self.__treeTable)

    @property
    def randomState(self):
        return self.__randomState

    @property
    def speciesTree(self):
        return self.__speciesTree

    @property
    def coalescentProcess(self):
        return self.__coalescentProcess

    @property
    def treeTable(self):
        return self.__treeTable

    @property
    def eventRates(self):
        return self.__eventRates
    @eventRates.setter
    def eventRates(self, eventRates):
        self.__eventRates = eventRates

    @property
    def recombination(self):
        return self.__recombination
    
    @property
    def hemiplasy(self):
        return self.__hemiplasy

    def getSkbioTree(self):
        return self.__treeTable.skbioTree
    def setSkbioTree(self, skbioTree):
        self.__treeTable.skbioTree = skbioTree

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

    def initialize(self, locusTree, coalescentProcess=None, rename=True):
        if not coalescentProcess:
            coalescentProcess, cladeSetIntoRoot = locusTree.coalescent(
                distanceAboveRoot=float('inf'))
        self.__coalescentProcess = coalescentProcess
        print('coalescent process:')
        print(coalescentProcess)
        print()

        timeSequences = locusTree.getTimeSequences(
            coalescentProcess=coalescentProcess)
        print('time sequences:')
        print(timeSequences)
        print()

        skbioTree = self.createSkbioTree(timeSequences)
        self.readFromSkbioTree(skbioTree, rename)

    def setEventRates(self, duplicationPrmt, transferPrmt, lossPrmt):
        self.__eventRates['d'] = self.randomState.gamma(
            shape=duplicationPrmt['shape'], scale=duplicationPrmt['scale'],
            size=len(self.getLeaves()))
        self.__eventRates['t'] = self.randomState.gamma(
            shape=transferPrmt['shape'], scale=transferPrmt['scale'],
            size=len(self.getLeaves()))
        self.__eventRates['l'] = self.randomState.gamma(
            shape=lossPrmt['shape'], scale=lossPrmt['scale'],
            size=len(self.getLeaves()))

    def setRecombination(self, recombination):
        self.__recombination = recombination

    def setHemiplasy(self, hemiplasy):
        self.__hemiplasy = hemiplasy

    def createSkbioTree(self, timeSequences):
        skbioTree = skbio.tree.TreeNode()
        tempTimeSequences = timeSequences.copy()
        for k, v in timeSequences.items():
            if not v: del tempTimeSequences[k]
        if len(tempTimeSequences) > 0:
            skbioTree.name = next(iter(tempTimeSequences.values()))[-1][0]
            self.__createSkbioTreeRecurse(
                skbioTree=skbioTree, timeSequences=timeSequences)
            skbioTree.length = None
        else:
            skbioTree.name = str(next(iter(timeSequences))) + '*' # speciesid*
            skbioTree.length = None
        return skbioTree

    def __createSkbioTreeRecurse(self, skbioTree, timeSequences):
        # one node (leaf)
        if skbioTree.name.count('*') == 1:
            skbioTree.length = self.__distanceToParent(
                nodeName=skbioTree.name, parentName=skbioTree.parent.name, 
                timeSequences=timeSequences)
            return
        # two nodes
        elif len(skbioTree.name) == 4:
            childLName = skbioTree.name[:2]
            childRName = skbioTree.name[2:]
            childL = skbio.tree.TreeNode(
                name=childLName, 
                length=self.__distanceToParent(
                    nodeName=childLName, parentName=skbioTree.name, 
                    timeSequences=timeSequences), 
                parent=skbioTree)
            childR = skbio.tree.TreeNode(
                name=childRName, 
                length=self.__distanceToParent(
                    nodeName=childRName, parentName=skbioTree.name, 
                    timeSequences=timeSequences),
                parent=skbioTree)
            skbioTree.children = [childL, childR]
            return
        # otherwise
        else:
            isFound = False
            for _, sequence in timeSequences.items():
                prevPair = None
                for pair in sequence:
                    if (prevPair != None and skbioTree.name == pair[0]):
                        childLName = prevPair[0]
                        childRName = self.__starReplace(
                            skbioTree.name, prevPair[0])
                        childL = skbio.tree.TreeNode(
                            name=childLName, 
                            length=self.__distanceToParent(
                                nodeName=childLName, 
                                parentName=skbioTree.name, 
                                timeSequences=timeSequences), 
                            parent=skbioTree)
                        childR = skbio.tree.TreeNode(
                            name=childRName, 
                            length=self.__distanceToParent(
                                nodeName=childRName, 
                                parentName=skbioTree.name, 
                                timeSequences=timeSequences),
                            parent=skbioTree)
                        self.__createSkbioTreeRecurse(
                            skbioTree=childL, timeSequences=timeSequences)
                        self.__createSkbioTreeRecurse(
                            skbioTree=childR, timeSequences=timeSequences)
                        skbioTree.children = [childL, childR]
                        isFound = True
                        break
                    prevPair = pair
                if isFound: break

    def readFromSkbioTree(self, skbioTree, rename=True):
        self.__treeTable = TreeTable()
        self.__treeTable.createFromSkbioTree(skbioTree, rename)

    def dtlProcess(self, distance, event=None):
        events = []
        if len(self.getNodes()) == 1:
            # trivial case
            distance = event['distanceToGeneNode']
        self.__dtlProcessRecurse(
            skbioTreeNode=self.getSkbioTree(), distance=distance, events=events)
        return events

    def __dtlProcessRecurse(self, skbioTreeNode, distance, events):
        node = self.getNodeByName(skbioTreeNode.name)

        distanceD = self.randomState.exponential(
            scale=1.0 / self.__getEventRateInAncestralBranch(
                eventType='d', clade=node.name))
        distanceT = self.randomState.exponential(
            scale=1.0 / self.__getEventRateInAncestralBranch(
                eventType='t', clade=node.name))
        distanceL = self.randomState.exponential(
            scale=1.0 / self.__getEventRateInAncestralBranch(
                eventType='l', clade=node.name))

        # duplication happens first
        if (distanceD < min(distanceL, distanceT) and distanceD < distance):
            eventHeight = self.getDistanceToLeaf(node.id, 0) + distance - distanceD
            speciesId, distanceAboveSpeciesNode = self.__mapEventToSpeciesTree(
                geneId=node.id, eventHeight=eventHeight, speciesId=None)
            events.append({
                'type': 'duplication',
                'geneNodeId': node.id, 
                'geneNodeName': node.name, 
                'distanceToGeneNode': distance - distanceD,
                'eventHeight': eventHeight,
                'speciesNodeId': speciesId,
                'distanceToSpeciesNode': distanceAboveSpeciesNode,
                'index': -1
            })
            # looking for more events on the same branch
            self.__dtlProcessRecurse(
                skbioTreeNode=skbioTreeNode, 
                distance=distance - distanceD, events=events)
        elif (distanceT <= min(distanceD, distanceL) and distanceT < distance):
            eventHeight = self.getDistanceToLeaf(node.id, 0) + distance - distanceT
            speciesTreeHeight = self.speciesTree.getTreeHeight()
            if eventHeight < speciesTreeHeight:
                target, originalSpeciesId = self.__findTransferTarget(
                    eventHeight=eventHeight, geneId=node.id)
                speciesId, distanceAboveSpeciesNode = self.__mapEventToSpeciesTree(
                    geneId=node.id, eventHeight=eventHeight, speciesId=None)
                if target:
                    events.append({
                        'type': 'transfer',
                        'geneNodeId': node.id, 
                        'geneNodeName': node.name, 
                        'distanceToGeneNode': distance - distanceT,
                        'targetSpeciesId': target,
                        'eventHeight': eventHeight,
                        'speciesNodeId': speciesId,
                        'distanceToSpeciesNode': distanceAboveSpeciesNode,
                        'index': -1
                    })
            self.__dtlProcessRecurse(
                skbioTreeNode=skbioTreeNode, 
                distance=distance - distanceT, events=events)
        elif (distanceL <= min(distanceD, distanceT) and distanceL < distance):      
            # loss happens first, the seaching process stops at the loss point
            eventHeight = self.getDistanceToLeaf(node.id, 0) + distance - distanceL
            speciesId, distanceAboveSpeciesNode = self.__mapEventToSpeciesTree(
                geneId=node.id, eventHeight=eventHeight, speciesId=None)
            # speciesId = self.node_by_id(speciesId).parent_id
            events.append({
                'type': 'loss',
                'geneNodeId': node.id, 
                'geneNodeName': node.name, 
                'distanceToGeneNode': distance - distanceL,
                'eventHeight': eventHeight,
                'speciesNodeId': speciesId,
                'distanceToSpeciesNode': distanceAboveSpeciesNode,
                'index': -1
            })
        else:
            # reach the end the current branch, looking for events in the 2 children branches
            # print('nothing happened at node ' + str(node.id) + ' (' + node.name + ')' + '\n')
            if (node.children):     # if children branches exist
                childL = skbioTreeNode.children[0]
                childR = skbioTreeNode.children[1]
                distanceToChildL = node.distanceToChildren[0]
                distanceToChildR = node.distanceToChildren[1]
                self.__dtlProcessRecurse(
                    skbioTreeNode=childL, 
                    distance=distanceToChildL, events=events)
                self.__dtlProcessRecurse(
                    skbioTreeNode=childR, 
                    distance=distanceToChildR, events=events)
            # else:       # if not exist, reach the leaves of the tree, searching process stops
            #     print('reach the end of node ' + str(node.id) + ' (' + node.name + ')' + '\n')

    def __mapGeneIdToSpeciesId(self, geneId):
        speciesId = None
        geneName = self.getNodeById(geneId).name
        if self.coalescentProcess:
            # non-trivial case
            for speciesNodeId, mergingSets in self.coalescentProcess.items():
                for mergingSet in mergingSets:
                    if (geneName in mergingSet['toSet'] 
                        and geneName not in mergingSet['fromSet']):
                        speciesId = speciesNodeId
            if speciesId == None:
                speciesId = int(geneName[:-1])
        else:
            # trivial case
            speciesId = int(geneName[:-1])     
        return speciesId

    # M function
    def __mapEventToSpeciesTree(self, geneId, eventHeight, speciesId=None):
        if speciesId == None:
            speciesId = self.__mapGeneIdToSpeciesId(geneId=geneId)
        distanceAboveSpeciesNode = eventHeight - self.speciesTree.getDistanceToLeaf(speciesId, 0)
        if speciesId == self.speciesTree.getRoot().id:
            return speciesId, distanceAboveSpeciesNode
        else:
            speciesIdParent = self.speciesTree.getNodeById(speciesId).parent
            speciesDistanceParent = self.speciesTree.getDistanceToLeaf(speciesIdParent, 0)
            if speciesDistanceParent > eventHeight:
                return speciesId, distanceAboveSpeciesNode
            else:
                return self.__mapEventToSpeciesTree(
                    geneId=geneId, eventHeight=eventHeight, speciesId=speciesIdParent)

    def __findTransferTarget(self, eventHeight, geneId):
        speciesNodes = self.speciesTree.getNodes()
        originSpeciesId = self.__mapGeneIdToSpeciesId(geneId=geneId)
        nodesList = []
        for node in speciesNodes:
            if node.id == originSpeciesId:
                continue
            if node.id == self.speciesTree.getRoot().id:
                continue
            parentHeight = self.speciesTree.getDistanceToLeaf(
                self.speciesTree.getNodeById(node.parent).id, 0)
            if parentHeight > eventHeight:
                nodeHeight = self.speciesTree.getDistanceToLeaf(node.id, 0)
                if nodeHeight <= eventHeight:
                    nodesList.append(node.id)
        return self.randomState.choice(nodesList), originSpeciesId

    def __getEventRateInAncestralBranch(self, eventType, clade):
        indices = []
        splited = clade.split('*')[:-1]
        for index in splited:
            indices.append(int(index))
        return mean(self.eventRates[eventType][indices])

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

    def dtSubtree(self, coalescentProcess, events, haplotypeTree, level):
        """
        1. find all the duplication points on the coalescent tree 
        2. find the corresponding duplicaion subtree
        3. do subtree coalescence to obtain the sub_coalescent_tree
        4. find all the duplication points on the sub_coalescent_tree
        5. recurse
        fixation: at each level, update the parent haplotype tree finishing each event
        """
        eventIndex = -1
        for event in events:
            if (event['type'] == 'duplication'):
                eventIndex = eventIndex + 1
                speciesId, distanceAboveSpeciesNode = self.__mapEventToSpeciesTree(
                    geneId=event['geneNodeId'], eventHeight=event['eventHeight'], 
                    speciesId=None)
                newHaplotypeTree = self.__dtSubtreeRecurse(
                    event=event, newLocusRootId=speciesId, 
                    distanceAboveRoot=distanceAboveSpeciesNode, level=level)

                for node in newHaplotypeTree.getSkbioTree().traverse():
                    node.name = node.name + '_lv=' + str(level) + '_id=' + str(eventIndex)
                
                print('new haplotype tree:')
                print(newHaplotypeTree.getSkbioTree().ascii_art())

                # insert newHaplotypeTree
                print('haplotype tree before:')
                print(haplotypeTree.getSkbioTree().ascii_art())

                geneNodeName = event['geneNodeName']
                geneNode = haplotypeTree.getSkbioTree().find(geneNodeName)
                geneNodeParent = geneNode.parent
                # 1. create new node
                newNode = skbio.tree.TreeNode()
                newNode.name = 't' + '_lv=' + str(level) + '_id=' + str(eventIndex)
                # 4. change length
                newHaplotypeTree.getSkbioTree().length = event['eventHeight'] - newHaplotypeTree.getTreeHeight()
                print('newHaplotypeTree length: ' + str(newHaplotypeTree.getSkbioTree().length))
                print('nht tree height: ' + str(newHaplotypeTree.getTreeHeight()))
                newNode.length = 0 if not geneNodeParent else geneNode.length - event['distanceToGeneNode']
                print('newNode.length: ' + str(newNode.length))
                geneNode.length = event['distanceToGeneNode']
                print('geneNode.length: ' + str(geneNode.length))
                # 2. change children
                newNode.children = [geneNode, newHaplotypeTree.getSkbioTree()]
                if geneNodeParent:
                    for i in range(len(geneNode.parent.children)):
                        if geneNode.parent.children[i] == geneNode:
                            geneNode.parent.children[i] = newNode
                            break
                # 3. change parent
                geneNode.parent = newNode
                newHaplotypeTree.getSkbioTree().parent = newNode

                # check root
                if not geneNodeParent:
                    haplotypeTree.setSkbioTree(newNode)
                else:
                    newNode.parent = geneNodeParent

                print('haplotype tree after:')
                print(haplotypeTree.getSkbioTree().ascii_art())

            elif (event['type'] == 'transfer'):
                eventIndex = eventIndex + 1
                transferTargetId = event['targetSpeciesId']
                targetHeight = self.speciesTree.getDistanceToLeaf(transferTargetId, 0)
                distanceAboveTarget = event['eventHeight'] - targetHeight
                newHaplotypeTree = self.__dtSubtreeRecurse(
                    event=event, newLocusRootId=transferTargetId, 
                    distanceAboveRoot=distanceAboveTarget, level=level)
                
                for node in newHaplotypeTree.getSkbioTree().traverse():
                    node.name = node.name + '_lv=' + str(level) + '_id=' + str(eventIndex)
                
                print('new haplotype tree:')
                print(newHaplotypeTree.getSkbioTree().ascii_art())

                # insert newHaplotypeTree
                print('haplotype tree before:')
                print(haplotypeTree.getSkbioTree().ascii_art())
                
                geneNodeName = event['geneNodeName']
                geneNode = haplotypeTree.getSkbioTree().find(geneNodeName)
                geneNodeParent = geneNode.parent
                # 1. create new node
                newNode = skbio.tree.TreeNode()
                newNode.name = 't' + '_lv=' + str(level) + '_id=' + str(eventIndex)
                # 4. change length
                newHaplotypeTree.getSkbioTree().length = event['eventHeight'] - newHaplotypeTree.getTreeHeight()
                print('newHaplotypeTree length: ' + str(newHaplotypeTree.getSkbioTree().length))
                print('nht tree height: ' + str(newHaplotypeTree.getTreeHeight()))
                newNode.length = 0 if not geneNodeParent else geneNode.length - event['distanceToGeneNode']
                print('newNode.length: ' + str(newNode.length))
                geneNode.length = event['distanceToGeneNode']
                print('geneNode.length: ' + str(geneNode.length))
                # 2. change children
                newNode.children = [geneNode, newHaplotypeTree.getSkbioTree()]
                if geneNodeParent:
                    for i in range(len(geneNode.parent.children)):
                        if geneNode.parent.children[i] == geneNode:
                            geneNode.parent.children[i] = newNode
                            break
                # 3. change parent
                geneNode.parent = newNode
                newHaplotypeTree.getSkbioTree().parent = newNode
                
                # check root
                if not geneNodeParent:
                    haplotypeTree.setSkbioTree(newNode)
                else:
                    newNode.parent = geneNodeParent

                print('haplotype tree after:')
                print(haplotypeTree.getSkbioTree().ascii_art())

            elif (event['type'] == 'loss'):
                print('loss')

                # cut tree bug here...
                # do not cut it for now, lable the loss points and delete all of them all at once when everything done
                # haplotypeSkbioTree = haplotypeTree.getSkbioTree()
                # haplotypeSkbioTree.remove_deleted(
                #     lambda x: x.name == event['geneNodeName'])
                # haplotypeSkbioTree.prune()

        return haplotypeTree

    def __dtSubtreeRecurse(self, event, newLocusRootId, distanceAboveRoot, level):
        if (event['type'] == 'duplication' or event['type'] == 'transfer'): 
            # nodeId = target_id
            speciesSkbioTree = self.speciesTree.getSkbioTree()
            newLocusRootName = self.speciesTree.getNodeById(newLocusRootId).name

            newLocusSkbioTree = speciesSkbioTree.find(newLocusRootName).deepcopy()
            newLocusTreeNames = [node.name for node in newLocusSkbioTree.traverse()]
            newLocusTreeNodes = [node for node in self.speciesTree.getNodes() 
                if node.name in newLocusTreeNames]
            newLocusTree = LocusTree(randomState=self.randomState)
            newLocusTree.initialize(nodes=newLocusTreeNodes, skbioTree=newLocusSkbioTree)
            newLocusTree.coalescentRate = self.speciesTree.coalescentRate

            locusTreeCoalescentProcess = None
            if self.hemiplasy == 1:
                locusTreeCoalescentProcess, chosenGeneName = \
                    newLocusTree.incompleteCoalescent(distanceAboveRoot)
            elif self.hemiplasy == 0:
                locusTreeCoalescentProcess = \
                    newLocusTree.boundedCoalescent(distanceAboveRoot)

            newHaplotypeTree = HaplotypeTree(
                randomState=self.randomState, speciesTree=self.speciesTree)
            newHaplotypeTree.initialize(
                locusTree=newLocusTree, 
                coalescentProcess=locusTreeCoalescentProcess, rename=False)
            newHaplotypeTree.eventRates = self.eventRates

            rootLength = event['eventHeight'] - newHaplotypeTree.getTreeHeight()
            newHaplotypeTree.getSkbioTree().length = rootLength
            
            newHaplotypeTreeEvents = newHaplotypeTree.dtlProcess(
                event=event, distance=rootLength)
            newHaplotypeTreeEvents.sort(reverse=True, key=lambda x: x['eventHeight'])

            print('new haplotype tree events level ' + str(level + 1) + ':')
            print(newHaplotypeTreeEvents)
            print()

            newHaplotypeTree.dtSubtree(
                coalescentProcess=locusTreeCoalescentProcess, 
                events=newHaplotypeTreeEvents, haplotypeTree=newHaplotypeTree, 
                level=level + 1)
            return newHaplotypeTree