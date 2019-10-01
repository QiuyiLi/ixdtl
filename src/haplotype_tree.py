import skbio
import numpy as np
from collections import defaultdict
from statistics import mean
from .tree_table import *
from .ixdtl_model import *


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

    @property
    def recombination(self):
        return self.__recombination
    
    @property
    def hemiplasy(self):
        return self.__hemiplasy

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

    def initialize(self, locusTree, coalescentProcess=None):
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
        self.readFromSkbioTree(skbioTree)

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

    def dtSubtree(self, coalescentProcess, events, path=None):
        """
        1. find all the duplication points on the coalescent tree 
        2. find the corresponding duplicaion subtree
        3. do subtree coalescence to obtain the sub_coalescent_tree
        4. find all the duplication points on the sub_coalescent_tree
        5. recurse
        """
        # if path:
        #     f = open(os.path.join(path, 'gene_tree.txt'), 'w')
        #     f.write(str(self.skbio_tree))
        #     f.close()
        #     f = open(os.path.join(path, 'species_tree.txt'), 'w')
        #     f.write(str(self.species_tree.skbio_tree))
        #     f.close()
        #     self.find_ils(path)

        # GeneTree.full_events += events

        for event in events:
            # a = event['geneNodeName'].split('*')[:-1]
            # for i in range(len(a)):
            #     a[i] = SpeciesTree.global_species_tree.get_fake_id_from_real_id(a[i])

            if (event['type'] == 'duplication'):
                speciesId, distanceAboveSpeciesNode = self.__mapEventToSpeciesTree(geneId=event['geneNodeId'],
                    eventHeight=event['eventHeight'], speciesId=None)
                self.__dtSubtreeRecurse(
                    event=event, nodeId=speciesId, 
                    coalescentDistance=distanceAboveSpeciesNode, path=path)

                # if coalescentProcess:
                #     # non-trivial case
                #     # find species node id for the event and the distance from the event point to the species node
                #     for speciesNodeId, mergingSets in coalescentProcess.items():
                #         for mergingSet in mergingSets:
                #             if (event['geneNodeName'] in mergingSet['toSet'] 
                #                 and event['geneNodeName'] not in mergingSet['fromSet']):
                #                 speciesId = speciesNodeId
                #                 coalescentDistance = mergingSet['distance']
                #     if speciesId == None:
                #         # does not coalese with others (geneNodeName = 1*)
                #         speciesId = int(event['geneNodeName'][:-1])
                #         coalescentDistance = 0
                #     self.__dtSubtreeRecurse(
                #         event=event, nodeId=nodeId, 
                #         coalescentDistance=coalescentDistance, path=path)
                # else:       
                #     # trivial
                #     nodeId = int(event['geneNodeName'][:-1])
                #     coalescentDistance = 0
                #     self.__dtSubtreeRecurse(
                #         event=event, nodeId=nodeId, 
                #         coalescentDistance=coalescentDistance, path=path)


            elif (event['type'] == 'transfer'):
                trans_target_id = event['target_species_id']
                target_height = SpeciesTree.global_species_tree.distance_to_leaf(trans_target_id, 0)
                distance_above_target = event['event_height'] - target_height
                self.__dtSubtreeRecurse(event=event, nodeId=trans_target_id, coalescentDistance=distance_above_target, path=path)
                # insert sub-tree
            elif (event['type'] == 'loss'):
                # index = Utility.increment()
                # event['index'] = index
                # Debug.event_count['L'] += 1
                # file_name = 'loss_' + str(event['distance_to_gene_node'])
                # f = open(os.path.join(path, file_name), 'w')
                # f.write(str(event['geneNodeName']) + ',' + str(event['distance_to_gene_node']) + ',' + str(index))
                # f.close()
                IxDTLModel.geneSkbioTree.remove_deleted(
                    lambda x: x.name == event['geneNodeName'])
                IxDTLModel.geneSkbioTree.prune()


    def __dtSubtreeRecurse(self, event, newLocusRootId, distanceAboveRoot, path=None):

        if (event['type'] == 'transfer'): 
            # nodeId = target_id
            speciesSkbioTree = self.speciesTree.getSkbioTree()
            newLocusRootName = self.speciesTree.getNodeById(newLocusRootId).name

            newLocusSkbioTree = speciesSkbioTree.find(newLocusRootName).deepcopy()
            newLocusTreeNames = [node.name for node in newLocusSkbioTree.traverse()]
            newLocusTreeNodes = [node for node in self.speciesTree.getNodes() 
                if node.name in newLocusTreeNames]
            newLocusTree = LocusTree(randomState=self.randomState)
            newLocusTree.initialize(nodes=newLocusTreeNodes, skbioTree=newLocusSkbioTree)

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
                coalescentProcess=locusTreeCoalescentProcess)

            rootLength = event['eventHeight'] - newHaplotypeTree.getTreeHeight()
            newHaplotypeTree.getSkbioTree().length = rootLength
            newHaplotypeTreeEvents = newHaplotypeTree.dtlProcess(
                event=event, distance=rootLength)
            
            # index = Utility.increment()
            # event['index'] = index
            # _id = 'trans_subtree_' + str(index)
            # next_dir = os.path.join(path, _id)
            # os.mkdir(next_dir)
            # file_name = 'event.txt'
            # f = open(os.path.join(next_dir, file_name), 'w')
            # f.write(str(event['geneNodeName']) + ',' + str(event['distance_to_gene_node']) + ',' + str(event['type']) + ',' + str(index))
            # f.close()
            newHaplotypeTree.dtSubtree(
                coalescentProcess=locusTreeCoalescentProcess, 
                events=newHaplotypeTreeEvents, path=path)

        if (event['type'] == 'duplication'):
            Debug.log(header='\n\n\n' + '='*80 + '\nCurrent event:' + '\n',
                      bodies=[event], pformat=True)
            species_skbio_tree = self.species_tree.skbio_tree
            name = self.species_tree.nodes_id_dict[nodeId].name

            subtree = species_skbio_tree.find(name).deepcopy()
            subtree_names = [node.name for node in subtree.traverse()]
            subtree_nodes = [node for node in self.species_tree.nodes if node.name in subtree_names]

            species_subtree = SpeciesTree(nodes=subtree_nodes)
            species_subtree.skbio_tree = subtree
            Debug.log(header='\nspecies_subtree_nodes:\n', bodies=species_subtree.nodes)

            distance_above_root = coalescentDistance
            sub_leaves = [int(nodeId) for nodeId in event['geneNodeName'].strip().split('*')[:-1]]
            Debug.log(header='\nspecies_subtree_coal:\n')
            if (GeneTree.recombination == 1):
                if (GeneTree.hemiplasy == 1):
                    species_subtree_coal_process, chosen_gene_name = species_subtree.incomplete_coalescent(distance_above_root, recombination=1)
                elif (GeneTree.hemiplasy == 0):
                    species_subtree_coal_process = species_subtree.bounded_coalescent(distance_above_root, recombination=1)       
                # original_gene_name = event['geneNodeName']
                # chosen_gene_set = set(chosen_gene_name.split('*')[:-1])
                # original_gene_set = set(original_gene_name.split('*')[:-1])
                # union_gene_set = original_gene_set.union(chosen_gene_set)
                # print(union_gene_set)
            elif (GeneTree.dup_recombination == 0):
                if (GeneTree.hemiplasy == 1):
                    species_subtree_coal_process, chosen_gene_name = species_subtree.incomplete_coalescent(distance_above_root, recombination=0, sub_leaves=sub_leaves)
                elif (GeneTree.hemiplasy == 0):
                    species_subtree_coal_process = species_subtree.bounded_coalescent(distance_above_root, recombination=0, sub_leaves=sub_leaves)        

            Debug.log(header='\nspecies_subtree_coal_process:\n',
                      bodies=[species_subtree_coal_process], pformat=True)

            species_subtree_time_seq = species_subtree.time_sequences(coalescentProcess=species_subtree_coal_process)
            Debug.log(header='\nspecies_subtree_time_seq:\n',
                      bodies=[species_subtree_time_seq], pformat=True)

            # save subtree
            # Debug.save_tree_nodes(nodes=species_subtree.nodes, 
            #                       path=Debug.subtree_file_name('output/subtrees', 'dup', nodeId, distance_above_root), 
            #                       distance=distance_above_root)
            
            gene_subtree = GeneTree(time_sequences=species_subtree_time_seq, 
                                    species_tree=species_subtree, 
                                    coalescentProcess=species_subtree_coal_process)
            gene_subtree.skbio_tree.length = event['event_height'] - gene_subtree.total_distance
            Debug.log(header='\ngene_subtree nodes:\n', bodies=gene_subtree.nodes)
            # Debug.save_tree_nodes(nodes=gene_subtree.nodes, 
            #                       path=Debug.subtree_file_name('output/subtrees', 'dup', nodeId, distance_above_root), 
            #                       mode='a')

            # Debug.save_output(contents=[gene_subtree.skbio_tree],
            #                   path=Debug.subtree_file_name('output/newick_gene_subtrees', 'dup', nodeId, distance_above_root))

            Debug.log(header='\ngene_subtree dlt_process:\n')
            gene_subtree_height = gene_subtree.total_distance
            gene_subtree_events = gene_subtree.dlt_process(event=event, distance=event['event_height'] - gene_subtree_height)
            Debug.log(header='\ngene_subtree events:\n',
                      bodies=[gene_subtree_events], pformat=True)

            index = Utility.increment()
            event['index'] = index
            _id = 'dup_subtree_' + str(index)
            next_dir = os.path.join(path, _id)
            os.mkdir(next_dir)
            file_name = 'event.txt'
            f = open(os.path.join(next_dir, file_name), 'w')
            f.write(str(event['geneNodeName']) + ',' + str(event['distance_to_gene_node']) + ',' + str(event['type']) + ',' + str(index))
            f.close()
            gene_subtree.dt_subtree(coalescentProcess=species_subtree_coal_process, events=gene_subtree_events, path=next_dir)