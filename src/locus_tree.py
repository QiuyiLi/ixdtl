from .species_tree import *

class LocusTree(SpeciesTree):

    def initialize(self, nodes, skbioTree):
        self.__tableTable = TreeTable()
        self.__tableTable.createFromEntries(
            entries=nodes, skbioTree=skbioTree)

    def boundedCoalescent(self, distanceAboveRoot):
        coalescentProcess, genesIntoRoot = self.coalescent(distanceAboveRoot)
        if len(genesIntoRoot) == 1:
            return coalescentProcess
        else:
            return self.boundedCoalescent(distanceAboveRoot)

    def incompleteCoalescent(self, distanceAboveRoot):
        coalescentProcess, genesIntoRoot = self.coalescent(distanceAboveRoot)
        chosenGene = np.random.choice(genesIntoRoot)
        selectedCoalescentProcess = self.__selectCoalescentProcess(
            coalescentProcess, chosenGene)
        return selectedCoalescentProcess, chosenGene

    def __selectCoalescentProcess(self, coalescentProcess, chosenGene):
        selectedCoalescentProcess = defaultdict(list)
        for speciesNodeId, mergingSets in coalescentProcess.items():
            for mergingSet in mergingSets:
                distance = mergingSet['distance']
                fromSet = []
                toSet = []
                for clade in mergingSet['fromSet']:
                    if self.__starInSet(target=clade, clade=chosenGene):
                        fromSet.append(clade)
                for clade in mergingSet['toSet']:
                    if self.__starInSet(target=clade, clade=chosenGene):
                        toSet.append(clade)
                if toSet:
                    selectedCoalescentProcess[speciesNodeId].append({
                        'fromSet': fromSet, 
                        'toSet': toSet,
                        'distance': distance
                    })
        return selectedCoalescentProcess
