import numpy as np
from .species_tree import *
from .exception import *


class HIDTLModel:
    def __init__(self, seed=0):
        self.__speciesTree = None
        self.__haplotypeTree = None
        self.__locusTrees = []
        self.__parameters = {}
        self.__randomState = np.random.RandomState(seed)

    @property
    def speciesTree(self):
        return self.__speciesTree

    @property
    def haplotypeTree(self):
        return self.__haplotypeTree

    @property
    def locusTrees(self):
        return self.__locusTrees

    @property
    def parameters(self):
        return self.__parameters

    def run(self, inputFile, coalescentArgs, duplicationArgs, transferArgs, 
            lossArgs, hemiplasy, recombination):
        # set parameters
        self.setParameters(
            coalescent=coalescentArgs, 
            duplication=duplicationArgs,
            transfer=transferArgs, 
            loss=lossArgs, 
            hemiplasy=hemiplasy,
            recombination=recombination)

        # read and create a species tree from input file
        self.createSpeciesTree(inputFile)

        # create a haplotype tree according to the species tree
        self.createHaplotypeTree()

    def setParameters(self, coalescent, duplication, transfer, loss, 
                      hemiplasy, recombination):
        if not coalescent:
            raise HIDTLError('missing coalescent parameter')
        self.__parameters['coalescent'] = coalescent

        if not duplication:
            raise HIDTLError('missing duplication parameter')
        self.__parameters['duplication'] = duplication

        if not transfer:
            raise HIDTLError('missing transfer parameter')
        self.__parameters['transfer'] = transfer

        if not loss:
            raise HIDTLError('missing loss parameter')
        self.__parameters['loss'] = loss

        if hemiplasy is None:
            raise HIDTLError('missing hemiplasy option')
        self.__parameters['hemiplasy'] = hemiplasy

        if recombination is None:
            raise HIDTLError('missing recombination option')
        self.__parameters['recombination'] = recombination

    def createSpeciesTree(self, path):
        self.__speciesTree = SpeciesTree(self.__randomState)
        self.__speciesTree.readFromNewickFile(path)
        self.__speciesTree.setLambdaCoalescent(self.__parameters['coalescent'])
            
    def createHaplotypeTree(self):
        coalescentProcess, cladeSetIntoRoot = self.__speciesTree.coalescent(
            distanceAboveRoot=10000)
        print(coalescentProcess)
