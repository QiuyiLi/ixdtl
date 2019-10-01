import numpy as np
from .species_tree import *
from .haplotype_tree import *
from .exception import *


class IxDTLModel:

    def __init__(self, seed=0):
        self.__randomState = np.random.RandomState(seed)

        self.__speciesTree = None
        self.__haplotypeTree = None
        self.__locusTrees = []
        self.__parameters = {}

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

    @property
    def randomState(self):
        return self.__randomState

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

        # read a species tree from input file
        self.readSpeciesTree(inputFile)

        # construct the original haplotype tree according to the species tree
        self.constructOriginalHaplotypeTree()

        # run dtl process
        events = self.__haplotypeTree.runDTLProcess()

    def setParameters(self, coalescent, duplication, transfer, loss, 
        hemiplasy, recombination):
        if not coalescent:
            raise IxDTLError('missing coalescent parameter')
        self.__parameters['coalescent'] = coalescent

        if not duplication:
            raise IxDTLError('missing duplication parameter')
        self.__parameters['duplication'] = duplication

        if not transfer:
            raise IxDTLError('missing transfer parameter')
        self.__parameters['transfer'] = transfer

        if not loss:
            raise IxDTLError('missing loss parameter')
        self.__parameters['loss'] = loss

        if hemiplasy is None:
            raise IxDTLError('missing hemiplasy option')
        self.__parameters['hemiplasy'] = hemiplasy

        if recombination is None:
            raise IxDTLError('missing recombination option')
        self.__parameters['recombination'] = recombination

    def readSpeciesTree(self, path):
        self.__speciesTree = SpeciesTree(self.__randomState)
        self.__speciesTree.readFromNewickFile(path)
        self.__speciesTree.setLambdaCoalescent(
            parameter=self.__parameters['coalescent'])
            
        print('species tree:')
        print(self.__speciesTree)
        print()
            
    def constructOriginalHaplotypeTree(self):
        self.__haplotypeTree = HaplotypeTree(self.__randomState)
        self.__haplotypeTree.initialize(locusTree=self.__speciesTree)
        self.__haplotypeTree.setEventRates(
            duplicationRate=self.__parameters['duplication'],
            transferRate=self.__parameters['transfer'],
            lossRate=self.__parameters['loss'])
        self.__haplotypeTree.setRecombination(
            recombination=self.__parameters['recombination'])
        self.__haplotypeTree.setHemiplasy(
            hemiplasy=self.__parameters['hemiplasy'])

        print('original haplotype tree:')
        print(self.__haplotypeTree)
        print()
        
        
