import numpy as np
from .species_tree import *
from .exception import *


class HIDTLModel:
    def __init__(self):
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

    def run(self, inputFile, coalescentArgs, duplicationArgs, transferArgs, lossArgs,
            hemiplasy, recombination):
        self.setParameters(coalescent=coalescentArgs, 
                           duplication=duplicationArgs, 
                           transfer=transferArgs, 
                           loss=lossArgs,
                           hemiplasy=hemiplasy,
                           recombination=recombination)

        self.readSpeciesTree(inputFile)

        self.createHaplotypeTree()

    def setParameters(self, coalescent, duplication, transfer, loss, hemiplasy, recombination):
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

    def readSpeciesTree(self, path):
        self.__speciesTree = SpeciesTree()
        self.__speciesTree.readFromNewickFile(path)


        self.__speciesTree.lambdaCoalescent = \
            np.random.gamma(shape=self.__parameters['coalescent']['shape'], 
                            scale=self.__parameters['coalescent']['scale'], 
                            size=len(self.__speciesTree.getLeaves()))
        

    def createHaplotypeTree(self):
        coalescent_process, clade_set_into_root = self.__speciesTree.coalescent(distance_above_root=10000)
        print(coalescent_process)