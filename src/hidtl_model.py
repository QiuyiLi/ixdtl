from .species_tree import *
from .exception import *


class HIDTLModel:
    def __init__(self):
        self.__speciesTree = None
        self.__haplotypeTree = None
        self.__locusTrees = []
        self.__parameter = {}
    
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
    def parameter(self):
        return self.__parameter

    def setParameters(self, coalescent, duplication, transfer, loss, hemiplasy, recombination):
        if not coalescent: 
            raise HIDTLError('missing coalescent parameter')
        self.__parameter['coalescent'] = coalescent

        if not duplication: 
            raise HIDTLError('missing duplication parameter')
        self.__parameter['duplication'] = duplication

        if not transfer: 
            raise HIDTLError('missing transfer parameter')
        self.__parameter['transfer'] = transfer

        if not loss:
            raise HIDTLError('missing loss parameter')
        self.__parameter['loss'] = loss
        
        if hemiplasy is None: 
            raise HIDTLError('missing hemiplasy option')
        self.__parameter['hemiplasy'] = hemiplasy

        if recombination is None: 
            raise HIDTLError('missing recombination option')
        self.__parameter['recombination'] = recombination

    def readSpeciesTree(self, path):
        self.__speciesTree = SpeciesTree()
        self.__speciesTree.readFromNewickFile(path)
        self.__speciesTree.parameter = self.__parameter
        print(self.__speciesTree.parameter['coalescent'])

    def createHaplotypeTree(self):
        coalescent_process, clade_set_into_root = self.__speciesTree.coalescent(distance_above_root=10000)
        print(coalescent_process)