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
        
        if hemiplasy != 0 and hemiplasy != 1: 
            raise HIDTLError('missing hemiplasy option')
        self.__parameter['hemiplasy'] = hemiplasy

        if recombination != 0 and recombination != 1: 
            raise HIDTLError('missing recombination option')
        self.__parameter['recombination'] = recombination

    def readSpeciesTree(self, path):
        self.__speciesTree = SpeciesTree()
        self.__speciesTree.readNewickFile(path)

    def createHaplotypeTree(self):
        pass