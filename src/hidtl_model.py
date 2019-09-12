from .species_tree import *


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

    @locusTrees.setter
    def locusTrees(self, locusTrees):
        self.__locusTrees = locusTrees

    @property
    def parameter(self):
        return self.__parameter

    @parameter.setter
    def parameter(self, parameter):
        self.__parameter = parameter

    def readSpeciesTree(self, path):
        self.__speciesTree = SpeciesTree()
        self.__speciesTree.readNewick(path)