import skbio
from .tree_table import *


class SpeciesTree:
    def __init__(self):
        self.__skbioTree = None
        self.__treeTable = None

    @property
    def skbioTree(self):
        return self.__skbioTree

    @property
    def treeTable(self):
        return self.__treeTable
    
    def readNewickFile(self, path):
        self.__treeTable = TreeTable()
        self.__skbioTree = self.__treeTable.createFromNewickFile(path)

        print(self.__treeTable.table)