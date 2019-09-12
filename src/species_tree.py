import skbio
from .tree_table import *


class SpeciesTree:
    def __init__(self):
        pass
    
    def readNewickFile(self, path):
        treeTable = TreeTable()
        self.skbio_tree = treeTable.createFromNewickFile(path)

        print(treeTable.table)