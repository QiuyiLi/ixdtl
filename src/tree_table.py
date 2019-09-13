import skbio

class TreeTableEntry:
    def __init__(self):
        self.__id = None
        self.__name = None
        self.__parent = None
        self.__distanceToParent = None
        self.__children = []
        self.__distanceToChildren = []

    def __repr__(self):
        return (f'<TreeTableEntry, id: {self.__id}, name: {self.__name}, ' \
                f'parent: {self.__parent}, distanceToParent: {self.__distanceToParent}, ' \
                f'children: {self.__children}, distanceToChildren: {self.__distanceToChildren})')
    
    @property
    def id(self):
        return self.__id
    @id.setter
    def id(self, id):
        self.__id = id

    @property
    def name(self):
        return self.__name
    @name.setter
    def name(self, name):
        self.__name = name

    @property
    def parent(self):
        return self.__parent
    @parent.setter
    def parent(self, parent):
        self.__parent = parent

    @property
    def distanceToParent(self):
        return self.__distanceToParent
    @distanceToParent.setter
    def distanceToParent(self, distanceToParent):
        self.__distanceToParent = distanceToParent

    @property
    def children(self):
        return self.__children
    @children.setter
    def children(self, children):
        self.__children = children

    @property
    def distanceToChildren(self):
        return self.__distanceToChildren
    @distanceToChildren.setter
    def distanceToChildren(self, distanceToChildren):
        self.__distanceToChildren = distanceToChildren


class TreeTable:
    def __init__(self):
        self.__table = []

    @property
    def table(self):
        return self.__table

    def createFromSkbioTree(self, skbioTree):
        self.__renameTreeNodes(skbioTree)
        nodesDict = self.__toDict(skbioTree)
        nodesList = self.__toList(nodesDict, root=nodesDict['skbio.TreeNode'])

        for i in range(len(nodesList)):
            parentIdx = -1
            for j in range(len(nodesList)):
                if nodesList[i]['parent'] is nodesList[j]['skbio.TreeNode']:
                    parentIdx = j
            entry = TreeTableEntry()
            entry.id = i
            entry.name = nodesList[i]['name']
            entry.parent = parentIdx
            entry.distanceToParent = -1.0 if not nodesList[i]['distance'] else nodesList[i]['distance']
            self.__table.append(entry)
        
        print(self.__table)
        for entry in self.__table:
            children, distanceToChildren = self.__getChildrenAndDistances(nodesDict['skbio.TreeNode'], entry.name)
            print(children)
            print(distanceToChildren)
            break

            # for child in children:
            #     node_id = self.nodes_name_dict[child.name].node_id
            #     node.children.append(node_id)

        return skbioTree

    def createFromNewickFile(self, path):
        f = open(path)
        skbioTree = skbio.read(f, format="newick", into=skbio.tree.TreeNode)
        f.close()

        return self.createFromSkbioTree(skbioTree)

    def __renameTreeNodes(self, skbioTree):
        if skbioTree.name:
            return skbioTree.name
        else:
            name = ''
            for child in skbioTree.children:
                name += self.__renameTreeNodes(child)
            skbioTree.name = name
            return skbioTree.name

    def __getChildrenAndDistances(self, skbioTree, name):
        print(skbioTree.children)
        for child in skbioTree.children:
            print(skbioTree.children)
            if child.name == name:
                distances = []
                print(skbioTree.children)
                for child in skbioTree.children:
                    distances.append(skbioTree.distance(child))
                return skbioTree.children, distances

    def __toDict(self, skbioTree):
        nodesDict = {
            'skbio.TreeNode': skbioTree,
            'name': skbioTree.name,
            'parent': skbioTree.parent,
            'children': [],
            'distance': skbioTree.length
        }
        if skbioTree.is_tip():
            return nodesDict
        for children in skbioTree.children:
            nodesDict['children'].append(self.__toDict(children))
        return nodesDict

    def __toList(self, node, root):
        nodesDict = node.copy()
        nodes = []
        for i in range(len(nodesDict['children'])):
            nodes += self.__toList(nodesDict['children'][i], root=root)
        del nodesDict['children']
        nodesDict['distanceToRoot'] = nodesDict['skbio.TreeNode'].distance(root)
        nodes.append(nodesDict)
        if nodesDict['skbio.TreeNode'] is root:
            nodes = sorted(nodes, key=lambda x: x['distanceToRoot'], reverse=True)
        return nodes