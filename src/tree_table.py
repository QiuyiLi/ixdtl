import skbio


class TreeTable:
    def __init__(self):
        self.__table = []

    @property
    def table(self):
        return self.__table

    def createFromSkbioTree(self, skbioTree):
        root = self.__parse(skbioTree)
        root = self.__rename(root)
        nodes = self.__toList(root, root=root['skbio.TreeNode'])

        for i in range(len(nodes)):
            parentIdx = None
            for j in range(len(nodes)):
                if (nodes[i]['parent'] is nodes[j]['skbio.TreeNode']):
                    parentIdx = j
            entry = {
                'id': i,
                'name': nodes[i]['name'],
                'parent': parentIdx,
                'd2p': nodes[i]['distance'],
            }
            self.__table.append(entry)

        return skbioTree

    def createFromNewickFile(self, path):
        f = open(path)
        skbioTree = skbio.read(f, format="newick", into=skbio.tree.TreeNode)
        f.close()

        return self.createFromSkbioTree(skbioTree)

    def __parse(self, skbioTree):
        currentNode = {
            'skbio.TreeNode': skbioTree,
            'name': skbioTree.name,
            'parent': skbioTree.parent,
            'children': [],
            'distance': skbioTree.length
        }
        if skbioTree.is_tip():
            return currentNode
        for children in skbioTree.children:
            currentNode['children'].append(self.__parse(children))
        return currentNode

    def __rename(self, node):
        currentNode = node.copy()
        name = ''
        for i in range(len(currentNode['children'])):
            currentNode['children'][i] = self.__rename(currentNode['children'][i])
            name += currentNode['children'][i]['name']
        if not currentNode['name']:
            currentNode['name'] = name
        return currentNode

    def __toList(self, node, root):
        currentNode = node.copy()
        nodes = []
        for i in range(len(currentNode['children'])):
            nodes += self.__toList(currentNode['children'][i], root=root)
        del currentNode['children']
        currentNode['distance_to_root'] = currentNode['skbio.TreeNode'].distance(root)
        nodes.append(currentNode)
        if currentNode['skbio.TreeNode'] is root:
            nodes = sorted(nodes, key=lambda x: x['distance_to_root'], reverse=True)
        return nodes