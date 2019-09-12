import skbio


class SpeciesTree:
    def __init__(self):
        pass
    
    def readNewick(self, path):
        self.skbio_tree = self.newickToTable(path)

    def newickToTable(self, path):
        f = open(path)
        tree = skbio.read(f, format="newick", into=skbio.tree.TreeNode)
        f.close()

        def parse(tree):
            node = {
                'skbio.TreeNode': tree,
                'name': tree.name,
                'parent': tree.parent,
                'children': [],
                'distance': tree.length
            }
            if tree.is_tip():
                return node
            for children in tree.children:
                node['children'].append(parse(children))
            return node

        def rename(node):
            ret = node.copy()
            name = ''
            for i in range(len(ret['children'])):
                ret['children'][i] = rename(ret['children'][i])
                name += ret['children'][i]['name']
            if (not ret['name']):
                ret['name'] = name
            return ret

        def toList(node, root):
            d = node.copy()
            ret = []
            for i in range(len(d['children'])):
                ret += toList(d['children'][i], root=root)
            del d['children']
            d['distance_to_root'] = d['skbio.TreeNode'].distance(root)
            ret.append(d)
            if (d['skbio.TreeNode'] is root):
                ret = sorted(ret, key=lambda x: x['distance_to_root'], reverse=True)
            return ret

        root = parse(tree)
        root = rename(root)
        nodes = toList(root, root=root['skbio.TreeNode'])

        # output to file
        # f = open(output_path, 'w')
        # f.write('id\tname\tparent\td2p\n')
        # for i in range(len(nodes)):
        #     parent_index = 'None'
        #     for j in range(len(nodes)):
        #         if (nodes[i]['parent'] is nodes[j]['skbio.TreeNode']):
        #             parent_index = str(j)
        #     f.write(str(i) + '\t' + nodes[i]['name'] + '\t' + parent_index + '\t' + str(nodes[i]['distance']) + '\n')
        # f.close()

        print(nodes)

        return tree