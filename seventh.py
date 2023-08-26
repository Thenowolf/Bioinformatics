class BinaryTreeNode:
    def __init__(self, data):
        self.data = data
        self.leftChild = None
        self.rightChild = None


treedata = ["A", "BC", "DE"]

def fitchAlg(tree):
    for i in range(len(tree)):
        #print(tree[i])
        #pomoc = i
        if len(tree[i]) == 2:
            #print(set(tree[i][0]).intersection(set(tree[i][1])))
            if(set(tree[i][0]).intersection(set(tree[i][1])) == set()):
                node = set(tree[i][0]).intersection(set(tree[i][1]))
            else:
                node = set(tree[i][0]).union(set(tree[i][1]))
                #print(1)

            #print(1)
            #print(tree[i][y])
    
fitchAlg(treedata)
