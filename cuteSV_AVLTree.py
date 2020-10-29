#from functools import total_ordering
import re

class Record(object):
    def __init__(self, record):
        self.name = record.ID
        self.type = record.INFO['SVTYPE']
        self.start = record.POS
        if record.INFO['SVTYPE'] == 'INS' and 'SVLEN' in record.INFO:
            self.end = record.INFO['SVLEN']
        elif record.INFO['SVTYPE'] != 'INS' and 'END' in record.INFO:
            self.end = record.INFO['END']
        elif record.INFO['SVTYPE'] == 'BND' or record.INFO['SVTYPE'] == 'TRA':
            # 若格式不正确可能下标越界
            self.end = re.search('[^N\[\]]+', str(record.ALT[0])).group().split(':')[1]
        self.end = record.INFO['SVLEN'] if record.INFO['SVTYPE'] == 'INS' else record.INFO['END']
        self.chrom1 = record.CHROM
        if record.INFO['SVTYPE'] == 'TRA':
            self.chrom2 = re.search('[^N\[\]]+', str(record.ALT[0])).group().split(':')[0]
        else:
            self.chrom2 = record.CHROM
        if 'STRANDS' in record.INFO:
            self.strand = record.INFO['STRANDS']
        else:
            self.strand = '.'


    def to_string(self):
        return 'id = ' + str(self.id) + ' start: ' + str(self.start)


#@total_ordering
class TreeNode(object):
    def __init__(self, id, record, left=None, right=None):
        self.variant_list = [Record(record)]
        self.start = variant_list[0].start
        self.end = variant_list[0].end
        self.vis = set()
        self.vis.add(id)
        self.left = left
        self.right = right
        self.height = 0
    def add(self, record):
        variant_list.append(record)
    '''
    def __eq__(self, other):
        return self.data.start == other.data.start
    def __lt__(self, other):
        return self.data.start < other.data.start
    '''


'''
    input, record -> Record
    return 0 -> input与record被认为可以merge
    return -1 -> 不可以merge
'''
def check_node(input, record):
    if input.type == record.type and input.strand == record.strand:
        if abs(input.start - record.start) < 1000 and abs(input.end - record.end) < 1000:
            return True, 0
        if abs(input.start - record.start) < 1000:
            return input.start - record.start
        return
    else:
        return False, -1


class AVLTree(object):
    def __init__(self):
        self.root = None

    def insert(self, id, record):
        if not self.root:
            self.root = TreeNode(id, record)
        else:
            self.root = self._insert(id, record, self.root)
    def _insert(self, id, record, node):  # id -> str, record -> Record, node -> TreeNode
        if node is None:
            node = TreeNode(id, record)
        if id not in node.vis:
            check_ans = check_node(node.variant_list[0], record)  # 默认与第一个比较
        if check_ans == 0 and id not in node.vis:  # 可以merge且id不重复
            node.vis.add(id)
            node.add(record)
        elif key < node.data:
            node.left = self._insert(id, record, node.left)
            if (self.height(node.left) - self.height(node.right)) == 2:
                if key < node.left.data:
                    node = self.singleLeftRotate(node)
                else:
                    node = self.doubleLeftRotate(node)
        elif key > node.data:
            node.right = self._insert(id, record, node.right)
            if (self.height(node.right) - self.height(node.left)) == 2:
                if key > node.right.data:
                    node = self.singleRightRotate(node)
                else:
                    node = self.doubleRightRotate(node)
        node.height = max(self.height(node.right), self.height(node.left)) + 1
        return node
    def delete(self, key):
        if self.root is None:
            raise KeyError('Error,empty tree')
        else:
            self.root = self._delete(key, self.root)
    def _delete(self, key, node):
        if node is None:
            raise KeyError('Error,key not in tree')
        elif key < node.data:
            node.left = self._delete(key, node.left)
            if (self.height(node.right) - self.height(node.left)) == 2:
                if self.height(node.right.right) >= self.height(node.right.left):
                    node = self.singleRightRotate(node)
                else:
                    node = self.doubleRightRotate(node)
            node.height = max(self.height(node.left), self.height(node.right)) + 1
        elif key > node.data:
            node.right = self._delete(key, node.right)
            if (self.height(node.left) - self.height(node.right)) == 2:
                if self.height(node.left.left) >= self.height(node.left.right):
                    node = self.singleLeftRotate(node)
                else:
                    node = self.doubleLeftRotate(node)
            node.height = max(self.height(node.left), self.height(node.right)) + 1
        elif node.left and node.right:
            if node.left.height <= node.right.height:
                minNode = self._findMin(node.right)
                node.key = minNode.key
                node.right = self._delete(node.key, node.right)
            else:
                maxNode = self._findMax(node.left)
                node.key = maxNode.key
                node.left = self._delete(node.key, node.left)
            node.height = max(self.height(node.left), self.height(node.right)) + 1
        else:
            if node.right:
                node = node.right
            else:
                node = node.left
        return node
    def inorder(self, start_node, node_list):
        if start_node is not None:
            self.inorder(start_node.left, node_list)
            node_list.append((start_node.variant_list, vis))
            #print(start_node.variant_list)
            self.inorder(start_node.right, node_list)

    def find(self, key):
        if not self.root:
            return None
        else:
            return self._find(key, self.root)
    def _find(self, key, node):
        if not node:
            return None
        elif key < node.data:
            return self._find(key, node.left)
        elif key > node.data:
            return self._find(key, node.right)
        else:
            return node
    def findMin(self):
        if self.root is None:
            return None
        else:
            return self._findMin(self.root)
    def _findMin(self, node):
        if node.left:
            return self._findMin(node.left)
        else:
            return node
    def findMax(self):
        if self.root is None:
            return None
        else:
            return self._findMax(self.root)
    def _findMax(self, node):
        if node.right:
            return self._findMax(node.right)
        else:
            return node
    def height(self, node):
        if node is None:
            return -1
        else:
            return node.height
    #在node节点的左孩子k1的左子树添加了新节点，左旋转
    def singleLeftRotate(self, node):
        k1 = node.left
        node.left = k1.right
        k1.right = node
        node.height = max(self.height(node.right), self.height(node.left)) + 1
        k1.height = max(self.height(k1.left), node.height) + 1
        return k1
    #在node节点的右孩子k1的右子树添加了新节点，右旋转
    def singleRightRotate(self, node):
        k1 = node.right
        node.right = k1.left
        k1.left = node
        node.height = max(self.height(node.right), self.height(node.left)) + 1
        k1.height = max(self.height(k1.right), node.height) + 1
        return k1
    #在node节点的左孩子的右子树添加了新节点，先左后右
    def doubleRightRotate(self, node):
        node.right = self.singleLeftRotate(node.right)
        return self.singleRightRotate(node)
    #在node节点的右孩子的左子树添加了新节点,先右后左
    def doubleLeftRotate(self, node):
        node.left = self.singleRightRotate(node.left)
        return self.singleLeftRotate(node)
