from cuteSV_AVLTree import TreeNode, AVLTree, Record
from pysam import VariantFile
from multiprocessing import Pool, Manager
import vcf
import sys


def test():
    vcf_reader = VariantFile('test1.vcf', 'r')
    for contig in vcf_reader.header.contigs:
        print(contig)


def solve_chrom(vcf_filenames, chrom):
    tree = AVLTree()
    for i in range(len(vcf_filenames)):
        vcf_reader = VariantFile(vcf_filenames[i], 'r')
        for record in vcf_reader.fetch(chrom):
            tree.insert(i, record)

    output = open(chrom + '_tmp.vcf', 'w')
    node_list = []
    tree.inorder(tree.root, node_list)
    for node in node_list:
        output.write(node)
        output.write('\n')



def main(argv):
    vcf_filenames = []
    with open(argv[0], 'r') as f:
        for line in f:
            vcf_filenames.append(line.strip())
    chrom_set = set()
    # 默认header的config中存放了所有chrom信息
    for vcf_filename in vcf_filenames:
        vcf_reader = VariantFile(vcf_filename, 'r')
        for contig in vcf_reader.header.contigs:
            chrom_set.add(contig)  # what is contig
    pool = Pool(processes = int(argv[1]))
    for chrom in chrom_set:
        # multi processes
        pool.apply_async(solve_chrom, args=(vcf_filenames, chrom))
        '''
        for record in vcf_reader:
            if record.CHROM in chrom_dict:
                chrom_dict[record.CHROM].append(record)
            else:
                chrom_dict[record.CHROM] = [record]
        '''


if __name__ == '__main__':
    #main(sys.argv[1:])  # filenames.txt, threads
    test()
