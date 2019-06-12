import sys


def genes_in_file(file):
    print("Opening {file}".format(file = file))
    genes = open(file, "r").readlines()
    print("{n} genes in ChIP-seq file".format(n = len(genes)))
    return genes


def rna_genes(rna_file):
    print("Opening {file}".format(file = rna_file))
    lines = open(rna_file, "r").readlines()
    genes = [line.split("\t")[1][:-1] for line in lines[1:]]
    print("{n} genes in RNA-seq file".format(n = len(genes)))
    return genes


def intersect_list(chip, rna):
    intersect_list = []
    for cgene in chip:
        for rgene in rna:
            if (cgene[:-1] == rgene):
                if (rgene not in intersect_list):
                    intersect_list.append(rgene)
    print("{n} gene(s) intersected".format(n = len(intersect_list)))
    return intersect_list


def print_genes(list):
    for gene in list:
        print("\t",gene)


def main(chip_seq_genes=sys.argv[1], rna_seq_genes=sys.argv[2]):
    """Take two lists of genes , and return intersected genes.
    
    chip_seq_genes: text file of ChIP-seq genes (first argument of the script)
    rna_seq_genes: text file of RNA-seq genes (second argument of the script)
    
    """
    print_genes(intersect_list(genes_in_file(chip_seq_genes), rna_genes(rna_seq_genes)))


main()
