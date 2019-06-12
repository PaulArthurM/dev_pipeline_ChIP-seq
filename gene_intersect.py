import sys


def genes_in_file(file):
    """"Take a gene list file and return a python list of all genes."""
    print("Opening {file}".format(file = file))
    lines = open(file, "r").readlines()
    genes = [line[:-1] for line in lines]
    print("{n} genes in ChIP-seq file".format(n = len(genes)))
    return genes


def rna_genes(rna_file):
    """Take tsv file from DEG analysis and return all genes."""
    print("Opening {file}".format(file = rna_file))
    lines = open(rna_file, "r").readlines()
    genes = [line.split("\t")[1][:-1] for line in lines[1:]]
    print("{n} genes in RNA-seq file".format(n = len(genes)))
    return genes


def intersect_list(chip, rna):
    """Take two lists of genes and return commons genes."""
    intersect_list = []
    for cgene in chip:
        for rgene in rna:
            if (cgene == rgene):
                if (rgene not in intersect_list):
                    intersect_list.append(rgene)
    print("{n} gene(s) intersected".format(n = len(intersect_list)))
    return intersect_list


def print_genes(list):
    """Take a list of genes and print them all."""
    for gene in list:
        print("\t",gene)


def main(chip_seq_genes=sys.argv[1], rna_seq_genes=sys.argv[2]):
    """Take two lists of genes , and return intersected genes.

    chip_seq_genes: text file of ChIP-seq genes (first argument of the script)
    rna_seq_genes: text file of RNA-seq genes (second argument of the script)

    """
    rna_seq_gene_list = rna_genes(rna_seq_genes)
    chip_seq_genes_list = genes_in_file(chip_seq_genes)
    intersected_gene_list = intersect_list(chip_seq_genes_list, rna_seq_gene_list)
    print_genes(intersected_gene_list)


if __name__ == '__main__':
    main()
