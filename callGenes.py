import sys


def get_genes(input_file):
    """
        Input: a bed file from SICER peakcalling annotate with closest-features (bedops)
        Output: a list of genes present in the file
    """
    opening = open(input_file, "r")
    lines = opening.readlines()
    genes = []
    for line in lines:
        #print(len(line.split("\t")))
        split = line.split("\t")
        length = len(split)
        if length == 24:
            genes.append(split[15])
        elif length == 35:
            genes.append(split[15])
            genes.append(split[26])
    #print(len(list(set(genes))))
    #print("\n".join(list(set(genes))))
    return list(set(genes))


def count_occ(conditions):
    genes = {}
    for condition_file in conditions:
        genes_list = get_genes(condition_file)
        for gene in genes_list:
            not_in_list = 1
            for key in genes.keys():
                if (key == gene):
                    genes[gene] += 1
                    not_in_list = 0
            if not_in_list:
                genes[gene] = 1
    #print(len(genes))
    return genes


def discard_genes(genes_dict):
    gene_set = [gene for gene in genes_dict if genes_dict[gene] == len(sys.argv)-2]
    return gene_set


def rna_genes(rna_file):
    opening = open(rna_file, "r")
    lines = opening.readlines()
    genes = [line.split("\t")[1][:-1] for line in lines]
    return genes



def overlap_rna_chip_genes(rna_genes, chip_genes):
    genes = []
    #not_in_rna = {}
    #for chip_gene in chip_genes:
        #not_in_rna[chip_gene] = 0
    for rna_gene in rna_genes:
        #occ = 0
        for chip_gene in chip_genes:
            if rna_gene == chip_gene:
                genes.append(rna_gene)
                #occ += 1
                #not_in_rna[chip_gene] += 1
        #if occ == 0:
            #print(rna_gene)
    """
    for chip_gene in not_in_rna:
        if not_in_rna[chip_gene] == 0:
            print(chip_gene)
    """
    return genes


def main():
    print("Retourne les gènes présents dans toutes les listes données en entrée, autre la liste RNA-seq. Retounr la longeur de la liste.")
    conditions = sys.argv[2:]
    genes_dict = count_occ(conditions)
    gene_set = discard_genes(genes_dict)
    print(len(gene_set))
    #rna_genes_list = rna_genes(sys.argv[1])
    #overlap_genes = overlap_rna_chip_genes(rna_genes_list, gene_set)
    #return overlap_genes


#out = list(set(main))
out = main()

#print("\n".join(out))
