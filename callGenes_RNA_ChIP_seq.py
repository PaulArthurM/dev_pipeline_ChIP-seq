import sys

start_msg = """ Inputs : a gene list from DEG RNA-seq data, a gene list from ChIP-seq data and a path where save results.
 Output: Overlapping genes
 > ChIP-seq data must come from ROSE_geneMapper"""

print(start_msg)

def gene_list_rna():
	opening_rna = open(sys.argv[1])

	lines_rna = opening_rna.readlines()
	genes_rna = [line.split("\t")[1][:-1] for line in lines_rna if line.split("\t")[1][:-1] not in lines_rna]

	return genes_rna


def gene_list_chip():
	opening_chip = open(sys.argv[2])

	lines_chip = opening_chip.readlines()

	gene_name_chip = []

	cols = [4, 5, 6]

	for line in lines_chip:
		for col in cols:
			gene_names = line.split("\t")[col].split(",")
			for gene_name in gene_names:
				if gene_name not in gene_name_chip:
					gene_name_chip.append(gene_name)

	return gene_name_chip


path = sys.argv[3]


chip_list = gene_list_chip()
rna_list = gene_list_rna()


file_to_write = open(path, "w")

genes_write = []
for gene_name_rna in rna_list:
	for gene_name_chip in chip_list:
		if (gene_name_rna == gene_name_chip) and (not gene_name_rna in genes_write):
			genes_write.append(gene_name_rna)
			gene_name = gene_name_rna + "\n"
			file_to_write.write(gene_name)
file_to_write.close()
