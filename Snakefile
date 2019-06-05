import json

# Configuration file location
configfile: "config.yaml"

# Load samples's analysis template
SAMPLES = json.load(open(config["SAMPLES"]))


TARGETS = []

ALL_BAM = []  # For BAM files
ALL_SORTED_BAM = []  # For sorted BAM files
ALL_BIGWIG = []  # For bigwig files
ALL_PEAKS = []  # For MACS2 peaks
ALL_STAT = []  # For alignments statistics
ALL_ANNOTATE = []  # For annotated files (ROSE_geneMapper output, ...)
ALL_RNA_DIFF = []  # DEG analysis results
ALL_GENE_LIST = []  # Final gene list
#ALL_HOMER = []  # HOMER not implemented yet
ALL_ROSE = []  # ROSE_main.py output
ALL_GENES = []
#ALL_SICER = []  # SICER is no longer implemented in the pipeline

# Determine if an INPUT should be use while performing the peak calling step.
CONTROL = config["CONTROL"]


# Create name templates based on SAMPLES["case:control"] data
for CASE in SAMPLES["case:control"]:

	if CONTROL == "Control":
		ALL_BAM.append("01aln/{case}.bam".format(case=SAMPLES["case:control"][CASE]))
		ALL_STAT.append("01stat/{case}_align_stat.txt".format(case=SAMPLES["case:control"][CASE]))
		ALL_SORTED_BAM.append("02aln/{case}_sorted.bam".format(case=SAMPLES["case:control"][CASE]))  # BAM
		ALL_SORTED_BAM.append("02aln/{case}_sorted.bai".format(case=SAMPLES["case:control"][CASE]))  # BAM Index
		ALL_BIGWIG.append("03bw/{case}.bigWig".format(case=SAMPLES["case:control"][CASE]))  # BIGWIG
	elif CONTROL == "NoControl":


	else: None
	#ALL_BAM.append("01aln/{case}.sam".format(case=CASE))
	#ALL_BAM.append("01aln/{case}.bam".format(case=CASE))
	#ALL_STAT.append("01stat/{case}_align_stat.txt".format(case=CASE))
	#ALL_SORTED_BAM.append("02aln/{case}_sorted.bam".format(case=CASE))
	#ALL_SORTED_BAM.append("02aln/{case}_sorted.bai".format(case=CASE))
	#ALL_SORTED_BAM.append("02aln/{case}_sorted.bed".format(case=CASE))
	#ALL_BIGWIG.append("03bw/{case}.bigWig".format(case=CASE))
	#ALL_PEAKS.append("06peaks/{case}_{control}_peaks.xls".format(case=CASE, control=CONTROL))
	#ALL_GENES.append("11diff_binding/{case}the_control_lambda.bdg".format(case=CASE))
	#ALL_GENES.append("11diff_binding/{case}_peaks.xls".format(case=CASE))



for SAMPLE in SAMPLES["diff"]:
	#ALL_PEAKS.append("06peaks/{alt_case}_vs_{cases}_noOverlap_{control}.bed".format(alt_case=SAMPLE, cases=":".join(SAMPLES["diff"][SAMPLE]), control=CONTROL))
	#ALL_ROSE.append("06peaks/{alt_case}_vs_{cases}_noOverlap_{control}_colSwap.bed".format(alt_case=SAMPLE, cases=":".join(SAMPLES["diff"][SAMPLE]), control=CONTROL))
	#ALL_ROSE.append("10ROSE/all_{alt_case}_vs_{cases}_peaks_sorted_merge_colswap_keep_cond1_ENHANCER_TO_GENE.txt".format(alt_case=SAMPLE, cases=":".join(SAMPLES["diff"][SAMPLE]), control=CONTROL))
	#ALL_ROSE.append("12final_output_genes/{alt_case}_vs_{cases}_{control}_gene_list_ROSE.txt".format(alt_case=SAMPLE, cases=":".join(SAMPLES["diff"][SAMPLE]), control=CONTROL))
	#ALL_GENES.append("11diff_binding/all_{alt_case}_vs_{cases}_peaks_sorted_merge_colswap_cond1.bed".format(alt_case=SAMPLE, cases=":".join(SAMPLES["diff"][SAMPLE])))
	#ALL_GENES.append("11diff_binding/all_{alt_case}_vs_{cases}_peaks_sorted_merge_colswap_keep_cond1.bed".format(alt_case=SAMPLE, cases=":".join(SAMPLES["diff"][SAMPLE])))
	#"11diff_binding/all_{alt_case}_vs_{case1}:{case2}:{case3}_peaks_sorted_merge_colswap_keep_cond1.bed"
	#12final_output_genes/{alt_case}_vs_{case1}:{case2}:{case3}_{control}_gene_list_ROSE.txt
	#ALL_HOMER.append("09homer/{alt_case}_vs_{cases}_noOverlap_{control}.txt".format(alt_case=SAMPLE, cases=":".join(SAMPLES["diff"][SAMPLE]), control=CONTROL))
	for COND in SAMPLES["diff"][SAMPLE]:
		#ALL_SICER.append("02aln/{case}_vs_{cond}_increased.bed".format(case=SAMPLE, cond=COND))
		#ALL_SICER.append("02aln/{case}_vs_{cond}_increased_peaks_annotation.bed".format(case=SAMPLE, cond=COND))
		#ALL_SICER.append("02aln/{case}_sorted-vs-{cond}_sorted-W{window_size}-G{gap_size}-decreased-islands-summary-FDR{FDR}".format(case=SAMPLE, cond=COND, window_size=config["window_size"], gap_size=config["gap_size"], FDR=config["FDR"]))
		#ALL_SICER.append("02aln/{case}_sorted-vs-{cond}_sorted-W{window_size}-G{gap_size}-decreased-islands-summary-FDR{FDR}-annotated".format(case=SAMPLE, cond=COND, window_size=config["window_size"], gap_size=config["gap_size"], FDR=config["FDR"]))
		None
ALL_SICER.append("02aln/{case}_sorted-vs-{cond}_sorted-W{window_size}-G{gap_size}-increased-islands-summary-FDR{FDR}".format(case="IP1_Ac_R1", cond="IP2_Ac_R1", window_size=config["window_size"], gap_size=config["gap_size"], FDR=config["FDR"]))
ALL_SICER.append("02aln/{case}_sorted-vs-{cond}_sorted-W{window_size}-G{gap_size}-increased-islands-summary-FDR{FDR}-annotated".format(case="IP1_Ac_R1", cond="IP2_Ac_R1", window_size=config["window_size"], gap_size=config["gap_size"], FDR=config["FDR"]))



"""
		ALL_ANNOTATE.append("06_2annotate/{alt_case}_vs_{case}_noOverlap_{control}_annotation_{regulation}.bed".format(alt_case=ALT_CASE, case=CASE, control=CONTROL))
		ALL_RNA_DIFF.append("07diff/{alt_case}_vs_{case}_{control}_significiant_genes_rnaseq_{regulation}.tsv".format(alt_case=ALT_CASE, case=CASE, control=CONTROL))
		ALL_GENE_LIST.append("08summary/{alt_case}_vs_{case}_{control}_gene_list.txt".format(alt_case=ALT_CASE, case=CASE, control=CONTROL))
"""

# Add all generated names to TARGETS
TARGETS.extend(ALL_BAM)
TARGETS.extend(ALL_STAT)
TARGETS.extend(ALL_SORTED_BAM)
TARGETS.extend(ALL_BIGWIG)
TARGETS.extend(ALL_PEAKS)
TARGETS.extend(ALL_GENES)
#TARGETS.extend(ALL_ANNOTATE)
#TARGETS.extend(ALL_RNA_DIFF)
TARGETS.extend(ALL_GENE_LIST)
TARGETS.extend(ALL_HOMER)
TARGETS.extend(ALL_ROSE)
TARGETS.extend(ALL_SICER)

#print(TARGETS)

# This rule tell snakemake wich files to create, based on samples.json data
rule all:
	input: TARGETS

#=========================================================
#=====================   Alignment   =====================
#=========================================================

# Map fastq data on hg19 with bowtie2. Return temporary sam files + alignment statistics
rule bowtie2_map:
	input:
		fastq="	Fastq_files/{sample}.fastq"
	params:
		index="/home/data/pameslin/data/ChIP_seq_projet/data/Index_hg19/genome",
		processors= config["bowtie2_processors"]
	output:
		sam=temp("01aln/{sample}.sam"),
		stat="01stat/{sample}_align_stat.txt"

	shell:
		"bowtie2 -p {params.processors} -x {params.index} -U {input.fastq} -S {output.sam} 2> {output.stat}"


# Translate a sam file in sorted bam files with samtools view and sort tools
rule samtools_view:
	input:
		"01aln/{sample}.sam"
	output:
		"02aln/{sample}_sorted.bam"
	shell:
		"samtools view -Sb {input} | samtools sort -@ 20 -o {output} - "


# Create a index file (bai) from a given sorted bam file with samtools index tool
rule samtools_index:
	input:
		"02aln/{sample}_sorted.bam"
	output:
		"02aln/{sample}_sorted.bai"
	shell:
		"samtools index {input} 02aln/{wildcards.sample}_sorted.bai"


# Create a bigWig file from a sorted bam and its index file (bai). Use bamCoverage from Deeptools
rule bamCoverage_bigWig:
	input:
		bam="02aln/{sample}_sorted.bam",
		bai="02aln/{sample}_sorted.bai"
	params:
		processors = config["bamCoverage_processors"]
	output:
		"03bw/{sample}.bigWig"
	shell:
		"source activate deeptools &&"
		"bamCoverage -b {input.bam} -o {output} -of bigwig -e 200 --normalizeUsing RPKM -p {params.processors} --smoothLength 150 --ignoreDuplicates &&"
		"source deactivate deeptools"


# BAM to BED conversion
rule bam_to_bed:
	input:
		"02aln/{sample}_sorted.bam"
	output:
		"02aln/{sample}_sorted.bed"
	shell:
		"bamToBed -i {input}>{output}"


#=========================================================
#===================   Peak calling   ====================
#=========================================================

# Perform peakcalling with macs2 for a IP file WITHOUT its Input file
rule call_peaks_macs2_nocontrol:
	input:
		case="02aln/{case}_sorted.bam",
	output:
		bed="06peaks/{case}_NoControl_summits.bed",
		xls="06peaks/{case}_NoControl_peaks.xls",
		bdg="06peaks/{case}_treat_pileup.bdg"
	params:
		name="{case}_NoControl"
	shell:
		"macs2 callpeak -B -t {input.case} -f BAM -g hs -n {params.name} --outdir 06peaks/"


#=========================================================
#==================  Annotate peaks   ====================
#=========================================================

# Perform columns swap needed for ROSE_geneMapper.py tool using awk
rule awk_bed_formating:
	input:
		"06peaks/{alt_case}_vs_{case}_noOverlap_{control}.bed"
	output:
		"06peaks/{alt_case}_vs_{case}_noOverlap_{control}_colSwap.bed"
	shell:
		"./awk_formating.sh {input} {output}"


# Use ROSE_geneMapper.py script to call genes overlapping or close to a ChIP-seq peak, plus the closest one.
# python /usr/local/bin/ROSE_geneMapper.py -i tmp/all_IP2_Ac_R1_peaks_vs_IP1_Ac_R1_peaks:IP3_Ac_R1_peaks:IP4_Ac_R1_peaks_sorted_merge_colswap_c3.0_cond1.bed -g HG19 -o tmp/
rule ROSE_geneMapper:
	input:
		"11diff_binding/all_{alt_case}_vs_{case1}:{case2}:{case3}_peaks_sorted_merge_colswap_keep_cond1.bed"
	output:
		"10ROSE/all_{alt_case}_vs_{case1}:{case2}:{case3}_peaks_sorted_merge_colswap_keep_cond1_ENHANCER_TO_GENE.txt"
	shell:
		"source activate ROSE_env &&"
		"python /usr/local/bin/ROSE_geneMapper.py -i {input} -g HG19 -o 10ROSE/ &&"
		"source deactivate ROSE_env"


#=========================================================
#===================   DEG results   ====================
#=========================================================

# Return significantly up-regulated or down-regulated genes between {case} and {alt_case}, up-regulated/down-regulated	 in alt_case so.
rule cummeRbund_regulated_genes:
	input:
		"06peaks/{alt_case}_vs_{case}_noOverlap_{control}.bed"
	output:
		"07diff/{alt_case}_vs_{case}_{control}_significiant_genes_rnaseq.tsv"
	shell:
		"Rscript ../ChIP_seq_projet/callUpRegulatedGenes.R {input} {output} {wildcards.alt_case} {wildcards.case}"




"""
# No longer implemented
rule peakcalling_SICER:
	input:
		bed="02aln/{case}_sorted.bed",
		condition="02aln/{cond}_sorted.bed"
	params:
		FDR=config["FDR"]
	output:
		"02aln/{case}_sorted-vs-{cond}_sorted-W{window_size}-G{gap_size}-increased-islands-summary-FDR0.05"#{params.FDR}"
	shell:
		"source activate SICER_env &&"
		"SICER-df-rb.sh {wildcards.case}_sorted.bed {wildcards.cond}_sorted.bed {wildcards.window_size} {wildcards.gap_size} 100 0.05 &&"#"{wildcards.FDR} &&"
		"source deactivate SICER_env"


# No longer implemented
rule annotate_peaks_bedops:
	input:
		"02aln/{case}_sorted-vs-{cond}_sorted-W{window_size}-G{gap_size}-increased-islands-summary-FDR0.05"#{FDR}"
	output:
		"02aln/{case}_sorted-vs-{cond}_sorted-W{window_size}-G{gap_size}-increased-islands-summary-FDR{FDR}-annotated"
	shell:
		"closest-features {input} /home/pameslin/data/ChIP_seq_projet/data/hg19_genes_symbol_sorted.bed > {output}"


# Diffenrial peak calling with MACS2
rule create_bdg_macs2:
	input:
		"02aln/{sample}_sorted.bam"
	output:
		ctr="11diff_binding/{sample}_control_lambda.bdg",
		tre="11diff_binding/{sample}_treat_pileup.bdg",
		xls="11diff_binding/{sample}_peaks.xls"
	shell:
		"macs2 callpeak -B -t {input} -n {wildcards.sample} --nomodel --extsize 200 --outdir 11diff_binding/"


rule call_diff_peaks_macs2:
	input:
		alt_case_treat="11diff_binding/{alt_case}_treat_pileup.bdg",
		alt_case_control="11diff_binding/{alt_case}_control_lambda.bdg",
		alt_case_xls="11diff_binding/{alt_case}_peaks.xls",

		case_treat_1="11diff_binding/{case1}_treat_pileup.bdg",
		case_control_1="11diff_binding/{case1}_control_lambda.bdg",
		case_xls_1="11diff_binding/{case1}_peaks.xls",

		case_treat_2="11diff_binding/{case2}_treat_pileup.bdg",
		case_control_2="11diff_binding/{case2}_control_lambda.bdg",
		case_xls_2="11diff_binding/{case2}_peaks.xls",

		case_treat_3="11diff_binding/{case3}_treat_pileup.bdg",
		case_control_3="11diff_binding/{case3}_control_lambda.bdg",
		case_xls_3="11diff_binding/{case3}_peaks.xls"

	output:
		"11diff_binding/all_{alt_case}_vs_{case1}:{case2}:{case3}_peaks_sorted_merge_colswap_cond1.bed"
	shell:
		"echo {input.case_control_2} &&"
		"./cmd.sh {input.alt_case_xls} {input.case_xls_1} {input.case_xls_2} {input.case_xls_3} {input.alt_case_treat} {input.alt_case_control} {input.case_treat_1} {input.case_control_1} {input.case_treat_2} {input.case_control_2} {input.case_treat_3} {input.case_control_3}"



rule genes_discard:
	input:
		"11diff_binding/all_{alt_case}_vs_{case1}:{case2}:{case3}_peaks_sorted_merge_colswap_cond1.bed"
	output:
		"11diff_binding/all_{alt_case}_vs_{case1}:{case2}:{case3}_peaks_sorted_merge_colswap_keep_cond1.bed"
	shell:
		#"echo {output} &&"
		"python3 onlyThreePeaks.py {input} {output}"


# Return a txt file with genes significantly up-regulated genes in RNA-seq experiment AND overlapping and/or near and the closest one to a peak in ChIP-seq experiment.
rule python3_ChIPseq_vs_RNAseq:
	input:
		tsv="07diff/{alt_case}_vs_{case1}:{case2}:{case3}_{control}_significiant_genes_rnaseq.tsv",
		txt="10ROSE/all_{alt_case}_vs_{case1}:{case2}:{case3}_peaks_sorted_merge_colswap_keep_cond1_ENHANCER_TO_GENE.txt"
	output:
		"12final_output_genes/{alt_case}_vs_{case1}:{case2}:{case3}_{control}_gene_list_ROSE.txt"
	shell:
		"python3 callGenes_RNA_ChIP_seq.py {input.tsv} {input.txt} {output}"


# Return a bed file with features present in IP2 and NOT in IP1. Use intersect tool from bedtools
rule bedtools_intersect	:
	input:
		alt_case="06peaks/{alt_case}_{control}_summits.bed",
		case1="06peaks/{bed1}_{control}_summits.bed",
		case2="06peaks/{bed2}_{control}_summits.bed",
		case3="06peaks/{bed3}_{control}_summits.bed"
	output:
		"06peaks/{alt_case}_vs_{bed1}:{bed2}:{bed3}_noOverlap_{control}.bed"
	shell:
		"bedtools intersect -wa -a {input.alt_case} -b {input.case1} -v -f 0.50 -r | bedtools intersect -wa -a stdin -b {input.case2} -v -f 0.50 -r | bedtools intersect -wa -a stdin -b {input.case3} -v -f 0.50 -r > {output}"

# Use HOMER tool to annotate a given bed file.
rule homer_annotatePeaks:
	input:
		"06peaks/{alt_case}_vs_{case}_noOverlap_{control}.bed"
	output:
		file="09homer/{alt_case}_vs_{case}_noOverlap_{control}.txt",
		dir=directory("09homer/{alt_case}_vs_{case}_noOverlap_{control}")
	shell:
		"source activate homer &&"
		"annotatePeaks.pl {input} hg19 -go {output.dir} > {output.file} &&"
		"source deactivate homer"

# Perform peakcalling with macs2 for a IP file WITH its proper Input file
rule call_peaks_macs2:
	input:
		case="02aln/{case}_sorted.bam",
		control="02aln/{control}_sorted.bam"
	output:
		bed="04peaks/{case}_vs_{control}_Control_summits.bed"
	params:
		name="{case}_vs_{control}"
	shell:
		"macs2 callpeak -t {input.case} -c {input.control} -f BAM -g hs -n {params.name} --outdir 04peaks/"



# Not used
rule bedtools_closest_features:
	input:
		bed="06peaks/{alt_case}_vs_{case}_noOverlap_{control}.bed",
		annotation="../ChIP_seq_projet/data/hg19_genes_symbol_sorted.bed"
	output:
		"06_2annotate/{alt_case}_vs_{case}_noOverlap_{control}_annotation.bed"
	shell:
		"bedtools closest -a {input.bed} -b {input.annotation} | cut -d$'\t' -f1,2,3,4,9 > {output}"


# Old script, see next rule.
rule python_gene_name:
	input:
		tsv="07diff/{alt_case}_vs_{case}_{control}_significiant_genes_rnaseq.tsv",
		bed="06_2annotate/{alt_case}_vs_{case}_noOverlap_{control}_annotation.bed"
	output:
		"08summary/{alt_case}_vs_{case}_{control}_gene_list.txt"
	shell:
		"python3 ../RNA-seq/main.py {input.tsv} {input.bed} {output}"

# Not yet implemented
rule call_peaks_macs2_broad_nomodel:
	input:
		case="02aln/{case}_sorted.bam",
		control="02aln/{control}_sorted.bam"
	output:
		bed="05peaks/{case}_vs_{control}_nomodel_broad_peaks.xls"
	params:
		name="{case}_vs_{control}_nomodel_broad"
	shell:
		"macs2 callpeak -t {input.case} -c {input.control} -f BAM -g hs -n {params.name} --outdir 05peaks/ --nomodel --broad"
"""
