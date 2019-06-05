import json

"""
============= Configuration =============
"""
# Configuration file location
configfile: "config.yaml"

# Load samples's analysis template
SAMPLES = json.load(open(config["SAMPLES"]))

# Determine if an INPUT should be use while performing the peak calling step.
CONTROL = config["CONTROL"]

TARGETS = []
"""
============= Files name =============
"""
for CASE in SAMPLES["case:control"]:

	if CONTROL == "Control":
		TARGETS.append("04peaks/{case}_vs_{control}_Control_summits.bed".format(case=CASE, control=SAMPLES["case:control"][CASE]))


rule all:
	input: TARGETS


"""
============= PeakCalling =============
"""
# Perform peakcalling with macs2 for a IP file WITH its proper Input file
rule call_peaks_macs2:
	input:
		case="02aln/{case}_sorted.bam",
		control="02aln/{control}_sorted.bam"
	output:
		bed="04peaks/{case}_vs_{control}_Control_summits.bed"
	params:
		name="{case}_vs_{control}_Control"
	shell:
		"macs2 callpeak -t {input.case} -c {input.control} -f BAM -g hs -n {params.name} --outdir 04peaks/"
