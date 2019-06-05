#!/usr/bin/zsh

# Arguments: BED file of super-enhancers and a BED file of black-listed regions by ENCODE.
SUPER_BED=$1
BLACK_LISTED="/home/data/pameslin/ENCFF001TDO.bed"
OUTPUT=$2
bedtools intersect -a $SUPER_BED -b $BLACK_LISTED -v > $OUTPUT
