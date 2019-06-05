#!/usr/bin/zsh

# Variables
WHITE_LISTED_SE=$1
RANKED_SE=$WHITE_LISTED_SE"_RANKED"

IP_SLOPED=$2


# Pour ranker les enhancers par rang dans le fichier bed:
sort -t$'\t' -k5 -n $WHITE_LISTED_SE.bed > $RANKED_SE.bed

# Calcul du coverage pour chaque features:
bedtools coverage -a $RANKED_SE.bed -b $IP_SLOPED".txt" > $IP_SLOPED"_processed.txt"
