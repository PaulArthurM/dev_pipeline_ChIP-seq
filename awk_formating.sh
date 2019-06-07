#!/usr/bin/zsh

# Perform columns swap with awk  
awk 'BEGIN{OFS=FS="\t"} {col1=$1;col2=$2;col3=$3;col4=$4;col5=$5;$1=col4;$2=col1;$3=col2;$4=col3;$5="";print}' $1 > $2
