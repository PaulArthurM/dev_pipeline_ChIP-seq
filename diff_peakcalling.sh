
echo "This script was for diffenrial peakcalling of ChIP-seq data"
alt_case="$(grep 'tags after filtering in treatment:' $1 | cut -d' ' -f7)"
case1="$(grep 'tags after filtering in treatment:' $2 | cut -d' ' -f7)"
case2="$(grep 'tags after filtering in treatment:' $3 | cut -d' ' -f7)"
case3="$(grep 'tags after filtering in treatment:' $4 | cut -d' ' -f7)"

alt_treat=$5
alt_contr=$6

treat_1=$7
contr_1=$8

treat_2=$9
#echo $treat_2 
contr_2=${10}
#echo $contr_2
treat_3=${11}
#echo $treat_3
contr_3=${12}
#echo $contr_3


#echo $treat_1
#echo $alt_contr
#echo $treat_1
#echo $contr_1
#echo $alt_case
#echo $case1

base_name_alt="$(echo $1 | cut -d'.' -f1 | cut -d'/' -f2 | sed 's/_peaks//')"
base_name_1="$(echo $2 | cut -d'.' -f1 | cut -d'/' -f2 | sed 's/_peaks//')"
base_name_3="$(echo $3 | cut -d'.' -f1 | cut -d'/' -f2 | sed 's/_peaks//')"
base_name_4="$(echo $4 | cut -d'.' -f1 | cut -d'/' -f2 | sed 's/_peaks//')"



#echo $contr_2

#echo $base_name_alt
#echo $base_name_1
#echo $base_name_3
#echo $base_name_4

prefix_1="diff_${base_name_alt}_vs_${base_name_1}"
prefix_3="diff_${base_name_alt}_vs_${base_name_3}"
prefix_4="diff_${base_name_alt}_vs_${base_name_4}"

#echo $contr_2

cmd1="macs2 bdgdiff --t1 $alt_treat --c1 $alt_contr --t2 $treat_1 --c2 $contr_1 --d1 $alt_case --d2 $case1 -g 60 -l 200 --o-prefix 11diff_binding/${prefix_1}"
cmd2="macs2 bdgdiff --t1 $alt_treat --c1 $alt_contr --t2 $treat_2 --c2 $contr_2 --d1 $alt_case --d2 $case2 -g 60 -l 200 --o-prefix 11diff_binding/${prefix_3}"
cmd3="macs2 bdgdiff --t1 $alt_treat --c1 $alt_contr --t2 $treat_3 --c2 $contr_3 --d1 $alt_case --d2 $case3 -g 60 -l 200 --o-prefix 11diff_binding/${prefix_4}"

#echo $cmd1
#echo $cmd2
#echo $cmd3

$cmd1
$cmd2
$cmd3

#echo "11diff_binding/all_${base_name_alt}_vs_${base_name_1}:${base_name_3}:${base_name_4}_peaks_sorted_merge_colswap_cond1.bed"


for f in 11diff_binding/*_c3.0_cond1.bed; do cat $f >> 11diff_binding/all_${base_name_alt}_vs_${base_name_1}:${base_name_3}:${base_name_4}_cond1.bed; done

bedtools sort -i 11diff_binding/all_${base_name_alt}_vs_${base_name_1}:${base_name_3}:${base_name_4}_cond1.bed > 11diff_binding/all_${base_name_alt}_vs_${base_name_1}:${base_name_3}:${base_name_4}_sorted_cond1.bed

bedtools merge -i 11diff_binding/all_${base_name_alt}_vs_${base_name_1}:${base_name_3}:${base_name_4}_sorted_cond1.bed  -c 4 -o distinct -delim "|" > 11diff_binding/all_${base_name_alt}_vs_${base_name_1}:${base_name_3}:${base_name_4}_sorted_merge_cond1.bed

./awk_formating.sh 11diff_binding/all_${base_name_alt}_vs_${base_name_1}:${base_name_3}:${base_name_4}_sorted_merge_cond1.bed 11diff_binding/all_${base_name_alt}_vs_${base_name_1}:${base_name_3}:${base_name_4}_peaks_sorted_merge_colswap_cond1.bed

