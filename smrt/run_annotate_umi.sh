#! /usr/bin/env bash
#
# Copyright (C) 2020 yech <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2020-01-11 20:32
# Discription: run_cluster_umi

# extract UMI sequence into fq file
declare -a Samples=("L5" "L6" "A5" "A7")

# Iterate the string array using for loop
for sample in ${Samples[@]}; do
  echo $sample
  awk 'NR==FNR{a[$5]=$1}NR!=FNR{if($1 in a){printf "@%s\t%s-%04i\n%s\n+\n%s\n", $1, $4, a[$1], $2, $3}}' \
    <(cat "./clustered_umi/${sample}_"*.tsv) \
    <(bioawk -c fastx '{print $name,$seq,$qual,$4}' "./matched_barcode/${sample}C1_ssccs_matched.fq" "./matched_barcode/${sample}C2_ssccs_matched.fq") >"./annnotated_umi/${sample}_umi_annotated.fq"
done
