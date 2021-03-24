#! /usr/bin/env bash
#
# Copyright (C) 2019 yech <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2018-04-17 11:06
# Discription: cluster UMI

function usage() {
  echo "Usage: $0 -i input_bam -o output_stat" 1>&2
}

while getopts ":i:o:h" opt; do
  case $opt in
  i)
    ┊ input_bam="${OPTARG}"
    ┊
    ;;
  o)
    ┊ output_stat="${OPTARG}"
    ┊
    ;;
  h)
    ┊ usage
    ┊ exit 0
    ┊
    ;;
  *)
    ┊ usage
    ┊ exit 1
    ┊
    ;;
  esac
done

if [ -z "${input_bam}" ] || [ -z "${output_stat}" ]; then
  usage
  exit 1
fi

threads=64
sample=$(basename "${input_bam}" "_match.bam")
cluster_dir=$(dirname "${output_stat}")"/_clustering"
mkdir -p "${cluster_dir}"

# generate umi fastq file
samtools view "${input_bam}" |
  awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^US:Z:|^UQ:Z:"){ td[substr($i,1,2)] = substr($i,6,length($i)-5); } }; print "@"$1"\n"td["US"]"\n+\n"td["UQ"] }' >"${cluster_dir}/${sample}_umi.fq"

# run clustring
usearch -cluster_fast "${cluster_dir}/${sample}_umi.fq" \
  -id 0.95 -gapopen 3.0I/2.0E -gapext 1.0I/0.5E -match +2.0 -mismatch -20.0 -sizeout \
  -uc "${cluster_dir}/${sample}_umi_UsearchClusters.uc" \
  -consout "${cluster_dir}/${sample}_umi_UsearchConsensus.fa" \
  -threads ${threads}

# make mapping table
# ==> clustering/tmp.consout.fasta <==
# >Cluster187422;size=29;
# ----
# H       93963   14      100.0   .       0       14      =       m54079_180801_182145/71500183/ccs       m54079_180801_182145/62391210/ccs
# S       8924    14      *       .       *       *       *       m54079_180801_182145/25821991/ccs       *
# ----

awk 'NR==FNR{m[$1]=$2"\t"$3}NR!=FNR{if($1~/S|H/) print $9"\t"$2"\t"$1"\t"m[$2]"\t"$4}' \
  <(bioawk -c fastx '{split($name,a,";");split(a[1],b,"Cluster");split(a[2],c,"=");print b[2]"\t"c[2]"\t"$seq}' "${cluster_dir}/${sample}_umi_UsearchConsensus.fa") \
  "${cluster_dir}/${sample}_umi_UsearchClusters.uc" >"${output_stat}"
