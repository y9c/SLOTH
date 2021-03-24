#! /usr/bin/env bash
#
# Copyright (C) 2019 yech <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2019-10-31 13:18
# Discription: count_base

input_bam=$1
output_tsv=$2

pysamstats -D 0 --fasta ./ref/hmf32.fa --type variation ${input_bam} |
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$10,$12,$14,$16,$18,$20,$22}' >${output_tsv}
