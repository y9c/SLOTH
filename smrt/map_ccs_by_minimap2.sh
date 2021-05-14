#! /usr/bin/env bash
#
# Copyright (C) 2019 yech <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2019-10-31 02:48
# Discription: run_mapping

# generate the index files first
# minimap2 -d ref.mmi ref.fa

input_bam=$1
output_bam=$2
threads=24
ref="./hmf32.mmi"

samtools fastq ${input_bam} |
  minimap2 -t ${threads} -O 2,12 -E 2,1 --score-N 0 --end-bonus 100 -a --MD -x map-pb ${ref} - |
  samtools sort -@ ${threads} >${output_bam}

samtools index -@ ${threads} ${output_bam}
