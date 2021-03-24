#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2020 yech <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2020-01-14 14:37

"""call mutaion from centroid seq."""


import multiprocessing as mp
import sys

from Bio import SeqIO

from ccc import ASeq

#  testing
SAMPLE = sys.argv[1]
#  SAMPLE = "L5"


def run_call(a):
    a.call_mutation()
    return a


ref = next(SeqIO.parse("../ref/hmf32_bioseq.fa", "fasta"))
align_list = [
    ASeq(record, ref)
    for record in SeqIO.parse(f"../selected_centroid/{SAMPLE}_ccc.fq", "fastq")
]

with mp.Pool(mp.cpu_count() // 2) as p:
    align_queue_called = p.map(run_call, align_list)

with open(f"../called_mutations/{SAMPLE}_mut.tsv", "w") as of:
    for align in align_queue_called:
        of.write(f"{align.record.id}\t{align.mut}\t{align.indel}\n")
