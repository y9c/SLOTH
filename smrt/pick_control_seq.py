#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2020 Ye Chang <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2020-05-20 18:35

""""""

import re
import sys

from Bio import SeqIO

# input_file = "./centroid_seq_archived/L5.fq"
input_file = sys.argv[1]

for r in SeqIO.parse(input_file, "fastq"):
    m = re.search("(\w+)-UMI\d+\((\d+)\)-CLUSTER\d+\((\d+)\)", r.id)
    if m is not None:
        org_name = m.group(1)
        umi_size = int(m.group(2))
        cluster_size = int(m.group(3))
        if org_name == "Ck" and umi_size > 3 and cluster_size / umi_size > 0.5:
            print(r.format("fastq"))
