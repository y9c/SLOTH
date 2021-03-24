#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2020 Ye Chang <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2020-06-04 11:35

"""tree dist vs mutation number.
"""

from ete3 import Tree
from Bio import SeqIO
from collections import Counter

t = Tree(
    "/data/resource/project/FlyHMF3k_YeChang2020/tree_archived/L5.bin.fa.contree"
)
for n in t.traverse("preorder"):
    if n.is_root():
        n.depth = n.dist
    else:
        n.depth = n.up.depth + n.dist


d = {
    s.name.replace("(", "_").replace(")", "_"): Counter(s.seq)["1"]
    for s in SeqIO.parse(
        "../backup/build_tree_with_append_legecy/L5.bin.fa",
        "fasta",
    )
}

"""
d2 = {
    s.name.replace("(", "_").replace(")", "_"): Counter(s.seq)["1"]
    for s in SeqIO.parse(
        "/data/resource/project/FlyHMF3k_YeChang2020/centroid_seq_archived/L5_core.fq",
        "fastq",
    )
}

d = {**d1, **d2}
"""

for n in t:
    if not n.name.startswith("root"):
        print(n.name, n.depth, d[n.name])

