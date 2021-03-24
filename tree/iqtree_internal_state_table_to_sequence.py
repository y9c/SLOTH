#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2020 Ye Chang <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2020-06-05 01:24

"""Convert iqtree .state ouput (tsv file) into sequence (fasta file)."""


import sys

import pandas as pd

# seq_file = "./annnotated_tree/A5.internal.fa"
# state_file = "./build_tree/A5.bin.fa.state"

SAMPLE = sys.argv[1]
seq_file = f"./annnotated_tree/{SAMPLE}.internal.fa"
state_file = f"./build_tree/{SAMPLE}.bin.fa.state"

df = pd.read_csv(state_file, sep="\t", comment="#", low_memory=False)


def f(row):
    if row["State"] == "-":
        if row["p_0"] > row["p_1"]:
            return "0"
        return "1"
    return str(row["State"])


# df['State'] = df.apply(f, axis=1)

df.loc[(df["State"] == "-") & (df["p_0"] > df["p_0"]), "State"] = "0"
df.loc[(df["State"] == "-") & (df["p_0"] <= df["p_0"]), "State"] = "1"


with open(seq_file, "w") as f:
    for n, m in df.groupby("Node")["State"].apply("".join).iteritems():
        f.write(f">{n}\n{m}\n")
