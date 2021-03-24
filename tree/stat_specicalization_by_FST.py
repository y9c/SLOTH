#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2020 Ye Chang <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2020-02-20 20:25

"""Stat specialization of organs in the tree."""

import sys
import math
from collections import Counter
from collections import defaultdict
from scipy.stats import binom
from ete3 import Tree
import pandas as pd
import numpy as np

# p = binom.pmf(c, s, v)

sample_name: str = sys.argv[1]

tree = Tree(f"./build_tree/{sample_name}.ultrametric.nwk", format=1)

d = defaultdict(lambda: defaultdict(float))

organ_counter = Counter(node.name.split("-")[0] for node in tree)
organ_freq = {k: v / len(tree) for k, v in organ_counter.items()}
organ_ht = {k: 2 * (1 - v) * v for k, v in organ_freq.items()}


dist2root = {tree: 0.0}
for node in tree.iter_descendants(strategy="levelorder"):
    dist2root[node] = dist2root[node.up] + node.dist
    if len(node) > 1:
        d[node.name]["depth_node"] = dist2root[node]
        d[node.name]["depth_up"] = dist2root[node.up]
        node_stat = Counter(l.name.split("-")[0] for l in node)
        for organ, organ_size in node_stat.items():
            # weighted Hs
            Htw = (
                2
                * (organ_size / len(node))
                * (1 - organ_size / len(node))
                * len(node)
                / len(tree)
            )

            d[node.name][organ] = Htw

df = pd.DataFrame(d).T

df_stat = 1- pd.concat(
    [
        df.query("depth_node > @i and depth_up < @i").sum().rename(i)
        for i in np.linspace(0, max(df["depth_node"].max(), df["depth_up"].max()), 100)
    ],
    axis=1,
).drop(index=["depth_node", "depth_up"]).T / pd.Series(organ_ht)

df_stat.to_pickle(f"./stat_specialization/FST_{sample_name}.pkl")
