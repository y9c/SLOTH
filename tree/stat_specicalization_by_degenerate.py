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

tree.resolve_polytomy(recursive=True)
# annot organ
for node in tree.traverse(strategy="postorder"):
    if node.is_leaf():
        node.organ = set([node.name.split("-")[0]])
    else:
        node.organ = set([o for c in node.get_children() for o in c.organ])

print(len(tree))
# combind organ clade
for node in tree.iter_descendants(strategy="postorder"):
    if len(node.organ) == 1 and len(node.up.organ) > 1:
        node.name = f"{list(node.organ)[0]}-{node.name}"
        for c in node.get_children():
            node.remove_child(c)
print(len(tree))

# calc
for node in tree.traverse(strategy="preorder"):
    if node.is_root():
        node.depth = 0
    else:
        node.depth = node.up.depth + node.dist


path = {}
for leaf in tree.get_leaves():
    node = leaf
    path_dict = {list(node.organ)[0]: node.depth}
    while node.up:
        for delta in node.up.organ - node.organ:
            path_dict[delta] = node.up.depth
        node = node.up
    path[leaf.name] = path_dict

df = pd.DataFrame(path).T
df.to_pickle("./test.pkl")
