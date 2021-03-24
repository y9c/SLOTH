#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2020 Ye Chang <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2020-02-20 01:35

"""Reset outgroup of contree."""

import statistics
import sys

from Bio import SeqIO
from ete3 import Tree


def convert_to_ultrametric(tree, tree_length=None):
    """Convert tree to ultrametric."""
    # get origin distance to root
    dist2root = {tree: 0.0}
    for node in tree.iter_descendants("levelorder"):
        dist2root[node] = dist2root[node.up] + node.dist

    # get tree length by the maximum
    if not tree_length:
        tree_length = max(dist2root.values())
    else:
        tree_length = float(tree_length)

    dist2root_updated = {tree: 0.0}
    for node in tree.iter_descendants("levelorder"):
        dist_left = tree_length - dist2root_updated[node.up]
        if node.is_leaf():
            node.dist = dist_left
        else:
            median_depth = statistics.mean(dist2root[l] for l in node)
            if median_depth == dist2root[node]:
                node.dist = dist_left
            else:
                node.dist = (
                    dist_left
                    * node.dist
                    / (median_depth - dist2root[node] + node.dist)
                )
        dist2root_updated[node] = dist2root_updated[node.up] + node.dist
    return tree


def modify_tree(tree, mut):
    """Modified tree."""
    tree_out = tree.copy()
    # reroot outgroup
    size = 0
    for n in tree_out.traverse("postorder"):
        if not n.is_leaf():
            if len(n) > size and sum(mut[l.name] for l in n) == 0:
                root = n
    tree_out.set_outgroup(root)

    # remove control sequence
    for l in tree_out:
        if l.name.startswith("Ck-"):
            l.delete(prevent_nondicotomic=True, preserve_branch_length=True)
    # fix root node
    if len(tree_out.get_children()) == 1:
        tree_out = tree_out.get_children()[0]
        tree_out.dist = 0

    ## ultrametric
    # tree_out.convert_to_ultrametric(strategy="balanced")
    # tree_out.convert_to_ultrametric(strategy="fixed")
    tree_out = convert_to_ultrametric(tree_out)

    # rename internal node
    index = 0
    for n in tree_out.traverse("postorder"):
        if not n.is_leaf():
            n.name = f"I{index:04}"
            index += 1

    return tree_out


if __name__ == "__main__":
    sys.setrecursionlimit(10000)

    # sample_name = "A5"
    SAMPLE_NAME = sys.argv[1]

    MUT_DICT = {
        record.id: record.seq.count("1")
        for record in SeqIO.parse(
            f"./build_tree/{SAMPLE_NAME}.fa", format="fasta"
        )
    }

    TREE_IN = Tree(f"./build_tree/{SAMPLE_NAME}.fa.contree")

    TREE_OUT = modify_tree(TREE_IN, MUT_DICT)

    TREE_OUT.write(
        # outfile=f"./build_tree/test_{SAMPLE_NAME}.weighted.nwk",
        outfile=f"./build_tree/{SAMPLE_NAME}.ultrametric.nwk",
        # outfile=f"./build_tree/{SAMPLE_NAME}.nwk",
        format=1,
        format_root_node=False,
    )
