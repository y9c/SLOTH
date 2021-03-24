#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2019 yech <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2019-07-25 16:05

"""count and average mutation of barcodes.

taxa    sample  count   proportion      event
N0      7B_Brain        5441    0.16485880499333413     12D+136_46D+157_46D+157_98D+210_98D+210_98D+210_98D+210_27D+331_27D+331_1D+382
N1      7B_Brain        4802    0.1454975154526724      7D+135_9D+160_NONE_5D+220_44D+230_44D+230_21D+281_76D+306_76D+306_76D+306
"""

import re

import pandas as pd

df = pd.read_csv(
    #  "./GSE81713_fish_ADR1_PHYLIP_MIX_gte5_input.annotations.txt", sep="\t"
    "./GSE81713_fish_ADR2_PHYLIP_MIX_gte5_input.annotations.txt",
    sep="\t",
)


def count_event(l):
    event = []
    for s in l:
        m = re.search("\+(\d+)", s)
        if m:
            event.append(m.group(1))
    return len(set(event))


mut_sum = (
    df["event"].str.split("_").apply(lambda x: count_event(x)) * df["count"]
).sum()

read_sum = df["count"].sum()

print(mut_sum / read_sum)
