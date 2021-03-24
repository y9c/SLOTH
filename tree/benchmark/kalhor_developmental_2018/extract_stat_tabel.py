#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2019 yech <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2019-07-25 13:26

"""extract distinct barcode number into table."""

import re

from scipy.special import binom

con_dict = {}
num_dict = {}

with open("./aat9804-Table-S4.txt", "r") as f:
    for l in f.readlines():
        m = re.match(
            "Mutant observed in (\d+) mice out of (\d+) with ([ATGC]+) \(", l
        )
        if m:
            barcode_num, group_id = int(m.group(1)), m.group(3)
            com = binom(barcode_num, 2)
            if group_id in con_dict:
                num_dict[group_id] += barcode_num
                con_dict[group_id] += com
            else:
                num_dict[group_id] = barcode_num
                con_dict[group_id] = com


res = 1
for k, v in con_dict.items():
    res *= num_dict[k] * (num_dict[k] - 1) / v + 1

print(res)
