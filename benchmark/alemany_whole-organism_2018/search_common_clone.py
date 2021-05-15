#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2019 yech <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2019-07-25 22:18

"""find common mutation profile clone."""


import pandas as pd

df1 = (
    pd.read_csv("./R1_scarclones.txt", sep="\t")
    .drop(columns=["oclust", "hclust"])
    .melt(id_vars="CellID", var_name="Site", value_name="Freq")
    .mask(lambda x: x["Freq"] == 0)
    .dropna()
    .groupby("CellID")["Site"]
    #  .apply(len)
    .apply(tuple)
    .drop_duplicates()
)


df2 = (
    pd.read_csv("./R3.txt", sep="\t")
    .drop(columns=["oclust", "hclust", "organ"])
    .rename(columns={"Unnamed: 0": "CellID"})
    .melt(id_vars="CellID", var_name="Site", value_name="Freq")
    .mask(lambda x: x["Freq"] == 0)
    .dropna()
    .groupby("CellID")["Site"]
    #  .apply(len)
    .apply(tuple)
    .drop_duplicates()
)

#  print(df1.mean())

print(df1.shape)
print(df2.shape)
print(pd.concat([df1, df2]).value_counts().value_counts())
