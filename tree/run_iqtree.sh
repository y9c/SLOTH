#! /usr/bin/env bash
#
# Copyright (C) 2020 Ye Chang <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2020-06-06 16:33
# Discription: run_iqtree

input_fa="./build_tree/A5.bin.fa"
THREAD=48

iqtree2 -nt ${THREAD} -s ${input_fa} -m GTR2+FO+I+R10 -alrt 1000 -bb 1000 -asr --mlrate -wslmr -wspmr
