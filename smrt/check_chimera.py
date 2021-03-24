#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2019 yech <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2019-11-06 00:51

"""check chimera."""


def get_diff(seq1, seq2):
    """
    @return dict of position (0-based) and bases
    """
    assert len(seq1) == len(seq2), "The length of sequence is not equal"
    return {
        i: (b1, b2) for i, (b1, b2) in enumerate(zip(seq1, seq2)) if b1 != b2
    }


def get_info(site, info):
    """
    site: (ref_pos, query_base)
    """
    ref_pos, query_base = site
    if query_base in info[ref_pos]:
        return info[ref_pos].index(query_base)
    return None


if __name__ == "__main__":
    SEQ1 = "ATCG"
    SEQ2 = "AAAA"
    diff_info = get_diff(SEQ1, SEQ2)
    # certain site
    SITE = (2, "A")
    tmp = get_info(SITE, diff_info)
    print(tmp)
