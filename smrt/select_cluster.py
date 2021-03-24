#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2020 yech <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2020-01-13 16:29

"""select UMIT cluster from fq file."""

from collections import defaultdict
from typing import Dict
import sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def filter_clustered_read(input_file):
    """Filter cluster read in fastq file.

    cluster proportion > x
    cluster read number > x
    read indel proportion < x
    """
    cluster_counter = defaultdict(dict)
    record_dict: Dict[str, SeqRecord] = {}

    for record in SeqIO.parse(input_file, "fastq"):
        cluster_name, cluster_samples = record.description.split()
        record_dict[cluster_name] = record
        umi_id, cluster_id = cluster_name.rsplit("-", 1)
        cluster_size = len(cluster_samples.split(","))
        cluster_counter[umi_id][cluster_id] = cluster_size

    for umi_id, cluster_dict in cluster_counter.items():
        umi_size = sum(cluster_dict.values())
        for cluster_id, cluster_size in cluster_dict.items():
            if cluster_size / umi_size >= 0.75 and cluster_size >= 3:
                # umi_id, cluster_id, cluster_size, umi_size, cluster_size / umi_size,
                res = record_dict[f"{umi_id}-{cluster_id}"]
                if res.seq.count("-") <= 3:
                    print(res.format("fastq"))


if __name__ == "__main__":
    #  "../centroid_seq/A5_ccc.fq"
    INPUT_FILE = sys.argv[1]
    filter_clustered_read(INPUT_FILE)
