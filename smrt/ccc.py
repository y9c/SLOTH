#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2018 yech <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2018-01-06 17:13


"""Circular Consensus Clustering.

# features

- reference guided (speed enchanced)
- uncompromising qual consolidation (posterior porb)
"""

import math
import multiprocessing as mp
import sys
from collections import defaultdict
from functools import reduce
from itertools import combinations
from typing import Dict, Generator, List, Tuple

import numpy as np
from Bio import SeqIO, pairwise2
from Bio.Alphabet import SingleLetterAlphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scipy.spatial.distance import squareform
from sklearn import cluster


class ASeq:
    """Enhance SeqRecord Object for reference based alignment."""

    def __init__(self, record: SeqRecord, ref: SeqRecord):
        # How to initialize the base class with abc in here?
        self.record: SeqRecord = record
        self.ref = ref
        # init align
        self.align = SeqRecord(Seq("", SingleLetterAlphabet))
        # output
        self.qual = record.letter_annotations["phred_quality"]
        self.mut: List[str] = []
        self.indel: List[str] = []

    def call_mutation(self) -> Tuple[List, List]:
        """Call mutation.

        - one based
        - only report subsitution mutation
        """
        record = self.record
        # ref query match_score mismatch_score gap_open_score gap_extend_score
        alignments = pairwise2.align.globalms(
            str(self.ref.seq), str(record.seq), 2, -3, -2, -1
        )
        #  print(pairwise2.format_alignment(*alignments[0]))
        if len(self.mut) != 0:
            self.mut = []
        if len(self.indel) != 0:
            self.indel = []
        aligned_seq = ""
        aligned_qual = []
        ref_index = 0
        seq_index = 0
        for ref_site, query_site in zip(alignments[0][0], alignments[0][1]):
            if query_site != "-":
                seq_index += 1

            if ref_site != "-":
                ref_index += 1
                aligned_seq += query_site
                # use the prevouse qual as the quality for a gap base
                aligned_qual.append(self.qual[seq_index - 1])
                if (
                    ref_site != query_site
                    and ref_site in "ATGC"
                    and query_site in "ATGC"
                ):
                    self.mut.append(f"{ref_site}{ref_index}{query_site}")

            if ref_site == "-":
                self.indel.append(f"{seq_index}I")

            if query_site == "-":
                self.indel.append(f"{seq_index}D")

        self.align.seq = Seq(aligned_seq, SingleLetterAlphabet)
        self.align.letter_annotations["phred_quality"] = aligned_qual
        return self.mut, self.indel


def get_ref_seq(sample: str) -> str:
    """get reference seqeunce from file."""
    if sample.startswith("flyS2_"):
        return next(SeqIO.parse("../ref/flyS2.fa", "fasta")).seq
    if sample.startswith("flyS4_"):
        return next(SeqIO.parse("../ref/flyS4.fa", "fasta")).seq
    raise FileNotFoundError("no ref file...")


def compare_mutation(a1: ASeq, a2: ASeq) -> float:
    """compare similarity score of two mutation list.

    0: no share
    1: identical
    """

    assert (
        a1.ref.seq == a2.ref.seq
    ), "Sequences are not mapped to the same ref."
    ref_seq = a1.ref.seq
    score = 0
    for ref_base, a1_base, a2_base in zip(ref_seq, a1.align.seq, a2.align.seq):
        if a1_base == a2_base == ref_base:
            score += 1
        elif a1_base == a2_base != ref_base:
            score += 10
        elif a1_base != a2_base:
            if a1_base == ref_base or a2_base == ref_base:
                score -= 3
            else:
                score -= 6
    return score / len(ref_seq)


def shift_record(record: SeqRecord, indel: List) -> SeqRecord:
    """Shift the record position in the coordination of reference seq."""
    ins_record = record[:1]
    ins_record.seq = "-"
    ins_record.letter_annotations["phred_quality"] = [0]

    new_record = record[:0]

    previous_index = 0
    for tag in indel:
        # tag look like '2345D'
        mut_index, mut_type = int(tag[:-1]) - 1, tag[-1]
        mut_index = 0 if mut_index == -1 else mut_index
        if mut_type == "I":
            new_record += record[previous_index:mut_index]
            previous_index = mut_index + 1
        if mut_type == "D":
            new_record += record[previous_index:mut_index] + ins_record
            previous_index = mut_index
    new_record += record[previous_index:]
    return new_record


def generate_cluster(align_list: List[ASeq]) -> Dict[int, List[ASeq]]:
    """Cluster by mutation."""
    aff_matrix = squareform(
        [compare_mutation(a1, a2) for a1, a2 in combinations(align_list, 2)]
    )
    cluster_fit = cluster.AffinityPropagation(
        damping=0.8,
        preference=None,
        max_iter=1000,
        convergence_iter=30,
        affinity="precomputed",
    ).fit(aff_matrix)
    labels = cluster_fit.labels_
    if labels.ndim != 1:
        labels = np.arange(0, len(labels))

    cluster_dict: Dict[int, List[ASeq]] = defaultdict(list)
    for k, align in zip(labels, align_list):
        cluster_dict[k].append(align)
    return cluster_dict


def ungap_record(record: SeqRecord) -> SeqRecord:
    """Drop gap in SeqRecord object."""
    ungapped_record = SeqRecord(
        Seq(str(record.seq).replace("-", ""), SingleLetterAlphabet)
    )
    ungapped_record.letter_annotations["phred_quality"] = [
        record.letter_annotations["phred_quality"][i]
        for i, b in enumerate(record.seq)
        if b != "-"
    ]
    ungapped_record.id = record.id
    ungapped_record.name = record.name
    ungapped_record.description = record.description

    return ungapped_record


def iter_record(record: SeqRecord) -> Generator[Tuple, None, None]:
    """Read nucletide type and error prob base by base from a seq record."""
    for nuc, qual in zip(record, record.letter_annotations["phred_quality"]):
        # prob of error
        prob = 10 ** -(qual / 10)
        yield nuc, prob


def merge_base(
    b1: Tuple[str, float], b2: Tuple[str, float]
) -> Tuple[str, float]:
    """Combind two base with quality stat into a new one.

    b1 = ('A', 0.001)
    b2 = ('T', 0.023)
    """
    n1, p1 = b1
    n2, p2 = b2
    ps = p1 + p2 - 4 * p1 * p2 / 3
    if n1 == n2:
        return n1, (p1 * p2 / 3) / (1 - ps)

    if p1 < p2:
        return n1, p1 * (1 - p2 / 3) / ps
    if p1 > p2:
        return n2, p2 * (1 - p1 / 3) / ps
    # n1 and n2 is equal prob
    return n1, p1 * (1 - p2 / 3) / ps


def merge_identical_bases(*p_list: float) -> float:
    """Combind list of identical base with different quality stat into one.

    b1 = ('A', 0.001)
    b2 = ('A', 0.023)
    """
    # base the order of identical base does not matters
    return reduce(
        lambda p1, p2: (p1 * p2 / 3) / (1 - p1 - p2 + 4 * p1 * p2 / 3), p_list
    )


def merge_different_bases(*p_list: float) -> float:
    """Combind two different base with different quality stat into a new one.

    keep the first one !!!!!! not very corrent here !!!
    """
    return reduce(
        lambda p1, p2: p1 * (1 - p2 / 3) / (p1 + p2 - 4 * p1 * p2 / 3), p_list
    )


def consolidate_base(base_list: List[Tuple[str, float]]) -> Tuple[str, float]:
    """Combind a list of bases with quality stats into a new one."""
    base_dict: Dict[str, Tuple[str, float]] = defaultdict(list)
    for n, p in base_list:
        base_dict[n].append(p)
    base_merged_list = [
        (n, merge_identical_bases(*pl)) for n, pl in base_dict.items()
    ]
    base_type_list, base_prob_list = zip(
        *sorted(base_merged_list, key=lambda x: x[1])
    )
    # the order of identical base does not matters
    return base_type_list[0], merge_different_bases(*base_prob_list)


def consolidate_cluster(align_list: List[ASeq]) -> SeqRecord:
    """Combind a list of reads with quality stats into a new one."""
    pmin = 10 ** (-93 / 10)
    assert all(
        len(r.align) == len(align_list[0].align) for r in align_list
    ), "The Sequence are not the same in length!!!"

    seq = ""
    qual = []
    for b in zip(*[iter_record(r.align) for r in align_list]):
        s, p = consolidate_base(b)
        if p < pmin:
            p = pmin
        q = int(-10 * math.log10(p))
        seq += s
        qual.append(q)
    cons_record = SeqRecord(Seq(seq, SingleLetterAlphabet))
    cons_record.letter_annotations["phred_quality"] = qual
    #  cons_record = ungap_record(cons_record)
    return cons_record


def run_call(a):
    """Run mutation calling."""
    a.call_mutation()
    print(a.record.id.split("/")[1:], a.mut)
    return a


def read_group(ref_file: str, *query_files: str) -> Dict[str, List[ASeq]]:
    """Read UMI group.

    - query file is in fastq format
    - ref file is in fasta format
    """
    min_size = 2
    max_size = 100
    ref = next(SeqIO.parse(ref_file, "fasta"))
    record_group: Dict[str, List] = defaultdict(list)
    for query in query_files:
        for record in SeqIO.parse(query, "fastq"):
            # name, organ umi_seq, umi_qual
            align = ASeq(record, ref)
            #  a.call_mutation()
            #  group_name = "-".join(align.record.description.split("\t")[1:3])
            group_name = align.record.description.split("\t")[1]
            record_group[group_name].append(align)

    name_queue: List[str] = []
    align_queue: List[SeqRecord] = []
    # run filter
    for name, align_list in record_group.items():
        if min_size <= len(align_list) <= max_size:
            name_queue += [name] * len(align_list)
            align_queue += align_list
    # run call in parallel
    with mp.Pool(mp.cpu_count() // 2 - 2) as p:
        align_queue_called = p.map(run_call, align_queue)

    record_group_filtered: Dict[str, List] = defaultdict(list)
    for k, v in zip(name_queue, align_queue_called):
        record_group_filtered[k].append(v)

    return record_group_filtered


if __name__ == "__main__":
    #  testing
    SAMPLE = sys.argv[1]
    SEQ_GROUP = read_group(
        "./ref/hmf32_bioseq.fa",
        #  "../test/A7C1_CGGAGGGACGCCGGCGGG.fq"
        f"./annnotated_umi/{SAMPLE}_umi_annotated.fq",
    )

    OUTFILE = open(f"./centroid_seq/{SAMPLE}_centroid.fq", "w")
    for group_name, group_record_list in SEQ_GROUP.items():
        for cluster_name, cluter_record_list in generate_cluster(
            group_record_list
        ).items():
            cons_record = consolidate_cluster(cluter_record_list)
            cons_record.id = f"{group_name}-{cluster_name:02}"
            cons_record.name = f"{group_name}-{cluster_name:02}"
            cons_record.description = (
                ",".join([record.record.id for record in cluter_record_list])
                + f"\tGROUP_SIZE={len(group_record_list)}"
                + f"\tCLUSTER_SIZE={len(cluter_record_list)}"
            )
            OUTFILE.write(cons_record.format("fastq"))
    OUTFILE.close()
