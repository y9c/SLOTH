#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright © 2019 yech <yech1990@gmail.com>
#
# Distributed under terms of the MIT license.

"""Match sequence with sample barcode and UMI index.

- update in 2019-11-03: use mapped result
"""

import multiprocessing as mp
import os
import sys
from dataclasses import dataclass
from typing import Tuple

import edlib
import pysam

from utils import get_logger, seq_revcom

LOGGER = get_logger(__name__)

#  SAMPLE_NAME = "A5"
SAMPLE_NAME = sys.argv[3]

BARCODE_DICT = {}
with open(
    "ref/sample_barcodes.tsv",
    "r",
) as f_barcode:
    for line in f_barcode.readlines():
        ind, org, _, seq = line.replace("\n", "").split("\t")
        if ind == SAMPLE_NAME:
            BARCODE_DICT[org] = seq


if SAMPLE_NAME == "A5":
    _LIB5_SEQ = "TACATGATAGACGACA"  # A5
elif SAMPLE_NAME == "A7":
    _LIB5_SEQ = "ATGTAGTAGTGAGCAT"  # A7
elif SAMPLE_NAME == "L5":
    _LIB5_SEQ = ""  # L5
elif SAMPLE_NAME == "L6":
    _LIB5_SEQ = "AGTATCATGTGTATCT"  # L6
_LIB3_SEQ = seq_revcom(_LIB5_SEQ)


@dataclass
class Seq:
    seq: str = ""
    qual: str = ""


@dataclass
class MatchedSeq:
    lib5: Seq
    sam5: Seq
    umi: Seq
    bio: Seq
    sam3: Seq
    lib3: Seq


def check_library_barcode(b1: str, b2: str) -> bool:
    """only some library with libary barcode."""
    #  if len(b1) == 0 or len(b2) == 0:
    #  return False
    d1 = edlib.align(b1, _LIB5_SEQ, mode="NW", task="distance")["editDistance"]
    d2 = edlib.align(b2, _LIB3_SEQ, mode="NW", task="distance")["editDistance"]
    return d1 <= 3 and d2 <= 3


def edit_barcode(barcode_seq: str) -> Tuple[str, int]:
    """map the nearest sample barcode."""
    # correct barcode
    min_distance = len(barcode_seq)
    organ_matched = ""
    if min_distance == 0:
        return "", 0
    for organ, barcode in BARCODE_DICT.items():
        result = edlib.align(barcode_seq, barcode, mode="NW", task="distance")
        dis = result["editDistance"]
        if dis <= min_distance:
            min_distance = dis
            organ_matched = organ

    return organ_matched, min_distance


def check_sample_barcode(b1, b2, dis_cutoff=0.1):
    """check pair end barcode.

    dis_cutoff: proportion of mistach
    """
    # too strick filter
    if len(b1) == 0 or len(b2) == 0:
        return None
    org1, dis1 = edit_barcode(b1)
    org2, dis2 = edit_barcode(b2)
    if (
        org1 == org2 != ""
        and dis1 / len(b1) < dis_cutoff
        and dis2 / len(b2) < dis_cutoff
    ):
        return org1
    return None


def lookup_pos(pos, span_dict=None):
    """
    [0:16] : library forward barcode =>
    [21:37] : sample barcode =>
    [57:71] : UMI
    [79:3019] : HMF =>
    [95:3038] : BioSeq =>
    [3062:3078] : sample barcode =>
    [3083:3099] : library barcode <=
    """
    if span_dict is None:
        if _LIB5_SEQ or 1 == 1:
            span_dict = {
                range(0, 16): "lib5",
                range(21, 37): "sam5",
                range(57 - 2, 71 + 2): "umi",
                range(95, 3038): "bio",
                range(3062, 3078): "sam3",
                range(3083, 3099): "lib3",
            }
        else:
            span_dict = {
                range(5, 21): "sam5",
                range(41, 55): "umi",
                range(79, 3022): "bio",
                range(3046, 3062): "sam3",
            }

    for span in span_dict:
        if pos in span:
            return span_dict[span]


def parse_read(read):
    """parse bioseq to fastq."""
    #  dele = "-"
    dele = ""

    # check if there is MD tag
    # check whether query sequence is exist
    if not read.has_tag("MD") or (query_seq := read.query_sequence) is None:
        return None

    query_qual = read.query_qualities

    # init matched result object
    matched = MatchedSeq(
        lib5=Seq(seq="", qual=""),
        sam5=Seq(seq="", qual=""),
        umi=Seq(seq="", qual=""),
        bio=Seq(seq="", qual=""),
        sam3=Seq(seq="", qual=""),
        lib3=Seq(seq="", qual=""),
    )

    ref_pos_pointer = 0
    for query_pos, ref_pos, ref_base in read.get_aligned_pairs(with_seq=True):
        if ref_pos is not None:
            ref_pos_pointer = ref_pos
        feature_type = lookup_pos(ref_pos_pointer)
        if feature_type:
            getattr(matched, feature_type).seq += (
                query_seq[query_pos] if query_pos is not None else dele
            )
        if feature_type in ["bio", "umi"]:
            getattr(matched, feature_type).qual += (
                chr(query_qual[query_pos] + 33)
                if query_pos is not None
                else dele
            )

    #  print(
    #  len(matched.bio.seq) > 2000,
    #  check_library_barcode(matched.lib5.seq, matched.lib3.seq),
    #  (organ := check_sample_barcode(matched.sam5.seq, matched.sam3.seq)) is not None
    #  )
    if (
        len(matched.bio.seq) > 2000
        and check_library_barcode(matched.lib5.seq, matched.lib3.seq)
        and (organ := check_sample_barcode(matched.sam5.seq, matched.sam3.seq))
        is not None
    ):
        return (
            organ,
            matched.umi.seq,
            matched.umi.qual,
            matched.bio.seq,
            matched.bio.qual,
        )
    return None


def run_iter(input_bam_file, output_fq_file):
    """test run."""
    with pysam.AlignmentFile(
        input_bam_file, "rb", check_sq=False
    ) as infile, open(output_fq_file, "w") as outfile:
        for read in infile.fetch(until_eof=True):
            result = parse_read(read)
            if result:
                outfile.write(
                    f"@{read.qname}\t{result[0]}\t{result[1]}\t{result[2]}\n{result[3]}\n+\n{result[4]}\n"
                )


if __name__ == "__main__":
    #  SAMPLE_NAME = "flyS4_runA_ssccs"
    #  INPUT_BAM_FILE = f"../pacbio_ccs/{SAMPLE_NAME}.bam"
    #  OUTPUT_BAM_FILE = f"./matched_sequence/{SAMPLE_NAME}_match.bam"
    INPUT_BAM_FILE = sys.argv[1]
    OUTPUT_BAM_FILE = sys.argv[2]
    if not os.path.exists(INPUT_BAM_FILE):
        raise ValueError("BAM FILE of input sample does not exist!")
    # run the matching
    run_iter(INPUT_BAM_FILE, OUTPUT_BAM_FILE)
