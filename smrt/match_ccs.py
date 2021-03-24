#!/usr/bin/env python3

"""match sequence with sample barcode and UMI index.

- update in 2018-08-06: support reverse complement
    - method
    one step match
    - pattern
    # 454UPR-(UMI)-3kF-(sequence)-3kR-(sampleBarcode)
"""

import logging
import multiprocessing as mp
import os
import sys

import pysam
import regex

LOGGER: logging.Logger = logging.getLogger()
if not LOGGER.handlers:
    HANDLER: logging.StreamHandler = logging.StreamHandler()
    FORMATTER: logging.Formatter = logging.Formatter(
        "%(asctime)s %(name)-12s %(levelname)-8s %(message)s"
    )
    HANDLER.setFormatter(FORMATTER)
    LOGGER.addHandler(HANDLER)
    LOGGER.setLevel(logging.DEBUG)
    #  LOGGER.setLevel(logging.INFO)


_SEQ_COMPLEMENTS = str.maketrans("ACTG", "TGAC")


def reverse_complement(seq):
    return seq.translate(_SEQ_COMPLEMENTS)[::-1]


def find_dna(p, seq):
    matches = [(m, 1) for m in regex.finditer(p, seq, overlapped=False)] + [
        (m, -1)
        for m in regex.finditer(p, reverse_complement(seq), overlapped=False)
    ]
    if len(matches) == 1:
        return matches[0]
    return (None, None)


def match_read(read):
    seq = read.query_sequence

    # filter 1: check length
    if len(seq) < 2848 or len(seq) > 3448:
        LOGGER.info(f"{read.query_name}:\tlength is not in range")
        return None

    # filter 2: check barcode
    p_barcode = regex.compile(
        r"(?:GGAACTAGTTCCTAAGCTAGTAGGTTAGTA){e<=5}(?:TCTACACAAGGAACAAACACTG){e<=4}((?:GATG){e<=1}(?:[ATGC]{14,18})(?:ACCC){e<=1})"
    )
    barcode_match, barcode_orientation = find_dna(p_barcode, seq)
    if not barcode_match:
        LOGGER.info(f"{read.query_name}:\tbarcode not match")
        return None

    # filter 3: check umi
    p_umi = regex.compile(
        r"(?:TCGCCTCCCTCGCGCCA){e<=3}((?:TCAG){e<=1}(?:[ATGC]{12,16})(?:AGAA){e<=1})(?:GCTTACTAACCAGCCAACTAGC){e<=4}(?:TGGCTAGCAGGTAAACCTGCCAGCCTGCCG){e<=5}"
    )
    umi_match, umi_orientation = find_dna(p_umi, seq)
    if not umi_match:
        LOGGER.info(f"{read.query_name}:\tumi not match")
        return None

    # filter 4: check target
    p_target = regex.compile(
        r"(?:GCTTACTAACCAGCCAACTAGC){e<=4}((?:TGGCTAGCAGGTAAACCTGCCAGCCTGCCG){e<=5}[ATGC]{2780,3180}(?:GGAACTAGTTCCTAAGCTAGTAGGTTAGTA){e<=5})(?:TCTACACAAGGAACAAACACTG){e<=4}((?:GATG){e<=1}(?:[ATGC]{14,18})(?:ACCC){e<=1})"
    )
    target_match, target_orientation = find_dna(p_target, seq)
    if not target_match:
        LOGGER.info(f"{read.query_name}:\ttarget not match")
        return None

    if not (
        umi_orientation == target_orientation == barcode_orientation
        and umi_match.end(1) < target_match.start(1)
        and barcode_match.start(1) > target_match.end(1)
    ):
        LOGGER.debug(f"{read.query_name}:\tmatch in not in correct order")
        return None
    LOGGER.debug(f"{read.query_name}:\tpass match filter")
    orientation_str = "+" if target_orientation == 1 else "-"

    quality = (
        read.query_qualities
        if orientation_str == "+"
        else read.query_qualities[::-1]
    )

    umi_span = f"{umi_match.start(1)}-{umi_match.end(1)}"
    umi_seq = umi_match.group(1)
    umi_qual = quality[umi_match.start(1) : umi_match.end(1)]

    target_span = f"{target_match.start(1)}-{target_match.end(1)}"
    target_seq = target_match.group(1)
    target_qual = quality[target_match.start(1) : target_match.end(1)]

    barcode_span = f"{barcode_match.start(1)}-{barcode_match.end(1)}"
    barcode_seq = barcode_match.group(1)
    barcode_qual = quality[barcode_match.start(1) : barcode_match.end(1)]

    # output match only when fit this step
    read_matched = pysam.AlignedSegment()

    read_matched.query_name = read.query_name + "/matched"
    read_matched.query_sequence = target_seq
    read_matched.query_qualities = target_qual

    read_matched.flag = read.flag
    read_matched.reference_id = read.reference_id
    read_matched.reference_start = read.reference_start
    read_matched.mapping_quality = read.mapping_quality
    read_matched.cigar = read.cigar
    read_matched.next_reference_id = read.next_reference_id
    read_matched.next_reference_start = read.next_reference_start
    read_matched.template_length = read.template_length
    read_matched.set_tags(
        read.get_tags()
        + [
            ("RN", read.query_name),
            ("RN", orientation_str),
            ("TR", target_span),
            ("UR", umi_span),
            ("US", umi_seq),
            ("UQ", umi_qual),
            ("BR", barcode_span),
            ("BS", barcode_seq),
            ("BQ", barcode_qual),
        ]
    )

    return read_matched


def match_read_dict(read_dict):
    """
    {'name': 'm54079_180817_091252/32440632/ccs/rev', '
    flag': '4', '
    ref_name': '*', '
    ref_pos': '0', '
    map_quality': '255', '
    cigar': '*', '
    next_ref_name': '*', '
    next_ref_pos': '0', '
    length': '0', '
    seq': 'GCCTCCCTCGCGCCATCCGTTAGGATT...
    qual': ")+5=/G.8=B5=3AA0,)81%A;D?....
    tags': ['np:i:2', '
    rq:f:0.936408', '
    rs:B:i,6,0,0,0,0,0', '
    sn:B:f,6.47793,12.5605,6.01262,9.98906', '
    za:f:2.75598', '
    zm:i:32440632', '
    zs:B:f,4.01855,4.07533,-1.52781', '
    RG:Z:a751be35']}
    """
    seq = read_dict["seq"]
    name = read_dict["name"]

    # filter 1: check length
    if len(seq) < 2848 or len(seq) > 3448:
        LOGGER.info(f"{name}:\tlength is not in range")
        return None

    # filter 2: check barcode
    p_barcode = regex.compile(
        r"(?:GGAACTAGTTCCTAAGCTAGTAGGTTAGTA){e<=5}(?:TCTACACAAGGAACAAACACTG){e<=4}((?:GATG){e<=1}(?:[ATGC]{14,18})(?:ACCC){e<=1})"
    )
    barcode_match, barcode_orientation = find_dna(p_barcode, seq)
    if not barcode_match:
        LOGGER.info(f"{name}:\tbarcode not match")
        return None

    # filter 3: check umi
    p_umi = regex.compile(
        r"(?:TCGCCTCCCTCGCGCCA){e<=3}((?:TCAG){e<=1}(?:[ATGC]{12,16})(?:AGAA){e<=1})(?:GCTTACTAACCAGCCAACTAGC){e<=4}(?:TGGCTAGCAGGTAAACCTGCCAGCCTGCCG){e<=5}"
    )
    umi_match, umi_orientation = find_dna(p_umi, seq)
    if not umi_match:
        LOGGER.info(f"{name}:\tumi not match")
        return None

    # filter 4: check target
    p_target = regex.compile(
        r"(?:GCTTACTAACCAGCCAACTAGC){e<=4}((?:TGGCTAGCAGGTAAACCTGCCAGCCTGCCG){e<=5}[ATGC]{2780,3180}(?:GGAACTAGTTCCTAAGCTAGTAGGTTAGTA){e<=5})(?:TCTACACAAGGAACAAACACTG){e<=4}((?:GATG){e<=1}(?:[ATGC]{14,18})(?:ACCC){e<=1})"
    )
    target_match, target_orientation = find_dna(p_target, seq)
    if not target_match:
        LOGGER.info(f"{name}:\ttarget not match")
        return None

    if not (
        umi_orientation == target_orientation == barcode_orientation
        and umi_match.end(1) < target_match.start(1)
        and barcode_match.start(1) > target_match.end(1)
    ):
        LOGGER.debug(f"{name}:\tmatch in not in correct order")
        return None
    LOGGER.debug(f"{name}:\tpass match filter")
    orientation_str = "+" if target_orientation == 1 else "-"

    quality = (
        read_dict["qual"]
        if orientation_str == "+"
        else read_dict["qual"][::-1]
    )

    umi_span = f"{umi_match.start(1)}-{umi_match.end(1)}"
    umi_seq = umi_match.group(1)
    umi_qual = quality[umi_match.start(1) : umi_match.end(1)]

    target_span = f"{target_match.start(1)}-{target_match.end(1)}"
    target_seq = target_match.group(1)
    target_qual = quality[target_match.start(1) : target_match.end(1)]

    barcode_span = f"{barcode_match.start(1)}-{barcode_match.end(1)}"
    barcode_seq = barcode_match.group(1)
    barcode_qual = quality[barcode_match.start(1) : barcode_match.end(1)]

    # output match only when fit this step
    read_dict["name"] = name + "/matched"
    read_dict["seq"] = target_seq
    read_dict["qual"] = target_qual
    read_dict["tags"] += [
        f"RN:Z:{name}",
        f"RN:Z:{orientation_str}",
        f"TR:Z:{target_span}",
        f"UR:Z:{umi_span}",
        f"US:Z:{umi_seq}",
        f"UQ:Z:{umi_qual}",
        f"BR:Z:{barcode_span}",
        f"BS:Z:{barcode_seq}",
        f"BQ:Z:{barcode_qual}",
    ]

    return read_dict


def run_iter(input_bam_file, output_bam_file):
    """test run."""
    infile = pysam.AlignmentFile(input_bam_file, "rb", check_sq=False)
    with pysam.AlignmentFile(
        output_bam_file, "wb", template=infile
    ) as outfile:
        for read in infile.fetch(until_eof=True):
            read_matched = match_read(read)
            if read_matched:
                outfile.write(read_matched)


def run_pool(input_bam_file, output_bam_file):
    infile = pysam.AlignmentFile(input_bam_file, "rb", check_sq=False)
    _HEADER = infile.header
    read_list = [read.to_dict() for read in infile.fetch(until_eof=True)]
    LOGGER.info("finish reading file into list")

    with mp.Pool(processes=mp.cpu_count() - 10) as pool:
        matches_list = pool.map(match_read_dict, read_list)
        with pysam.AlignmentFile(
            output_bam_file, "wb", template=infile
        ) as outfile:
            for read_dict in matches_list:
                if read_dict:
                    read_matched = pysam.AlignedSegment.from_dict(
                        read_dict, _HEADER
                    )
                    read_matched.from_dict(read_dict, _HEADER)
                    outfile.write(read_matched)


if __name__ == "__main__":
    #  sample_name = "flyS4_runA_ssccs"
    SAMPLE_NAME = sys.argv[1]
    INPUT_BAM_FILE = f"../pacbio_data/ccs/{SAMPLE_NAME}.bam"
    OUTPUT_BAM_FILE = f"./sequence_matched/{SAMPLE_NAME}_match.bam"

    if not os.path.exists(INPUT_BAM_FILE):
        raise ValueError("BAM FILE of input sample does not exist!")

    # test run
    #  run_iter(INPUT_BAM_FILE, OUTPUT_BAM_FILE)

    # parallel run
    run_pool(INPUT_BAM_FILE, OUTPUT_BAM_FILE)
