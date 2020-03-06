# -*- coding: utf-8 -*-
"""
    48.single_bp_split_beds.single_bp_split_bed
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Purpospe of this script is to take in a bed file and generate a seperate
    bed file of intervals at sinlge NT resultion 500 BP UPSTREAM of the gene.
    
    This has been developed to bed used for R2C2 which grants us single bp
    resolution of TSS.

    :copyright: (c) 2020 by YOUR_NAME
    :license: LICENSE_NAME, see LICENSE for more details.
"""





from sys import argv
import copy

import argparse
import sys
import os
import itertools 


def remove_if_exists(filename):
    """TODO: Docstring for remove_if_exists.
    :returns: TODO
    """
    try:
        os.remove(filename)
    except OSError:
        pass


def read_in_bed(file_name):
    """TODO: Docstring for read_in_bed.
    :returns: TODO

    """
    line_list = []
    with open(file_name, 'r') as f:
        for line in f:
            clean_line = line.strip().split()
            line_list.append(clean_line)
    return line_list


def generate_range(nested_bed):
    """TODO: Docstring for generate_single_bp_regions(.

    """
    single_bp_beds = []
    strand = nested_bed[5]
    chrom = nested_bed[0]
    start = nested_bed[1]
    stop = nested_bed[1]
    transcript_name = nested_bed[3]

    counter = 1
    if strand == '+':
        if int(start) - 1000 < 0:
            for bp_val in range(0, int(start), 1):
                generate_name = transcript_name + "__" + "bp" + "__" + str(counter)
                micro_bed = [chrom, str(bp_val), str(bp_val + 1), generate_name,
                        ".", strand]
                single_bp_beds.append(micro_bed)
                counter += 1
        elif int(start) - 1000 > 0:
            for bp_val in range(int(start) - 1000, int(start), 1):
                generate_name = transcript_name + "__" + "bp" + "__" + str(counter)
                micro_bed = [chrom, str(bp_val), str(bp_val + 1), generate_name,
                        ".", strand]
                single_bp_beds.append(micro_bed)
                counter += 1

    elif strand == '-':
        for bp_val in range(int(stop), int(stop) + 1000, 1):
            generate_name = transcript_name + "__" + "bp" + "__" + str(counter)
            micro_bed = [chrom, str(bp_val), str(bp_val + 1), generate_name,
                    ".", strand]
            single_bp_beds.append(micro_bed)
            counter += 1

    return single_bp_beds







bed_file = read_in_bed(argv[1])
remove_if_exists(argv[2])

with open(argv[2], 'a+') as f:
    for item in bed_file:
        single_bp_res = generate_range(item)
        for single_bp_bed in single_bp_res:
            f.write('\t'.join(single_bp_bed))
            f.write('\n')
