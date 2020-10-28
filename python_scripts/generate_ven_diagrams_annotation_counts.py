import matplotlib
matplotlib.use('Agg')
import pybedtools
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import argparse
import pybedtools
import argparse
import os
import itertools
import copy


def read_in_bed_file(bed_file):
    """
    Reads in a bed file using the pybedtools standard appraoch.
    """
    def intersecting_feature(feature, index):
        """TODO: Docstring for intersecting_feature(feature, L.
        :returns: TODO
        """
        return(feature[index]) != '.'

    try:
        read_file = pybedtools.BedTool(bed_file)
    except: 
        print("File DOES NOT EXISTS")
        exit(1)

    return read_file 


def shared_peaks(bed_file_1, bed_file_2):
    """
    :returns: TODO

    """
    intersection = bed_file_1 + bed_file_2
    middle_number = intersection.count()
    return middle_number


def unique_peaks(bed_file_1, bed_file_2):
    """TODO: Docstring for unique_peaks.

    :arg1: TODO
    :returns: TODO

    """
    bed_1_only = bed_file_1 - bed_file_2
    bed_2_only = bed_file_2 - bed_file_1

    count_1 = bed_1_only.count()
    count_2 = bed_2_only.count()
    
    return (count_1, count_2)



def get_parser():
    parser = argparse.ArgumentParser(description='Finds peaks shared between \
            replicate peak calls, as well as unqiue peaks to each replicate and \
            outputs said peaks. ')
    parser.add_argument('-bed1','--bed_file1', help='bed1 file',\
            nargs="+", required=True, dest='bed')

    parser.add_argument('-header_name','--headers', help='bed3 file', \
          nargs="+", required=False, dest='heads', type=str)

    parser.add_argument('-title','--title', help='bed3 file', \
            required=False, dest='title', type=str)

    parser.add_argument('-o','--output_name', help='output', \
            required=False, dest='o')

    args = vars(parser.parse_args())
    return parser

#python generate_ven_diagrams.py -bed1 00.data/histone_mods/tis_root_mod_H3K36me3_group_broad_merged.bed 00.data/histone_mods/tis_root_mod_H3K4me1_group_broad_merged.bed -header_name H3K36me3 H3K4me1 > root_extension.txt



if __name__ == "__main__":
    args = get_parser().parse_args()
    
    #Read in tissue names
    tissue_list = args.heads

    #Read Files
    bed_list = [read_in_bed_file(i) for i in list(args.bed)]

    bed_1 = bed_list[0]
    bed_2 = bed_list[1]
    bed_3 = bed_list[2]

    header_1_2 = args.heads[0] + "_" + args.heads[1]
    header_1_3 = args.heads[0] + "_" + args.heads[2]
    header_2_3 = args.heads[1] + "_" + args.heads[2]

    intersection_all = str((bed_1 + bed_2 + bed_3).count())
    intersection_1_2  = str((bed_1 + bed_2).count())
    intersection_1_3 = str((bed_1 + bed_3).count())
    intersection_2_3 = str((bed_2 + bed_3).count())
    b1_only = str((bed_1).count())
    b2_only = str((bed_2).count())
    b3_only = str((bed_3).count())
    
    print('\t'.join(["all", header_1_2, header_1_3, header_2_3, args.heads[0], args.heads[1], args.heads[2]]))
    print('\t'.join([intersection_all, intersection_1_2, intersection_1_3, intersection_2_3, b1_only, b2_only, b3_only]))

