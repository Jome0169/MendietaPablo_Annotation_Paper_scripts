"""
Loads in BED files, and generates the overlap between a given set of THREE Bed
files

"""
from datetime import datetime
import argparse
import sys
import os
import pybedtools


def read_bed_file(bed_file):
    """Takes in bed file and reads it in using pybedtools. Will make life
    easier longitudinally.

    :bed_file: TODO
    :returns: TODO

    """
    print("Loading Bed file %s" % bed_file)

    try:
        os.path.isfile(bed_file)
    except:
        FileNotFoundError
    bed_file_load = pybedtools.BedTool(bed_file).sort()
    return bed_file_load


def get_parser():
    parser = argparse.ArgumentParser(description='Reads in 3 bed files and generates values for a ven diagram')
    parser.add_argument('-beds','--bed_files', help='Complted pasa genes from complete_gene_retriever.py ', \
            nargs="+", required=True, dest='beds', type=str)
    parser.add_argument('-names','--file_names', help='Names to be added into the Ven diagram', \
        nargs="+", required=True, dest='names', type=str)
    parser.add_argument('-script_name','--script_name', help='Names to be added into the Ven diagram', \
        required=True, dest='script', type=str)
    parser.add_argument('-dir','--dir_name', help='Names to be added into the Ven diagram', \
        required=True, dest='dir', type=str)
    parser.add_argument('-o','--output', help='Base output file name',\
            required=False,dest='o') 
    
    args = vars(parser.parse_args())
    return parser

if __name__ == "__main__":
    args = get_parser().parse_args()
    bed_files = args.beds

    read_in_bed_files = [read_bed_file(x) for x in bed_files]

    output_script_name = os.path.join(args.dir, args.script)
    
    if (args.o).endswith(".png"):
        output_file_name_final = os.path.join(args.dir, args.o)

    elif (args.o).endswith(".png") != True:
        output_file_name = args.o + ".png"
        output_file_name_final = os.path.join(args.dir, output_file_name)

    if args.names != None:
        pybedtools.contrib.venn_maker.venn_maker(read_in_bed_files, names =args.names,
            figure_filename=output_file_name_final,
            script_filename=output_script_name, run=True)
