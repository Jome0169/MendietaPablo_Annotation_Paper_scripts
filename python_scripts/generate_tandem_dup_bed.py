import argparse
import pybedtools
import argparse
import os
import itertools
import copy






def generate_dict_bed(tandem_dup_csv):
    """Takes in tandem duplicates CSV and puts tandem dup name as value, gene
    name as key

    :tandem_dup_csv: TODO
    :returns: TODO
    """
    gene_names = {}

    with open(tandem_dup_csv, 'r') as f:
        for line in f:
            clean_line = line.replace("\"","").strip().split(",")
            for name in clean_line[1:]:
                final_name = "gene:" + name
                gene_names[final_name] = clean_line[0]
    return gene_names


def read_bed_file(bed_file):
    """Reads in bed file putting gene name into the dict

    """
    dict_bed = {}
    with open(bed_file,'r') as f:
        for line in f:
            clean_line = line.strip().split('\t')
            dict_bed[clean_line[3]] = clean_line
    return dict_bed



def generate_tandem_dup_bed(tandem_dup_dict, bed_dict):
    """TODO: Docstring for generate_tandem_dup_bed.
    :returns: TODO

    """
    final_list = []

    for tandem_name, tandem_group in tandem_dup_dict.items():
        try:
            pull_bed_region = bed_dict[tandem_name]
            updated_name = pull_bed_region[3] + "_" + tandem_group
            pull_bed_region[3] = updated_name
            final_list.append(pull_bed_region)
        except KeyError:
            pass
    return final_list




def get_parser():
    parser = argparse.ArgumentParser(description='Finds peaks shared between \
            replicate peak calls, as well as unqiue peaks to each replicate and \
            outputs said peaks. ')
    parser.add_argument('-tandem_dup','--tandem', help='bed1 file',\
            required=True, dest='tandem')
    parser.add_argument('-bed','--bed_file', help='bed3 file', \
          required=False, dest='bed', type=str)
    args = vars(parser.parse_args())
    return parser



if __name__ == "__main__":
    args = get_parser().parse_args()

    gene_name_dict = generate_dict_bed(args.tandem)
    gene_bed_file = read_bed_file(args.bed) 

    final_list = generate_tandem_dup_bed(gene_name_dict, gene_bed_file)
    
    for name in final_list:
        print('\t'.join(name))
