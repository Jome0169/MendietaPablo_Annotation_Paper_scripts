# -*- coding: utf-8 -*-
"""

    Mon May 20 13:54:57 EDT 2019
    gtf_test.check_for_multiple_assemblies
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Purpose of this script is to read in a gtf file, and look for possible
    multiple assemblies that need to be merged. This could have a large effect in that a region could
    contain multiple assembled transctips

    :copyright: (c) 2019 by Pablo Mendieta SON
    :license: LICENSE_NAME, see LICENSE for more details.
"""

from random import choice
from string import ascii_uppercase
import itertools
import copy
import argparse
import sys
import os
import gffutils
import pybedtools
from gffutils import helpers
from gffutils.pybedtools_integration import to_bedtool, featurefuncs


def remove_file(file_name):
    """removes file if it does exist. Don't want to mess with appending garbage
    to wrong file

    :file_name: TODO
    :returns: TODO

    """
    try:
        os.remove(file_name)
    except OSError:
        pass

def check_empty_file(gtf_file):
    """Checks if the file is empty. Returns 1 if fine, returns 0 if not fine

    :gtf_file: TODO
    :returns: TODO

    """
    try:
        gtf_file_load = gffutils.create_db(gtf_file, database_random_name, \
            merge_strategy="merge", infer_gene_extent=True, force=True, \
            disable_infer_transcripts=True, id_spec="gene_id")
        return 1 
    except:
        return 0 

def generate_directory(output_name):
    """ Been having tons of issues if the output directory isn't made. So
    attempting to bypass thi be foolproofing directory generation and ensuring
    non-zero exit
    """

    split_output_name = output_name.split('/')

    if len(split_output_name) > 1:
        #Don't take file name
        split_output_name_remove_file = split_output_name[:-1]
        re_join_nest = '/'.join(split_output_name_remove_file)
        try:
            os.makedirs(re_join_nest)
        except FileExistsError:
        # directory already exists
            pass

    #No directory structure, just pass
    else:
        pass



def read_gtf_file(gtf_file):
    """reads in the gtf file using gffutils

    :gtf_file: TODO
    :returns: TODO

    """
    gtf_file_load = gffutils.create_db(gtf_file, database_random_name, \
            merge_strategy="merge", infer_gene_extent=True, force=True, \
            disable_infer_transcripts=True, id_spec="gene_id")
    gtf_file_load = gffutils.FeatureDB(database_random_name)

    return(gtf_file_load)


def check_multiple_genes(gtf_file_load):
    """checks how many "gene" features there are. If there is exactly one,
    return and write merged isoforms to file

    If greater than one we'll return a value that will cause downstream overlap
    checks

    :read_gtf_file: TODO
    :returns: TODO

    """
    gene_count = gtf_file_load.count_features_of_type("transcript")

    if gene_count == 1:
        return 0
    elif gene_count > 1:
        return 1

def get_feature(type_feature, gtf_file_load):
    """brief function to retrive all features of a specific type

    :gtf_file_load: TODO
    :returns: TODO

    """
    feature_group = gtf_file_load.features_of_type(type_feature)
    return(feature_group)


def check_overlaps(gtf_file_load_types):
    """Reads in the transcripts given, converts them to a bedtools interval
    format, and checks for overlap. If no overlap is foudn that's fine.

    If overlap found, will return a list of items to merge.

    :gtf_file_load_types: TODO
    :returns: TODO

    """
    number_total_tran = gtf_file_load_types.count_features_of_type('transcript')
    #Retrieve for merge
    all_transcripts = get_feature("transcript", loaded_gtf_file)
    #merge, see if different
    merged_featuers = gtf_file_load_types.merge(all_transcripts, ignore_strand=True)
    number_merged_tran = sum(1 for x in merged_featuers)

    if number_total_tran == number_merged_tran:
        return 0

    elif number_total_tran != number_merged_tran:
       return 1

def check_merging_pairs(transcript1, transcript2, gtf_file_load_types):
    """TODO: Docstring for check_merging_pairs.

    :transcript1: TODO
    :transcript2: TODO
    :returns: TODO

    """
    merge_attempt = gtf_file_load_types.merge([transcript1, transcript2], ignore_strand = True)
    number_merged_tran = sum(1 for x in merge_attempt)

    if 2 == number_merged_tran:
        return 0
    elif 2 != number_merged_tran:
       return 1

def return_copy_to_delete(transcript1, transcript2, gtf_file_load_types):
    """Compare coverage of hte two overlapping isoforms. return the name of the
    transcript with the LOWEST isoform, and delete that sucker.
    :returns: transcript to be deleted.

    """
    tran1_start = int(transcript1.start)
    tran2_start = int(transcript2.start)

    tran1_stop = int(transcript1.stop)
    tran2_stop = int(transcript2.stop)

    tran1_total = tran1_stop - tran1_start
    tran2_total = tran2_stop - tran2_start


    if tran1_total > tran2_total:
        return transcript2
    elif tran1_total < tran2_total:
        return transcript1
    elif tran1_total == tran2_total:
        return transcript2
    else:
        pass


def merge_exons(gtf_file_load):
    """TODO: Docstring for merge_exons.

    :gtf_file_load): TODO
    :returns: TODO

    """

    all_exons = get_feature("exon", gtf_file_load)
    all_transcripts= get_feature("transcript", gtf_file_load)
    #merge, see if different
    merged_transcripts = gtf_file_load.merge(all_transcripts, ignore_strand=False)
    merged_exons = gtf_file_load.merge(all_exons, ignore_strand=False)
    
    final_list = []
    for i in merged_transcripts:
        i.source = "Stringtie"
        i.score = "1000"
        i.frame = "."
        i["transcript_id"] = "MSTRG.1.1"
        i["gene_id"] = "MSTR.1"
        final_list.append(i)
    
    counter = 1
    for x in merged_exons:
        x.source = "Stringtie"             
        x.score = "1000"
        x.frame = "."
        x["exon_number"] = str(counter)
        x["transcript_id"] = "MSTRG.1.1"
        x["gene_id"] = "MSTRG.1"
        gtf_file_load.add_relation(final_list[0], x, level = 1)
        final_list.append(x)
        counter += 1
    return final_list

def write_to_seperate_files(gtf_file_load, output_file):
    """Writes the non overlapping features to seperate files. These appear to
    be real regoins. The outputfile names will be stuctured in the way of 

    output_file_name_1
    output_file_name_2... etc...

    :iterative_parent_child: 
    :output_file): TODO
    :returns: TODO

    """

    def merge_exons_list(gtf_file_load, list_elements, counter):
        """TODO: Docstring for merge_exons_list.
        :list_elements: TODO
        :returns: TODO
        """

        final_list = []

        parent = list_elements[0]
        counter_element = "MSTRG." + str(counter)
        merged_exons = gtf_file_load.merge(list_elements[2:], ignore_strand = False)

        parent["transcript_id"] = counter_element + ".1"
        parent["gene_id"] = counter_element
    
        final_list.append(parent)
        exon_counter = 1
        for item in merged_exons:
            item.source = "Stringtie"             
            item.score = "1000"
            item.frame = "."
            item["exon_number"] = str(exon_counter)
            item["transcript_id"] = counter_element + ".1"
            item["gene_id"] = counter_element
            gtf_file_load.add_relation(final_list[0], item, level = 1)
            final_list.append(item)
            exon_counter += 1
        return final_list


    list_test = gtf_file_load.iter_by_parent_childs("transcript")
    counter = 1 

    for item in list_test:
        final_list = merge_exons_list(gtf_file_load,item, counter)
        output_file_name = str(output_file) + '_' + str(counter) + ".gtf"
        remove_file(output_file_name)
        with open(output_file_name, 'w') as fout:
            for features in final_list:
                fout.write(str(features) + '\n')
        counter += 1



def get_parser():
    parser = argparse.ArgumentParser(description='GTF transcript merge and\
            seperator. Purpose of this script is to take in a gtf file and\
            merge transcripts that are of the same isoform, as well as when\
            transcripts overlap, take the longest transcript.')
    parser.add_argument('-gtf', "--gtf_file", help="Bed file to load intervals. \
            These will be the intervals that will have their PKM value calcualte \
            ", required=True, dest='gtf')
    parser.add_argument('-o', "--output", help="Bed file to load intervals. \
            These will be the intervals that will have their PKM value calcualte \
            ", required=False, dest='o')

    args = vars(parser.parse_args())    
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    
    global database_random_name
    database_random_name = "gff_utils_db_" + (''.join(choice(ascii_uppercase) for i in range(12)))
    
    generate_directory(args.o)
    #check if file is empty. If it is, just touch the output
    file_val = check_empty_file(args.gtf)
    if file_val == 0 :
        output_file_name = str(args.o) + '_' + str(1) + ".gtf"
        open(output_file_name, 'a').close()
        sys.exit()
    elif file_val == 1:
        pass
    
    #Load GTF file if file is NOT empty
    loaded_gtf_file = read_gtf_file(args.gtf)
    multi_gene_return_val = check_multiple_genes(loaded_gtf_file)


    if multi_gene_return_val == 0:
       #Write the merged lncRNA to a file. Merged transcripts are good to go. 
        output_file_name = str(args.o) + '_' + str(1) + ".gtf"
        remove_file(output_file_name)
        
        #Basically merge all Exons and generate the final transcript.
        merged_exons = merge_exons(loaded_gtf_file)
        with open(output_file_name, 'w') as fout:
            for f in merged_exons:
                fout.write(str(f) + '\n')
    
    #At some point this should all be refactored to a series of functions.
    elif multi_gene_return_val == 1:
        #get all transcripts, count number

        check_val = check_overlaps(loaded_gtf_file)
        if check_val == 0:
            #no overlap, write all seperate 'genes' to files
            write_to_seperate_files(loaded_gtf_file, args.o)
        
        elif check_val == 1:
            #overlap, generate pairs to ID overlaps
            all_transcripts = get_feature("transcript", loaded_gtf_file)

            #Create pairs, try and find overlapping pairs
            for a,b in itertools.combinations(all_transcripts, 2):
                returned_val = check_merging_pairs(a,b,loaded_gtf_file)
                if returned_val == 1:
                    delete_this_tran = return_copy_to_delete(a,b,loaded_gtf_file)
                    #get_children_to_delete = loaded_gtf_file.children(delete_this_tran)
                    loaded_gtf_file.delete(delete_this_tran)
                else:
                    pass
            write_to_seperate_files(loaded_gtf_file, args.o)

