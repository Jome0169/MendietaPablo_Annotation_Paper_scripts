"""
Purpose of htis script is to iterate through each region of a bed file - take
the RNA-s1eq reads from ALL tissues types that align to this regions, and
assemble them. This way we can attempt to posotivly ID whether or not these
regions have putatative lncRNA sequences in them.
"""

import argparse
import sys
import os
import pysam 
import pybedtools

def read_bam_file(bam_file):
    """Takes in BAM files and reads it in using pysam

    :bam_file: TODO
    :returns: TODO

    """
    try:
        os.path.isfile(bam_file)
    except:
        FileNotFoundError

    bam_file_load = pysam.AlignmentFile(bam_file, 'rb')
    return bam_file_load

def read_bed_file(bed_file):
    """Takes in bed file and reads it in using pybedtools. Will make life
    easier longitudinally.

    :bed_file: TODO
    :returns: TODO

    """
    try:
        os.path.isfile(bed_file)
    except:
        FileNotFoundError
    bed_file_load = pybedtools.BedTool(bed_file)
    return bed_file_load

def read_bam_file_list(bam_file_list):
    """TODO: Docstring for read_bam_file.
    :returns: TODO

    """
    file_list = []

    with open(bam_file_list, 'r') as f:
        for line in f:
            bam_file = line.strip()
            file_list.append(bam_file)

    return file_list

def count_number_start_stop_reads(intersec_start, intersec_stop, bam_region):
    """Takes in start and stop of region. Then returns the number of posotive
    and negative reads that align to each region. 

    :(arg1: TODO
    :returns: TODO

    """

    def wrong_interval_intersect(intersect_start_pos, intersect_end_pos, read_pos,):
        """Takes in an interval and tests whether it falls at ALL into the region
        of intersection. In case the read is just a MASSIve split read, we don't
        really give a shit.
    
        :inersect_start: TODO
        :intersect_end: TODO
        :read_start: TODO
        :read_end: TODO
        :returns: TODO
    
        """
        if intersect_start_pos <= read_pos <= intersect_end_pos:
            return True
        else:
            return False
    
    final_list = [] 
    for read in bam_region:
        ciagar_string = read.cigarstring
        
        #Odd exception where sometimes reads seem to not generate iterable. In
        #that case break
        if ciagar_string == None:
            break
        else:
            pass

        #If read is split, test whether it falls wtihin intersection
        if "N" in ciagar_string:
            check_read_start = wrong_interval_intersect(intersec_start, intersec_stop, read.reference_start)
            check_read_end = wrong_interval_intersect(intersec_start, intersec_stop, read.reference_end)
            # IF read doesn't fall within interval, fuck it. 
            if check_read_start == True and check_read_end == True:
                final_list.append(read)

        #Read falls in intersection and not split
        elif "N" not in ciagar_string:
                final_list.append(read)
    return(final_list)


def generate_bam_file_load(item,region_string, bam_file):
    """TODO: Docstring for generate_bam_file_load.
    :returns: TODO

    """
    all_reads = []
    bam_load = read_bam_file(bam_file) 
    read_number = bam_load.count(region = region_string)
    fetch_region = bam_load.fetch(region=region_string)
    if read_number != 0:
        retrieve_reads = count_number_start_stop_reads(item.start, item.end, fetch_region)
        all_reads.append(retrieve_reads)
    elif read_number == 0:
        pass

    return all_reads


def bam_writer(bam_file_name, final_list):
    """TODO: Docstring for bam_writer.

    :bam_file_name: TODO
    :final_list): TODO
    :returns: TODO

    """
    #header = { 'HD': {'VN': '1.0'},'SQ': [{'LN': 1575, 'SN': 'chr1'},{'LN': 1584, 'SN': 'chr2'}] }
    header = { "HD" : {"VN":"1.5"}, \
            "SQ":[{"SN":"chr1","LN":308452471},  \
                {"SN":"chr2","LN":243675191}, \
                {"SN":"chr3","LN":238017767}, \
                {"SN":"chr4","LN":250330460}, \
                {"SN":"chr5","LN":226353449}, \
                {"SN":"chr6","LN":181357234}, \
                {"SN":"chr7","LN":185808916}, \
                {"SN":"chr8","LN":182411202}, \
                {"SN":"chr9","LN":163004744}, 
                {"SN":"chr10","LN":152435371}] } 

    output_file = pysam.AlignmentFile(bam_file_name, 'wb', header=header)
    for all_reads in final_list:
        for read in all_reads:
            for thing in read:
                output_file.write(thing)


def parse_bed_intersections(bed_load, bam_list, output_directory):
    """
    Takes in loaded BAM and bedtools file.

    For each bed intersetoin, iterates through, and ask the question, are the
    reads found in this region on the posotive strand, or the negative strand.
    After this, it returns a formatted list htat will then be written to an
    output file
    """

    for item in bed_load:
        all_reads = []
        location_string = item.chrom + ":" + str(item.start) + '-' + str(item.end)
        true_dif = item.end - item.start
        for bam_file in bam_list:
            reads_from_bam = generate_bam_file_load(item,location_string, bam_file)
            all_reads.append(reads_from_bam)
        create_region_name = output_directory + '/' + location_string.replace(":",'_') + ".bam"
        bam_writer(create_region_name, all_reads)





def read_bam_file_list(bam_file_list):
    """TODO: Docstring for read_bam_file_list.

    :bam_file_list: TODO
    :returns: TODO

    """
    final_list = []

    with open(bam_file_list, 'r') as f:
        for line in f:
            clean_line = line.strip()
            final_list.append(clean_line)
    return final_list

def generate_directory(dir_name):
    """TODO: Docstring for generate_directory.

    :dir_name: TODO
    :returns: TODO

    """
    try:
        os.mkdir(dir_name)
    except OSError:
        print("failed to create directory")

def get_parser():
    parser = argparse.ArgumentParser(description='Pull our reads aligning to a\
            region from multiple list of BAM files, puts them into a BAM file\
            for later assembly.')
    parser.add_argument('-bam', "--bam_file", help="Bam file list of files to \
            pull reads from.", required=True, dest='bam', # nargs="+", type=str\
            )
    parser.add_argument('-bed', "--bed_file", help="Bed file to load intervals. \
            These will be the intervals that will have their PKM value calcualte \
            ", required=True, dest='bed')
    parser.add_argument('-o', "--output", help="Bed file to load intervals. \
            These will be the intervals that will have their PKM value calcualte \
            ", required=False, dest='o')

    args = vars(parser.parse_args())    
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()

    generate_directory(args.o)
    bam_files_to_iter = read_bam_file_list(args.bam)
    bed_file = read_bed_file(args.bed)

    parse_bed_intersections(bed_file, bam_files_to_iter, args.o)


    #Calc things. This value is passed the thing you are calculating

    #if args.o != None:
    #    calculation.saveas(args.o)
    #elif args.o == None:
    #    pass
