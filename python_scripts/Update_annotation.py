"""
Purpose of this script is to update the annotation of genes using chromatin
mods. Many of the annotations in maize appear to be lacking, so this is simply
a way of updating these annotations in a manner to see how many more "mystery"
peaks we actually end up picking.  
Note - this script utilizes pybedtools a TON. If you're not familiar with this
package i would highly reccomend reading the DOCs
"""

#Mon Aug 26 17:14:26 EDT 2019
#Current Issue Reigion 2:197184803..197245902 (61.1 Kb)

import argparse
import sys
import os
import pysam 
import pybedtools
from pybedtools.featurefuncs import extend_fields
from pybedtools.featurefuncs import greater_than

def read_bed_file(bed_file, type_of_bed):
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

    if type_of_bed == "annot":
        bed_file_load = pybedtools.BedTool(bed_file).sort().each(keep_nfeatures, 6).saveas()
    elif type_of_bed == "hist":
        bed_file_load = pybedtools.BedTool(bed_file).sort().each(keep_nfeatures, 3).saveas()
    
    return bed_file_load


def merge_intersect_features(feature, number_index_1, score_action, string_name):
    """Instead of hardcoding the merge based off the index, instead take
    the length of each list, and identicy the start and stop regions based
    off the length of the original list

    """
    
    #Update the score values based off the flag given
    if score_action == "add":
        score_off_index = feature[int(number_index_1) + 4]
        old_score = int(feature.score)
        new_score = int(old_score) + int(score_off_index)
        feature.score = int(new_score)
    elif score_action == "replace":
        score_off_index = feature[int(number_index_1) + 4]
        new_score = score_off_index
        feature.score = int(new_score)
    elif score_action == None:
        pass

    #gene_strand_info = feature[11]
    take_start_index_intersect = int(number_index_1) + 1
    take_stop_index_intersect = int(number_index_1) + 2

    take_starts =  [feature.start, int(feature[take_start_index_intersect])]
    take_stops = [feature.stop, int(feature[take_stop_index_intersect])]

    take_new_start = min(take_starts)
    take_new_stop = max(take_stops)

    #Generate the New Feature
    feature.start = take_new_start
    feature.stop = take_new_stop
    feature.name = string_name
    return feature


def keep_nfeatures(feature, n):
    """Given a list of features, keep only 
    n number of fields

    """
    new_feature = feature[:int(n)]
    feature = new_feature

    return feature


def replace_feature_name(feature, string_name):
    """ Given a bedtool feature, go through and replace the feature name. All
    Features will have the same name.
    """
    feature.name = string_name
    return feature

def alter_score(feature, action, val):
    """Alters the score feature of bedtool item. Will wither replace, or sum
    depending on the 'action' value given

    """
    if action == 'r':
        feature.score = str(val)
    elif action == 'a':
        new_score = int(feature.score) + val
        feature.score = new_score
    return feature


def merge_histone_marks(*args):
    """Takes two peak files of multiple histone mods, intersects the set,
    merging those that intersect, and adding "2" to the score column. This
    alue will be used in further analysis to assign the "confidence" of a
    histone modification overlapping a lncRNA

    Expects 3 args MAX

    hist_mod1: Histone mod bedTOols 
    hist_mod2: BedTool of the correlating histone mod
    mod_name: String, base name being affixed to end of the bedfile

    """
    #Given only a single histone mod file corresponding to those regions, read
    #it in, give all peaks a corresponding score of 1
    if len(args) == 2:
        mod_name = args[1]
        hist_mod_1 = read_bed_file(args[0], "hist").merge(d=100)
        hist_mod_rename = hist_mod_1.each(extend_fields, 5).each(alter_score, "r", \
                        1).each(replace_feature_name, mod_name)
        final_set = hist_mod_rename.saveas()


    elif len(args) == 3:
        #Record mod type going at the end of the file
        #REad in the BED files
        mod_name = args[2]
        hist_mod_1 = read_bed_file(args[0], "hist").merge(d=100)
        hist_mod_2 = read_bed_file(args[1], "hist").merge(d=100)

        #Take lengths for later use
        hist_len = hist_mod_1.field_count()

        #Intersect keeping values
        intersecting_hist_mods = hist_mod_1.intersect(hist_mod_2, wa = True, wb =
                True)

        #Merge those regions which intersect, giving them an inital valyue of 2 in
        #their score category, filter down to only 5 locations
        intersecting_hist_mods_merged = intersecting_hist_mods.each(merge_intersect_features, hist_len, None, mod_name).saveas()
        intersecting_hist_mods_merged_altered_score = intersecting_hist_mods_merged.each(alter_score, "r", 2).saveas()
        intersecting_hist_mods_merged_trimmed = intersecting_hist_mods_merged_altered_score.each(keep_nfeatures, 5).each(extend_fields, 6).saveas()

        #Take regions that are unique to each, and assing them the score of '1' for
        #less confidence overall. Trim these to 5 columns, same as above, assign score as 1
        unique_to_hist_mod_1 = hist_mod_1.intersect(hist_mod_2, wa = True, v = True).saveas()
        unique_to_hist_mod_1_extended_field = unique_to_hist_mod_1.each(extend_fields, 6).saveas()
        unique_to_hist_mod_1_extended_field_updated_score = unique_to_hist_mod_1_extended_field.each(alter_score, "r", 1).saveas()
        unique_to_hist_mod_1_final = unique_to_hist_mod_1_extended_field_updated_score.each(replace_feature_name, mod_name).saveas()


        unique_to_hist_mod_2 = hist_mod_2.intersect(hist_mod_1, wa = True, v = True).saveas()
        unique_to_hist_mod_2_extended_field = unique_to_hist_mod_2.each(extend_fields, 6).saveas()
        unique_to_hist_mod_2_extended_field_updated_score = unique_to_hist_mod_2_extended_field.each(alter_score, "r", 1).saveas()
        unique_to_hist_mod_2_final = unique_to_hist_mod_2_extended_field_updated_score.each(replace_feature_name, mod_name).saveas()
        
        #Merge all the regions together
        final_set = intersecting_hist_mods_merged_trimmed.cat(*[unique_to_hist_mod_1_final,unique_to_hist_mod_2_final],
                postmerge = False,force_truncate=False).sort().saveas()
        
    else:
        print("ERROR, TOO MANY ARGUMENTS PASSED")
        sys.exit(-1)
    
    
    #return the final set for use
    return(final_set)
    
    


def filter_genes_with_no_expression(gene_bed_list, RNA_seq, strand):
    """We are getting some weird marks due to histone mods overlapping genes
    that are overlapping other genes, causing incorret annotations. THis way
    we're filtereing at least for the genes taht are acutally ON

    :gene_bed_list: TODO
    :RNA_seq): TODO
    :returns: TODO

    """
    if strand == True:
        #Count number of Reads
        gene_annotation_expressed = gene_bed_list.intersect(RNA_seq, S=True, c=True)
        gene_annotation_expressed.saveas()


        #Filter for greater than 10 reads
        filtered_gene_annotation = gene_annotation_expressed.filter(lambda x: \
                int(x[-1]) > 5 ).saveas().each(remove_last_nfeature, 1).saveas()
    elif strand == False:
        #Count number of Reads
        gene_annotation_expressed = gene_bed_list.intersect(RNA_seq, c=True)
        
        #Filter for greater than 10 reads
        filtered_gene_annotation = gene_annotation_expressed.filter(lambda x: \
                int(x[-1]) > 5 ).saveas().each(remove_last_nfeature, 1).saveas()

       #Assign sign strandedness here. Since it's intersecting reads we can
       #count number of reads from each strand, take greater number, and assign

    return(filtered_gene_annotation)





def remove_last_nfeature(feature, feature_num=0):
    """Removes the laste feature found in each object. In this case
    stiping the RNA-seq counts as it's altering downstream things

    :feature: TODO
    :returns: TODO

    """
    feature = feature[:-feature_num]
    return feature


def generate_directory(dir_name):
    """TODO: Docstring for generate_directory.

    :dir_name: TODO
    :returns: TODO

    """
    try:
        os.mkdir(dir_name)
    except OSError:
        print("failed to create directory")


def find_closest_broad_with_correct_promoter(gene_annotation, peak_of_type):
        """Find the closes promoter to each gene. Initally split by strand
        because that dictates what 'upsteam' and downstream are 
        :returns: TODO

        """
        #Upstream/downstream is different and shit
        posotive_genes = pybedtools.BedTool(x for x in gene_annotation if x.strand == '+')
        negative_genes = pybedtools.BedTool(z for z in gene_annotation if z.strand == '-')

        gene_closest_promoter_pos = posotive_genes.closest(peak_of_type, id=True, D="a", t="first")
        gene_closest_promoter_neg = negative_genes.closest(peak_of_type, id=True, D="a", t="last")


        combined_genes = gene_closest_promoter_pos.cat(gene_closest_promoter_neg, postmerge=False).sort()

        return(combined_genes)

 
def merge_gene_feature_promoter(feature, number_index_1):
    """BEDTOOLS FUNCTION - takes in a featurea and takes the location of the
    promoter sequence farthest - after this boots back a new feature with
    different start sites based off of promoter mark
    """

    def select_new_start_stop(list_val, start_or_stop):
        """This was throwing an issues earlier as we had to defend against values
        being completely off the genomes
        """
        if start_or_stop == "start":
            #Take the min value for the start sequence
            take_min_value = min(list_val)
            if take_min_value == 0 or take_min_value == -1:
                take_min_value = max(list_val)
            return take_min_value
        elif start_or_stop == "stop":
            #Take the Max value for where it's stopiing
            take_min_value = max(list_val)
            return take_min_value

    #Take the absolute value of distance, and use that value to assign which
    #'category'  this is
    strand_info = feature.strand
    take_total_distance = int(feature[-1])
    #Twice the value of an average exon
    if abs(take_total_distance) >= 5000 and abs(take_total_distance) <= 20000:
        feature.name = feature.name + '_major_extended_hyper_large'
        feature.score = feature[number_index_1 + 4]
    elif abs(take_total_distance) >= 500 and abs(take_total_distance) < 5000:
        feature.name = feature.name + '_major_extended'
        feature.score = feature[number_index_1 + 4]
    elif abs(take_total_distance) >= 100 and abs(take_total_distance) < 500:
        feature.name = feature.name + '_minor_extended'
        feature.score = feature[number_index_1 + 4]
    elif abs(take_total_distance) < 100:
        feature.score = feature[number_index_1 + 4]
    else:
        pass


    take_start_index_intersect = int(number_index_1) + 1
    take_stop_index_intersect = int(number_index_1) + 2

    #Adjust stop for negative gene
    if strand_info == '-':
        #Defense against last gene with no promoter till end of genome
        id_new_starts = [feature.start, int(feature[take_start_index_intersect])] 
        id_new_stops = [feature.stop, int(feature[take_stop_index_intersect])]

        feature.start = select_new_start_stop(id_new_starts, "start")
        feature.stop = select_new_start_stop(id_new_stops, "stop")

    #Adjust start for posotive gene
    elif strand_info == '+':

        #Defense against last gene with no promoter till end of genome
        id_new_starts = [feature.start, int(feature[take_start_index_intersect])] 
        id_new_stops = [feature.stop, int(feature[take_stop_index_intersect])]
        
        feature.start = select_new_start_stop(id_new_starts, "start")
        feature.stop = select_new_start_stop(id_new_stops, "stop")

    return feature



def merge_unkown_peaks_gene(feature):
    """

    :feature: TODO
    :returns: TODO

    """
    #gene_strand_info = feature[11]

    take_starts =  [feature.start, int(feature[11])]
    take_stops = [feature.stop, int(feature[12])]

    take_new_start = min(take_starts)
    take_new_stop = max(take_stops)

    #Alter Names
    add_start_stop_diff = gene_extension_val(take_starts, take_stops)
    if add_start_stop_diff >= 20000:
        return feature
    elif add_start_stop_diff >= 5000 and add_start_stop_diff < 20000:
        feature.name = feature.name + '_major_extended_hyper_large'
        feature.score = "3"
    if add_start_stop_diff >= 500 and  add_start_stop_diff < 5000:
        feature.score = "2"
        feature.name = feature.name + '_major_extended'
    elif add_start_stop_diff < 500 and add_start_stop_diff > 50:
        feature.score = "1"
        feature.name =  feature.name + '_minor_extended'

    return feature

#def merge_unkown_peaks_gene2(feature):
#    """
#
#    :feature: TODO
#    :returns: TODO
#
#    """
#
#    take_starts =  [feature.start, int(feature[7])]
#    take_stops = [feature.stop, int(feature[8])]
#
#    take_new_start = min(take_starts)
#    take_new_stop = max(take_stops)
#
#    add_start_stop_diff = gene_extension_val(take_starts, take_stops)
#
#    if add_start_stop_diff >= 20000:
#        return feature
#    elif add_start_stop_diff >= 5000 and add_start_stop_diff < 20000:
#        feature.name = feature.name + '_major_extended_hyper_large'
#        feature.score = "3"
#    if add_start_stop_diff >= 500 and  add_start_stop_diff < 5000:
#        feature.score = "2"
#        feature.name = feature.name + '_major_extended'
#    elif add_start_stop_diff < 500 and add_start_stop_diff > 50:
#        feature.score = "1"
#        feature.name =  feature.name + '_minor_extended'
#
#    return feature



def find_closest_gene(isoalted_peaks_no_promoter, merged_promoter_gene_features):
        """Find the closes promoter to each gene. Initally split by strand
        because that dictates what 'upsteam' and downstream are 
        :returns: TODO

        """

        def merge_unkown_peaks_gene(feature):
            """TODO: Docstring for merge_unkown_peaks_gene.

            :feature: TODO
            :returns: TODO

            """
            #gene_strand_info = feature[11]

            take_starts =  [feature.start, int(feature[6])]
            take_stops = [feature.stop, int(feature[7])]

            take_new_start = min(take_starts)
            take_new_stop = max(take_stops)


            #Generate the New Feature
            feature.start = take_new_start
            feature.stop = take_new_stop
            feature.strand = feature[10]
            feature.name = feature[8] + '_extended'
            feature.score = '3'
            return feature

        closest_grouping = isoalted_peaks_no_promoter.sort().closest(merged_promoter_gene_features,
                D='a').sort().each(merge_unkown_peaks_gene).each(score_value_edit,'a',1 )

        return(closest_grouping)

#def find_closest_gene2(isoalted_peaks_no_promoter, conserved_promoter):
#        """Find the closes promoter to each gene. Initally split by strand
#        because that dictates what 'upsteam' and downstream are 
#        :returns: TODO
#
#        """
#
#        def merge_unkown_peaks_gene(feature):
#            """TODO: Docstring for merge_unkown_peaks_gene.
#
#            :feature: TODO
#            :returns: TODO
#
#            """
#            #gene_strand_info = feature[11]
#
#            take_starts =  [feature.start, int(feature[6])]
#            take_stops = [feature.stop, int(feature[7])]
#
#            take_new_start = min(take_starts)
#            take_new_stop = max(take_stops)
#
#
#            #Generate the New Feature
#            feature.start = take_new_start
#            feature.stop = take_new_stop
#            feature.strand = "."
#            feature.name = "Possible_novel_gene"
#            feature.score = '4'
#            return feature
#
#        closest_grouping = isoalted_peaks_no_promoter.sort().closest(conserved_promoter,
#                D='a').sort().each(merge_unkown_peaks_gene).sort().saveas()
#
#        merged_closest = closest_grouping.merge(c=[4,5,6,9],
#                o=["distinct,distinct,distinct,min"])
#
#        return(merged_closest)

def assign_strand_novel(feature):
    """Assign strand based off which direction promoter was from region of
    interest


    """
    diff_value = feature[-1]
    if int(diff_value) > 0:
        feature.strand = '-'
    elif int(diff_value) < 0:
        feature.strand = '+'
    return feature

def assign_strand_orphan_peaks(bed_file, BAM_file):
    """
    Some of our regions appear to have segments that belong to one region
    rather than the other, so for instance


    gene:         $$$$----$$$$$>             <$$$$----$$$$$
    H3K36me3:     ###############      ####################
    RNA_Reads:    -------------------------- ++++++++++++++ 

    So in this example the region overlapping the second gene has reads that
    are assigning to the other gene

    ProblemTypeRegions: 1       12103819        12114999


    Note: This script utilized Pysam as pybedtools didn't want to play nice
    with this

    :returns: TODO
    """

    def add_strand(f):
        # Omitting this extend_fields call, or using a field count < 6, causes
        # a segfault for BedTool `b`.
        f = extend_fields(f, 6)
        f.strand = '+'
        return f
       
    def set_true_strand(feature):
        """uses the below intersection to identify the true strandedness of the
        library

        Note - developed for ANTI-SENSE library f-R

        :arg1: TODO
        :returns: TODO

        """

        posotive_strand = int(feature[-2])
        negative_strand = int(feature[-1])

        if posotive_strand > negative_strand:
            feature.strand = "-"
        elif posotive_strand < negative_strand:
            feature.strand = '+'
        return feature
 
    def alter_region_read_support(feature):
        """Takes in the feature, and replaces the BED object with the region
        that's only supported by reads

        :feature: TODO
        :returns: TODO

        """
        new_start = int(feature[-2])
        new_stop = int(feature[-1])

        feature.start = new_start
        feature.stop= new_stop
        return feature


    #Count reads, remove regions with less than 20 reads
    count_reads_intersecting = bed_file.intersect(BAM_file,
        c = True).filter(lambda x: int(x[-1]) > 5).each(remove_last_nfeature, 1).saveas()


    #Add fake strand to avoid seg-fault
    set_temp_strand = count_reads_intersecting.each(add_strand)

    #First is same strand, so posotive, which is a gene on the ANTI-sense strand
    #Second is OPPOSITE STRAND, so NEGATIVE, genes on SENSE strand
    add_posotive_read_count = set_temp_strand.intersect(BAM_file, s = True, c=True).saveas()
    add_negative_read_count = add_posotive_read_count.intersect(BAM_file, S = True, c=True).saveas()

    #Take intersection from above and assign true strandedness
    true_strand_assignment = add_negative_read_count.each(set_true_strand).sort().each(remove_last_nfeature, 2).saveas()

    #Intersect, then take the region supportede by RNA_seq data
    #NOTE THE INTERSECT HERE IS OPPOSITE STRAND DUE TO THE SEQUENCING LIBRARY
    #TYPE. At some point this will need to be modified to be more accomidating
    take_longest_supported_region = true_strand_assignment.intersect(BAM_file,
            wa = True, wb = True, S = True).saveas()
    
    
    if len(take_longest_supported_region) == 0:
        return None
    else:
        pass

    take_longest_supported_region_group = take_longest_supported_region.groupby(g="1-6", c="8,9", 
            o="min,max").each(alter_region_read_support).saveas()
    clean_strand = take_longest_supported_region_group.each(remove_last_nfeature, 2).saveas()


    return clean_strand


def find_closest_promoter(gene_body_peaks, promoter_regions):
        """Find the closes promoter to each gene. Initally split by strand
        because that dictates what 'upsteam' and downstream are 
        :returns: TODO

        """

        def take_first_20(feature):
            """Given a list of bed intervals, this takes in the list, and reports
            back the starting 20% of that region. 
            """

            #Original Sequence + 20 percent. The gene TSS should be in here
            if feature.strand == '+':
                feature_length = feature.stop - feature.start
                percent_20 = round(.2 * feature_length)
                new_stop = feature.start + percent_20

                feature.start = feature.start
                feature.stop = new_stop

            elif feature.strand == '-':
                feature_length = feature.stop - feature.start
                percent_20 = round(.2 * feature_length)
                new_stop = feature.stop - percent_20
                feature.start= new_stop
            return feature


        closest_sort = gene_body_peaks.sort().saveas()
        field_number = closest_sort.field_count()

        closest_sort_top_20 = gene_body_peaks.each(take_first_20).saveas()
        
        posotive_regions = pybedtools.BedTool(x for x in closest_sort_top_20 if x.strand == '+')
        negative_regions = pybedtools.BedTool(z for z in closest_sort_top_20 if z.strand == '-')

        closest_promoter_pos = posotive_regions.closest(promoter_regions,
                id=True, D="a", t="first")
        closest_promoter_neg = negative_regions.closest(promoter_regions,
                id=True, D="a", t="last")

        combined_promoter_regions = closest_promoter_pos.cat(closest_promoter_neg, \
                postmerge=False).sort().filter(lambda x: int(x[-1]) <= 20000).saveas()

        closest_grouping = combined_promoter_regions.each(merge_intersect_features, field_number,
                "add", "Possible_novel_gene").filter(lambda \
                x: x.start > 0 and x.stop > 0 ).saveas()

        #Filter step for incorrect inputes like 3       -1      235587066       Possible_novel_gene
        closest_grouping_filtered = closest_grouping.remove_invalid()

        #LEFTOVER REGIONS. UNSURE WHAT THE PURPOSE HERE WAS
        #original_plus_possible_promter = gene_body_peaks.intersect(closest_grouping_filtered, 
        #        wa = True, wb = True, s = True)#.each(merge_unkown_peaks_gene).filter(lambda x: \
                        #x.start > 0 and x.stop > 0).saveas()

        return(closest_grouping_filtered)


def gene_extension_val(start_list, stop_list):
    """Generates scores for merger of gene regions.


    :start_list: TODO
    :stop_list: TODO
    :returns: TODO

    """
    start_diff = abs(start_list[0] - start_list[1])
    stop_diff = abs(stop_list[0] - stop_list[1])

    add_start_stop_diff = start_diff + stop_diff
    return add_start_stop_diff

def merge_regions(feature, idx_other_strt, idx_other_stop):
    """TODO: Docstring for merge_unkown_peaks_gene.

    :feature: TODO
    :returns: TODO

    """

    take_starts =  [feature.start, int(feature[idx_other_strt])]
    take_stops = [feature.stop, int(feature[idx_other_stop])]

    take_new_start = min(take_starts)
    take_new_stop = max(take_stops)
    
    #Alter Names
    add_start_stop_diff = gene_extension_val(take_starts, take_stops)

    if add_start_stop_diff >= 20000:
        return feature
    elif add_start_stop_diff >= 5000 and add_start_stop_diff < 20000:
        feature.name = feature.name + '_major_extended_hyper_large'
        feature.score = "3"
    if add_start_stop_diff >= 500 and  add_start_stop_diff < 5000:
        feature.score = "2"
        feature.name = feature.name + '_major_extended'
    elif add_start_stop_diff < 500 and add_start_stop_diff > 50:
        feature.score = "1"
        feature.name =  feature.name + '_minor_extended'

    #Generate the New Feature start and stop
    feature.start = take_new_start
    feature.stop = take_new_stop

    return feature


def get_parser():
    parser = argparse.ArgumentParser(description='Pull our reads aligning to a\
        region from multiple list of BAM files, puts them into a BAM file\
        for later assembly.')
    #parser.add_argument('-H3K4me3', "--H3K4me3_file", help="Bam file list of files to \
    #    pull reads from.", required=True, dest='H3K4me3')
    #parser.add_argument('-H3K56ac', "--H3K56ac_file", help="Bam file list of files to \
    #    pull reads from.", required=True, dest='H3K56ac')
    parser.add_argument('-broad', "--broad_peaks", help="Bam file list of files to \
        pull reads from.", required=True, dest='broad', nargs='+',)
    parser.add_argument('-narrow', "--narrow_peaks", help="Bam file list of files to \
        pull reads from.", required=True, dest='narrow', nargs='+',)
    #parser.add_argument('-H3K36me3', "--H3K36me3_file", help="Bam file list of files to \
    #    pull reads from.", required=True, dest='H3K36me3')
    #parser.add_argument('-H3K4me1', "--H3K4me1_file", help="Bam file list of files to \
    #    pull read from.", required=False, dest='H3K4me1')
    parser.add_argument('-annotation', "--gene_annotation", help="Bed file to load intervals. \
        These will be the intervals that will have their PKM value calcualte \
        ", required=True, dest='annot')
    parser.add_argument('-RNA', "--RNA_data", help="BAM file that was hopefully \
        generate at the same time as the CHIP_seq data ", required=True, dest='bam')
    parser.add_argument('-lncRNA', "--lncRNA_file", help="File of lncRNAs to \
            remove novels with", required=True, dest='lnc')
    parser.add_argument('-o', "--output", help="Bed file to load intervals. \
            These will be the intervals that will have their PKM value calcualte \
            ", required=False, dest='o')

    args = vars(parser.parse_args())
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    
    #Load all Bed files
    print("Loading Bed Files")
    print("Taking intersections between broad and narrow marks")
    promoter_marks_merged = merge_histone_marks(*args.narrow, "promoter_marks")
    gene_body_marks_merged = merge_histone_marks(*args.broad, "gene_body_marks")
    
    annotation_bed = read_bed_file(args.annot, "annot")
    lncRNA_bed = read_bed_file(args.lnc, "annot")
    RNA_seq_file = pybedtools.BedTool(args.bam)

    #Look for genes with at least 5 RNA-seq reads aligning
    print("Filtering for Expressed Genes")
    expressed_gene_list = filter_genes_with_no_expression(annotation_bed, RNA_seq_file, True).saveas()
    expressed_lncRNA_list = filter_genes_with_no_expression(lncRNA_bed, RNA_seq_file, True).saveas()

    #We Only Care about genes that are 'on'
    #I'm unsure if this is really the best way of going about this. Should I
    #really be using stranded RNA-seq data to dictate this?
    print("Gathering Genes Overlapping Gene Body Modifications")
    gene_overlapping_marks = expressed_gene_list.intersect(gene_body_marks_merged, wa = True, wb = True, u=True, f=.1)

    #Intersect broad peaks with and w/o promoter marks 
    #Find closest promoter to each gene, retruns gene feature with promoter
    #like feature next to it
    print("Intersecting Broad peaks with Promoter Peaks")
    genes_with_promoter_pairs = find_closest_broad_with_correct_promoter(gene_overlapping_marks, promoter_marks_merged)

    #Merge the closest promoter to the gene feature, add value     
    print("Merging closest promoter to gene features")
    gene_len = gene_overlapping_marks.field_count()
    gene_promoters_merged_raw = genes_with_promoter_pairs.each(merge_gene_feature_promoter, gene_len).sort().saveas()

    #Remove the overlapping region information
    gene_promoters_merged = gene_promoters_merged_raw.each(remove_last_nfeature,7).saveas()
    
    #Same thing as above, but do this with lncRNAs
    print("Checking Previous lncRNAs for Incorrect Annotation")
    lncRNA_overlapping_marks = expressed_lncRNA_list.intersect(gene_body_marks_merged, wa = True, wb = True, u=True, f=.1)
    lncRNAs_with_promoter_pairs = find_closest_broad_with_correct_promoter(lncRNA_overlapping_marks, promoter_marks_merged)
    lncRNA_len = lncRNA_overlapping_marks.field_count()
    lncRNA_promoters_merged_raw = lncRNAs_with_promoter_pairs.each(merge_gene_feature_promoter, lncRNA_len).sort().saveas()
    lncRNA_promoters_merged = lncRNA_promoters_merged_raw.each(remove_last_nfeature,7).saveas()

    correct_lncRNAs = lncRNA_promoters_merged.intersect(gene_promoters_merged,
            wa = True, v = True, s = True).saveas("Correct_lncRNA.bed")
    incorrect_lncRNAs_overlap_gene = lncRNA_promoters_merged.intersect(gene_promoters_merged,
            wa = True, s = True).saveas("Incorrect_lncRNA.bed")


    #Since these are broad peak regions with no promoter, these are probably
    #regions that are additional LENGTH areas of the gene.
    print("Finding parent of Orphan Exression Regions")
    broad_peaks_no_gene_intersections = gene_body_marks_merged.intersect(gene_promoters_merged, v=True, wa = True).saveas()

    #Sometimes this will return None if no novel regions are found. 
    novel_peaks_stranded = assign_strand_orphan_peaks(broad_peaks_no_gene_intersections, RNA_seq_file)
    if novel_peaks_stranded == None:
        possible_novel_gene = None
    else:
        novel_peaks_stranded.saveas()
        #Merge with upstream promoter element
        possible_novel_gene_with_promoter = find_closest_promoter(novel_peaks_stranded, \
                promoter_marks_merged).sort().each(remove_last_nfeature, \
                        4).saveas()

        #Remove any overlapping features that already overal a known lncRNA
        possible_novel_gene_with_promoter_no_overlap_lncRNA = possible_novel_gene_with_promoter.intersect(correct_lncRNAs, v=True, wa=True, s=True)

        #Remove any overlapping features that already overal a known GENE 
        possible_novel_gene = possible_novel_gene_with_promoter_no_overlap_lncRNA.intersect( \
                gene_promoters_merged, v=True, wa=True, s=True).saveas()


    print("Merging Orphan Exression Regions")
    #Sometimes these possible novel genes overlap a known gene and are just an
    #extenstion, so if they are, add them.
    if possible_novel_gene != None:
        gene_index_count = gene_promoters_merged.field_count()
        gene_promoters_merged_overlapping_extension_region = gene_promoters_merged.intersect(possible_novel_gene, 
            wa = True, wb= True).sort().each(merge_gene_feature_promoter, gene_index_count).saveas()
    else:
        #If we find no novel genes, just pass
        gene_promoters_merged_overlapping_extension_region = gene_promoters_merged

    #Get list of genes with NO extension
    unmodified_genes = gene_promoters_merged.intersect(gene_promoters_merged_overlapping_extension_region,
            wa = True, v = True, s = True).sort().saveas()
    
    #Add extended and un-extended genes back into same list 
    modified_and_unmodified_genes = unmodified_genes.cat(gene_promoters_merged_overlapping_extension_region,
            postmerge = False).sort().saveas()
    
    ##Generate list of novel like gene
    if possible_novel_gene == None:
        #If no novel genes found 
        print("No Novel Transcript Found")
        updated_gene_promoters_merged = modified_and_unmodified_genes.sort().saveas()
    else:
        #If novel Genes found, ensure that they don't overlap known genes or
        #known ncRNAs
        print("Novel Transcript Identification")
        novel_gene_list = possible_novel_gene.intersect(gene_promoters_merged_overlapping_extension_region,
            wa = True, v = True, s = True).intersect(correct_lncRNAs, wa = True, 
                    v = True, s= True).saveas()

    #Add possible novel genes to list of extended and non-extended genes
        updated_gene_promoters_merged = modified_and_unmodified_genes.cat(novel_gene_list, 
                postmerge=False).sort().saveas()

    #Assign strand orientation to the abandoned region that don't have strand
    #based off reads. Allowing us to merge it to it's correct features

    #Remove small regions
    left_over_regions = gene_body_marks_merged.subtract(updated_gene_promoters_merged).filter(greater_than,
            2000).saveas()

    #Give strand to leftover regions that HAVE RNA_seq support. Filter again,
    #if not greater than an exon, not worth the time fixing downstream. Note
    #that we're only inputting RNA that is not already overlapping a gene. I
    #was having issues where single base pair overlaps were causing poor
    #annotation/assignment

    RNA_not_overlapping_gene = RNA_seq_file.intersect(updated_gene_promoters_merged,
            S=True, v=True, wa = True)
   
    assigned_strand_left_over_regions = assign_strand_orphan_peaks(left_over_regions,
            RNA_not_overlapping_gene).filter(greater_than,
                    500).intersect(correct_lncRNAs, wa = True, v = True, s = True).saveas()

    #Find closest promoter for leftover regions
    left_overs_with_promoters = find_closest_promoter(assigned_strand_left_over_regions,
            promoter_marks_merged).intersect(correct_lncRNAs, wa = True, v = True, s = True).saveas() 

    print("Identifying Novel genes")
    #Identify gene overlap of unknonw regions
    genes_merged_with_leftovers = updated_gene_promoters_merged.intersect(left_overs_with_promoters, 
            wa = True, wb = True, s = True).each(merge_regions, 7,
                    8).saveas()
    
    print("Merging Novel genes") 
    z = left_overs_with_promoters.intersect(genes_merged_with_leftovers, 
            wa = True, v = True).saveas().intersect(annotation_bed,
                    wa = True, v = True, s = True).saveas()

    print("Producing Final List") 
    final_set = genes_merged_with_leftovers.cat(z, postmerge = False)
    remove_old_features = updated_gene_promoters_merged.intersect(final_set, wa = True, f = 1.00, v = True, s = True)
    final_updated_list = remove_old_features.cat(final_set,  postmerge = False).sort()
    
    print("Filtering Step")
    final_updated_list_sense = final_updated_list.filter(lambda x: x.strand == "+").saveas()
    final_updated_list_antisense = final_updated_list.filter(lambda x: x.strand == "-").saveas()


    print("Merging Identical Outputs Sense")
    final_updated_list_sense_merged = final_updated_list_sense.groupby(g="1-2",
            c="3,4,5,6", o="max,collapse,mean,first").saveas()

    print("Merging Identical Outputs anti-Sense")
    final_updated_list_antisense_merged1 = final_updated_list_antisense.groupby(g="1-2", 
            c="3,4,5,6", o="max,collapse,mean,first").saveas()
    final_updated_list_antisense_merged2 = final_updated_list_antisense_merged1.merge(s = True, c='4,5,6',o="distinct,mean,distinct")

    print("Generating Final List")
    final_list = final_updated_list_sense_merged.cat(final_updated_list_antisense_merged2,
            postmerge=False).sort().saveas(args.o)

    print("Done")
