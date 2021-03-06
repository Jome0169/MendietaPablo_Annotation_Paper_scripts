import os 
from sys import argv



generate_list = []
with open(argv[1]) as f:
    for line in f:
        clean_line = line.strip().split()
        generate_list.append(clean_line)




def count_number_novel(seperate_list):
    """TODO: Docstring for count_number_novel.

    :seperate_list: TODO
    :returns: TODO

    """
    counter = 0
    for i in seperate_list:
        if "novel" in i:
            counter += 1
    return counter


list_of_novel = []
list_of_merged = []
list_of_minor_extension = []
list_of_major_extension = []
list_of_un_altered= []
list_of_hyper_large = []



for item in generate_list:
    if "," in item[3] and str(argv[2]) in item[3]:
        split_name_list = item[3].split(',')
        number_novel = count_number_novel(split_name_list)
        number_merged_items = len(split_name_list)
        take_diff = (number_merged_items - number_novel)
        if take_diff == 2:
            list_of_major_extension.append(item)
        elif take_diff == 0:
            list_of_novel.append(item)
        elif take_diff > 1:
            list_of_merged.append(item)
    elif "," in item[3] and  str(argv[2]) not in item[3]:
        list_of_novel.append(item)
    elif "," not in item[3] and "novel" in item[3]:
        list_of_novel.append(item)
    elif "," not in item[3] and "hyper_large" in item[3]:
        list_of_hyper_large.append(item)
    elif "," not in item[3] and "major_extended" in item[3] and "hyper_large" not in item[3]:
        list_of_major_extension.append(item)
    elif "," not in item[3] and "minor_extended" in item[3] and "hyper_large" not in item[3]:
        list_of_minor_extension.append(item)
    else:
        list_of_un_altered.append(item)




def write_output(list_of_bed, output_file):
    """TODO: Docstring for write_output.
    :returns: TODO

    """
    try:
        os.remove(output_file)
    except OSError:
        pass
    with open(output_file, 'a+') as f:
        for item in list_of_bed:
            tab_sep = '\t'.join(item)
            f.write(tab_sep)
            f.write('\n')

write_output(list_of_novel, argv[3] + '_' + "novel_genes.bed")
write_output(list_of_merged, argv[3] + '_' + "merged_genes.bed")
write_output(list_of_minor_extension, argv[3] + '_' + "minor_extension_genes.bed")
write_output(list_of_major_extension, argv[3] + '_' + "major_extension_genes.bed")
write_output(list_of_un_altered, argv[3] + '_' + "un_altered_genes.bed")
write_output(list_of_hyper_large, argv[3] + '_' + "hyper_large_genes.bed")
