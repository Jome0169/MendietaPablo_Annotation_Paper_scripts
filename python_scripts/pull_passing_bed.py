from sys import argv


def read_bed_file_pull_passing(arg1, name_list):
    """TODO: Docstring for read_bed_file_pull_passing.

    :arg1: TODO
    :returns: TODO

    """
    passing_lines = []
    with open(arg1, 'r') as f:
        for line in f:
            clean_line = line.strip().split()
            if clean_line[3] in name_list:
                passing_lines.append(clean_line)
    return(passing_lines)

def grab_original_annotations(annot_list):
    """TODO: Docstring for grap_original_annotations.
    :returns: TODO

    """
    original_annotations = []
    merged_annotations_split = []
    for item in annot_list:
        if ',' in item:
            split_list = item.split(',')
            [merged_annotations_split.append(x) for x in split_list]
        elif ',' not in item:
            split_list = item.split("_")
            original_annotations.append(split_list[0])

    if original_annotations == [] and merged_annotations_split != []:
        for item in merged_annotations_split:
            split_list = item.split("_")
            original_annotations.append(split_list[0])
    else:
        pass
    return original_annotations


def write_bed_file(bed_list, string_ext, base_name):
    """TODO: Docstring for write_bed_file.

    :bed_list: TODO
    :string_ext: TODO
    :base_name: TODO
    :returns: TODO

    """
    output_file_name = base_name + "_" + string_ext + ".bed"
    with open(output_file_name, 'a') as f:
        for item in bed_list:
            f.write('\t'.join(item))
            f.write('\n')


original_annotations = []
updated_annotation_names = []

with open(argv[1], 'r') as f:
    for line in f:
        cleaned_line = line.strip().split("\t")
        original_annotations.append(cleaned_line[0])
        updated_annotation_names.append(cleaned_line[1])

#Clean original annotations
cleaned_original_annotations = grab_original_annotations(original_annotations)

cleaned_original_annotations = read_bed_file_pull_passing(argv[2], cleaned_original_annotations)
#updated_annotation_names = read_bed_file_pull_passing(argv[3], updated_annotation_names)

write_bed_file(cleaned_original_annotations, "original", argv[4])
#write_bed_file(updated_annotation_names, "updated", argv[4])
