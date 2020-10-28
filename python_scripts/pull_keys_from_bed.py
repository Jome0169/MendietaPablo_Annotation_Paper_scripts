"""
Takes in BED file - pulls out passing KEY file. Inverse of pull_list_bed.py

Bed file looking like:
❯ head 00.data/counts_each_tis_bed/leaf_annotation_hyper_large_genes.bed
1       138377  155486  gene:Zm00001d027236_major_extended_hyper_large  2       -
1       2023484 2034168 gene:Zm00001d027293_major_extended_hyper_large  8       -
1       4723400 4736938 gene:Zm00001d027412_major_extended_hyper_large  8       -
1       5790072 5812799 gene:Zm00001d027460_major_extended_hyper_large  16      +
1       6239000 6246787 gene:Zm00001d027478_major_extended_hyper_large  8       -


❯ head 02.generate_distribution_charts/hyper_large_genes.passed.txt
gene:Zm00001d027235_major_extended_hyper_large  hyper_large_genes_1
gene:Zm00001d027236_major_extended_hyper_large  hyper_large_genes_2
gene:Zm00001d027276_major_extended_hyper_large|gene:Zm00001d027276_minor_extended_major_extended_hyper_large    hyper_large_genes_3
gene:Zm00001d027293_major_extended_hyper_large|gene:Zm00001d027293_minor_extended_major_extended_hyper_large    hyper_large_genes_4
gene:Zm00001d027412_major_extended_hyper_large  hyper_large_genes_5
gene:Zm00001d027460_major_extended_hyper_large  hyper_large_genes_7


prints out bed file locations found in the key file. Skips all others.
Example: python pull_list_bed.py
02.generate_distribution_charts/hyper_large_genes.passed.txt
00.data/counts_each_tis_bed/leaf_annotation_hyper_large_genes.bed
"""




from sys import argv


def read_bed_file_pull_passing(arg1):
    """TODO: Docstring for read_bed_file_pull_passing.

    :arg1: TODO
    :returns: TODO

    """
    passing_lines = []
    with open(arg1, 'r') as f:
        for line in f:
            clean_line = line.strip().split()
            passing_lines.append(clean_line[3])
    return(passing_lines)


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




updated_annotation_names = read_bed_file_pull_passing(argv[1])

final_passing_group = []
with open(argv[2], 'r') as f:
    for line in f:
        cleaned_line = line.strip().split()
        take_names = cleaned_line[0].split("|")
        for item in take_names:
            if item in updated_annotation_names:
                final_passing_group.append(cleaned_line)



for i in final_passing_group:
    print("\t".join(i))
