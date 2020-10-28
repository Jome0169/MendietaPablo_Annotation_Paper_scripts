from sys import argv

read_in_bed = []

with open(argv[1], 'r') as f:
    for line in f:
        cleaned_line = line.strip().split()
        read_in_bed.append(cleaned_line)



def replace_appended_values(string_item):
    """TODO: Docstring for replace_appended_values.

    :string_item: TODO
    :returns: TODO

    """
    if "major_extend" in string_item:
        split_gene_name = string_item.split("_")[0]
    elif "minor_extend" in string_item:
        split_gene_name = string_item.split("_")[0]
    elif "hyper_large" in string_item:
        split_gene_name = string_item.split("_")[0]
    else:
        split_gene_name = string_item
    return split_gene_name


final_correct_list = []
for item in read_in_bed:
    take_name = item[3]
    split_on_sep = take_name.split("|")
    split_on_comma = split_on_sep[0].split(',')
    correct_gene_name = []
    for gene_name in split_on_comma:
        correct_gene = replace_appended_values(gene_name)
        correct_gene_name.append(correct_gene)
    final_correct_list.append(list(set(correct_gene_name)))


for item_pairs in final_correct_list:
    print(','.join(item_pairs))
