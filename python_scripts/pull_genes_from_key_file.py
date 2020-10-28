from sys import argv




passing_genes = []

with open(argv[1], 'r') as f:
    for line in f:
        passing_genes.append(line.strip())



read_in_key_file = []
with open(argv[2], 'r') as z:
    for line in z:
        cleaned_line = line.strip().split()
        if cleaned_line[1] in passing_genes:
            read_in_key_file.append(cleaned_line)






if len(read_in_key_file) == len(passing_genes):
    for item in read_in_key_file:
        print('\t'.join(item))

elif len(read_in_key_file) != len(passing_genes):
    print("INCORRECT NUMBER OF KEY FILES RETRIEVED")
    print(len(read_in_key_file), len(passing_genes))
    sys.exit(-1)
