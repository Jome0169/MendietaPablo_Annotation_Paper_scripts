from sys import argv

#List of "Novel" regions found
good_list = []
with open(argv[1], 'r') as f:
    for line in f:
        clean_line = line.strip()
        good_list.append(clean_line)


good_lines = []
with open(argv[2], 'r') as z:
    for line in z:
        if line.startswith("#"):
            pass
        else:
            clean_line = line.strip().split('\t')
            clean_line_2 = line.strip()
            gene_ID_line = clean_line[8].split(';')[0].split(' ')[1].replace("\"", '')
            if gene_ID_line in good_list:
                good_lines.append(clean_line_2)
            else:
                pass

for item in good_lines:
    print(item)
