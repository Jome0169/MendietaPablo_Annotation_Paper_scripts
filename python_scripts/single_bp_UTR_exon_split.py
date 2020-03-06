from sys import argv
import copy




def process_GTF(read_gtf_count_file):
    """reads in gtf count file
    :returns: TODO

    """
    def grab_gene_name(ID_string):
        """Pulls out the acutal gene name for each region
        :returns: TODO
        """
        split_string = ID_string.split(';')
        take_parent = split_string[0]


        if "ID=" in take_parent and "ID=CDS" not in take_parent: 
            #ID=transcript:Zm00001e000001_T001;Parent=gene:Zm00001e000001;
            cleaned_gene_name = take_parent.replace("ID=transcript:","")
        elif "ID=CDS" in take_parent:
            #ID=CDS:Zm00001e000001_P001;Parent=transcript:Zm00001e000001_T001;
            take_parent2 = split_string[1]
            cleaned_gene_name = take_parent2.replace("Parent=transcript:","")
        elif "Parent=transcript" in take_parent:
            cleaned_gene_name = take_parent.replace("Parent=transcript:","")
        else:
            pass
        return(cleaned_gene_name)


    gtf_file_dict = {}
    nested_dict = {"CDS":[], "five_prime_UTR":[], "three_prime_UTR":[]}
    with open(read_gtf_count_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                pass
            else:
                clean_line = line.strip().split()
                if clean_line[2] == 'mRNA':
                    gene_name = grab_gene_name(clean_line[8])
                    gtf_file_dict[gene_name] = copy.deepcopy(nested_dict)
                elif clean_line[2] == 'CDS':
                    gene_name = grab_gene_name(clean_line[8])
                    gtf_file_dict[gene_name]["CDS"].append(clean_line)
                elif clean_line[2] == 'five_prime_UTR':
                    gene_name = grab_gene_name(clean_line[8])
                    gtf_file_dict[gene_name]["five_prime_UTR"].append(clean_line)
                elif clean_line[2] == 'three_prime_UTR':
                    gene_name = grab_gene_name(clean_line[8])
                    gtf_file_dict[gene_name]["three_prime_UTR"].append(clean_line)
                else:
                    pass
    return(gtf_file_dict)


def split_genomic_regions(gtf_file_read):
    """Counts the number of reads in each gene
    :gtf_file_read: TODO
    :returns: TODO

    """
    def generate_single_bp_regions(chrom, start, stop, strand, genome_type, transcript_name):
        """TODO: Docstring for generate_single_bp_regions(.

        """
        nested_bed = []

        counter = 1 
        for bp_val in range(int(start), int(stop), 1):
            generate_name = transcript_name + "__" + genome_type + "__" + "bp" + "__" + str(counter)
            micro_bed = [chrom, str(bp_val), str(bp_val + 1), generate_name,
                    ".", strand]
            nested_bed.append(micro_bed)
            counter += 1
        return nested_bed


    def split_region(exon_dict, transcript_name, genome_type):
        """
        Takes in the exon nested dict, isolates the first exon and counts the
        numer of reads found in it, as well as the nth exons, and the number of
        reads falling into thos regions.
        """
        split_region = []

        if genome_type == "CDS":
            strand = exon_dict[0][6]
            if strand == "+":
                chromosome = exon_dict[0][0]
                exon_start = exon_dict[0][3]
                exon_stop = exon_dict[0][4]
                strand = exon_dict[0][6]

                nested_beds = generate_single_bp_regions(chromosome,
                        exon_start, exon_stop, strand, genome_type,
                        transcript_name)
                split_region.append(nested_beds)
            elif strand == "-":
                chromosome = exon_dict[-1][0]
                exon_start = exon_dict[-1][3]
                exon_stop = exon_dict[-1][4]
                strand = exon_dict[-1][6]
                nested_beds = generate_single_bp_regions(chromosome,
                        exon_start, exon_stop, strand, genome_type,
                        transcript_name)
                split_region.append(nested_beds)

                
        elif genome_type == "fp_UTR":
            for item in exon_dict:
                    chromosome = item[0]
                    exon_start = item[3]
                    exon_stop = item[4]
                    strand = item[6]

                    nested_beds = generate_single_bp_regions(chromosome,
                            exon_start, exon_stop, strand, genome_type,
                            transcript_name)
                    split_region.append(nested_beds)
        else:
            pass
        return split_region



    for transcript, stuff in gtf_file_read.items():
       split_exons = split_region(stuff['CDS'], transcript, "CDS")
       split_5pUTR = split_region(stuff['five_prime_UTR'], transcript, "fp_UTR")

       combined_nested_list = split_exons + split_5pUTR

       for nested_list in combined_nested_list:
           for small_list in nested_list:
               print('\t'.join(small_list))

    


cleaned_gtf = process_GTF(argv[1])
split_genomic_regions(cleaned_gtf)
