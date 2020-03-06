"""
Purpose - takes in a GTF like bed file with number of counts associated with
it.

For example, this file takes in something similar to:

chr1    NAM     gene    34617   40204   .       +       .       ID=gene:Zm00001e000001;biotype=protein_coding;logic_name=mikado_gene    10
chr1    NAM     mRNA    34617   40204   .       +       .       ID=transcript:Zm00001e000001_T001;Parent=gene:Zm00001e000001;biotype=protein_coding;transcript_id=Zm00001e000001_T001;canonical_transcript=1;_AED=0.23;_QI=0|0.75|0.66|0.77|0|0|9|1668|648;Dbxref=InterPro:IPR003690,Pfam:PF02536;Ontology_term=GO
chr1    NAM     exon    34617   35318   .       +       .       Parent=transcript:Zm00001e000001_T001;Name=Zm00001e000001_T001.exon.1;ensembl_end_phase=0;ensembl_phase=0;exon_id=Zm00001e000001_T001.exon.1;rank=1     4
chr1    NAM     CDS     34617   35318   .       +       0       ID=CDS:Zm00001e000001_P001;Parent=transcript:Zm00001e000001_T001;protein_id=Zm00001e000001_P001 4
chr1    NAM     exon    36037   36174   .       +       .       Parent=transcript:Zm00001e000001_T001;Name=Zm00001e000001_T001.exon.2;ensembl_end_phase=0;ensembl_phase=0;exon_id=Zm00001e000001_T001.exon.2;rank=2     0
chr1    NAM     CDS     36037   36174   .       +       0       ID=CDS:Zm00001e000001_P001;Parent=transcript:Zm00001e000001_T001;protein_id=Zm00001e000001_P001 0
chr1    NAM     exon    36259   36504   .       +       .       Parent=transcript:Zm00001e000001_T001;Name=Zm00001e000001_T001.exon.3;ensembl_end_phase=0;ensembl_phase=0;exon_id=Zm00001e000001_T001.exon.3;rank=3     0
chr1    NAM     CDS     36259   36504   .       +       0       ID=CDS:Zm00001e000001_P001;Parent=transcript:Zm00001e000001_T001;protein_id=Zm00001e000001_P001 0
chr1    NAM     exon    36600   36713   .       +       .       Parent=transcript:Zm00001e000001_T001;Name=Zm00001e000001_T001.exon.4;ensembl_end_phase=0;ensembl_phase=0;exon_id=Zm00001e000001_T001.exon.4;rank=4     0
chr1    NAM     CDS     36600   36713   .       +       0       ID=CDS:Zm00001e000001_P001;Parent=transcript:Zm00001e000001_T001;protein_id=Zm00001e000001_P001 0
chr1    NAM     exon    36822   37004   .       +       .       Parent=transcript:Zm00001e000001_T001;Name=Zm00001e000001_T001.exon.5;ensembl_end_phase=0;ensembl_phase=0;exon_id=Zm00001e000001_T001.exon.5;rank=5     0
chr1    NAM     CDS     36822   37004   .       +       0       ID=CDS:Zm00001e000001_P001;Parent=transcript:Zm00001e000001_T001;protein_id=Zm00001e000001_P001 0
chr1    NAM     exon    37416   37633   .       +       .       Parent=transcript:Zm00001e000001_T001;Name=Zm00001e000001_T001.exon.6;ensembl_end_phase=2;ensembl_phase=2;exon_id=Zm00001e000001_T001.exon.6;rank=6     0
chr1    NAM     CDS     37416   37633   .       +       0       ID=CDS:Zm00001e000001_P001;Parent=transcript:Zm00001e000001_T001;protein_id=Zm00001e000001_P001 0

Important to note that this script requries that the ID attributes on the end
be formatted in a very specific way. This script also counts relative ISO-FORM
COUNTS. NOT GENE COUNTS. This is becasue R2C2 start sites are more than likely
iso-form specific.

RUNNING EXAMPLE: python R2C2_generate_counts.py > OUTPUT
"""


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


def count_number_of_reads(gtf_file_read):
    """Counts the number of reads in each gene
    :gtf_file_read: TODO
    :returns: TODO

    """
    def classify_first_exon(exon_list):
        """
        Takes in the exon nested dict, isolates the first exon and counts the
        numer of reads found in it, as well as the nth exons, and the number of
        reads falling into thos regions.
        
        """
        counts_first_exon = None
        length_first_exon = None
        counts_other_exon = []
        length_other_exons = 0
        other_exon_n = 0
        exon_counter = 0
        take_all_starts = [i[4] for i in exon_list]
        strand = list(exon_list)[0][6]
        
        if strand == '+':
            for item in exon_list:
                if exon_counter == 0:
                    counts_first_exon = int(item[-1])
                    length_first_exon = int(item[4]) - int(item[3])
                    exon_counter += 1 
                elif exon_counter > 11:
                    exon_other_read_count = int(item[-1])
                    counts_other_exon.append(exon_other_read_count)
                    length_other_exons += (int(item[4]) - int(item[3]))
                    #Record number other exons
                    other_exon_n += 1 
                    exon_counter += 1 
                else:
                    pass
        elif strand == '-':
            for item in exon_list[::-1]:
                if exon_counter == 0:
                    counts_first_exon = int(item[-1])
                    length_first_exon = int(item[4]) - int(item[3])
                    exon_counter += 1 
                elif exon_counter > 11:
                    exon_other_read_count = int(item[-1])
                    counts_other_exon.append(exon_other_read_count)
                    length_other_exons += (int(item[4]) - int(item[3]))
                    #Record number other exons
                    other_exon_n += 1 
                    exon_counter += 1 
                else:
                    pass


        final_exon_counts = [counts_first_exon, length_first_exon,
                sum(counts_other_exon), length_other_exons, other_exon_n]
        return(final_exon_counts)

    def classify_other(other_dict):
        """TODO: Docstring for classify_other.
        :returns: TODO

        """
        counts_other = 0
        count_n = 0
        length_other_region = 0
        if len(other_dict) != 0:
            for item in other_dict:
                counts_other += int(item[-1])
                length_other_region += (int(item[4]) - int(item[3]))
                count_n += 1
        else:
            pass
        return [counts_other, length_other_region, count_n]



    final_dict = {}
    header = ['transcript_name', 'first_exon_read_count', "length_first_exon", "read_count_other_exons",
            "length_other_exons", "number_other_exons",
            "three_prime_UTR_read_count", "length_three_prime",
            "number_three_prime_UTRs", "five_prime_UTR_read_count",
            "length_five_prime_UTR", 'number_five_prime_UTR']



    print('\t'.join(header))
    for transcript, stuff in gtf_file_read.items():
       exon_list_count = classify_first_exon(stuff['CDS'])
       three_prime_UTR_list = classify_other(stuff['three_prime_UTR'])
       five_prime_UTR_list = classify_other(stuff['five_prime_UTR'])

       final_counts = [transcript] + exon_list_count + three_prime_UTR_list + five_prime_UTR_list
       print('\t'.join([str(i) for i in final_counts]))




cleaned_gtf = process_GTF(argv[1])
count_number_of_reads(cleaned_gtf)





#print(cleaned_gtf["gene:Zm00001e000047"])
