library(dplyr)
library(tidyverse)

setwd("/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure3_5/04.merged_gene_class_analysis")


load_merged_genes <- readr::read_delim("merged_genes_overlap_aware_non_overlapping_list.txt", delim=",", col_names=c("gene1", "gene2", "gene3")) %>% 
  mutate(positionInCategory = 1:n()) %>% 
  mutate(gene_cluster = str_c("merge_group", positionInCategory, sep = '_')) %>% 
  select(gene1, gene2, gene3, gene_cluster) %>% 
  pivot_longer(c("gene1", "gene2", "gene3"), names_to = "gene_count", values_to = "gene", names_repair="minimal")


load_merged_genes %>% 
  select(gene_count) %>% 
  count(gene_count)




load_gene_annotations <- readr::read_delim("Zea_mays.AGPv4.38.gene_name.description.txt", delim="\t", col_names=c("gene", "description"))

cleaned_description_names <- load_gene_annotations %>% 
  separate(description, c("first", "second", "third"), sep = ";", extra = "drop") %>% 
  select(gene, first) %>% 
  filter(stringr::str_detect(first, 'description')) %>% 
  mutate(first = gsub("description=", "", first))
  

  
gather_all_genes <- left_join(load_merged_genes, cleaned_description_names, by = c("gene"))

length(unique(gather_all_genes$gene_cluster))

number_items_in_groups_after_drop <- gather_all_genes %>%
  group_by(gene_cluster) %>% 
  count()



number_items_in_groups_after_drop
add_group_counts <- left_join(gather_all_genes, number_items_in_groups_after_drop, by="gene_cluster") %>% 
  filter(n > 1) %>% 
  drop_na()

add_group_counts %>% 
  group_by(gene_cluster) %>% 
  select(first, n) %>% 
  summarise(number_go_terms = n_distinct(first)) %>% 
  count(number_go_terms) %>% 
  mutate(freq = n / length(unique(load_merged_genes$gene_cluster)))


length(unique(add_group_counts$gene_cluster))




# Compare Against Blast ---------------------------------------------------
blast_col_names <-  c("gene_1", "gene_2", "pident" ,"length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
read_in_blast <- readr::read_delim("mergers_pre_assembly_pep_all_vs_all.blast", delim = "\t", col_names = blast_col_names)

mean_blast_scores <- read_in_blast %>% 
  group_by(gene_1, gene_2) %>% 
  top_n(1, wt = "length") %>% 
  select(gene_1, gene_2, pident ,length, mismatch, evalue, bitscore)
  



            
previous_dup_names_cols  <- c("gene_cluster","gene_1", "gene_2", "gene_3", "gene_4")
read_in_known_tandem_dups <- readr::read_csv("previously_known_tandem_duplicates.csv", col_names = previous_dup_names_cols) %>% 
  select("gene_cluster","gene_1", "gene_2") %>% 
  mutate(grouping_name = "Known tadnem duplicates")


prev_known_tandem_blast_p <- left_join(read_in_known_tandem_dups, mean_blast_scores, by = c("gene_1", "gene_2")) %>% 
  arrange(gene_cluster,gene_1,gene_2, grouping_name, pident ,length, mismatch, evalue, bitscore)



#Load Pre Assembly Gene class - add numbers 
load_merged_genes_replaced_names <- readr::read_csv("merged_genes_overlap_aware_non_overlapping_list.txt", col_names=c("gene1", "gene2", "gene3")) %>% 
  mutate(positionInCategory = 1:n()) %>% 
  mutate(gene_cluster = str_c("merge_group", positionInCategory, sep = '_')) %>% 
  select(gene1, gene2, gene_cluster) %>% 
  mutate(grouping_name = "Merged genes") %>% 
  mutate(gene_1 = gene1) %>% 
  mutate(gene_2 = gene2) %>% 
  select(-gene1, -gene2)


merged_gene_blast_p <- left_join(load_merged_genes_replaced_names, mean_blast_scores, by = c("gene_1", "gene_2")) %>% 
  arrange(gene_cluster,gene_1,gene_2, grouping_name, pident ,length, mismatch, evalue, bitscore)


add_group_counts

rw_widen_merged_genes <- add_group_counts %>% 
  select(-n) %>% 
  pivot_wider(names_from = c(gene_count), values_from = c(gene, first))


widen_GO_terms_identified <- rw_widen_merged_genes %>% 
  select(-gene_gene3, -first_gene3) %>% 
  mutate(yes_no_identical_GO_term = case_when( first_gene1 == first_gene2 ~ 1, 
                                               first_gene1 != first_gene2 ~ 0)) %>% 
  mutate(gene_1 = gene_gene1) %>% 
  mutate(gene_2 = gene_gene2) 


View(widen_GO_terms_identified)
merged_genes_with_blastp <- left_join(widen_GO_terms_identified, mean_blast_scores, by = c("gene_1", "gene_2"))

View(merged_genes_with_blastp)
merged_genes_with_blastp %>% 
  mutate(BLAST_P_HIT = case_when(is.na(pident) == TRUE ~ 0, 
                                 is.na(pident) == FALSE & pident < 50  ~ 0, 
                                 is.na(pident) == FALSE & pident > 50  ~ 1, )) %>% 
  select(gene_cluster, gene_1, gene_2,  yes_no_identical_GO_term, BLAST_P_HIT, first_gene1, first_gene2) %>% 
  unique() %>% 
  group_by(yes_no_identical_GO_term) %>% 
  count(BLAST_P_HIT) 


merged_genes_with_blastp %>% 
  mutate(BLAST_P_HIT = case_when(is.na(pident) == TRUE ~ 0, 
                                 is.na(pident) == FALSE & pident < 50  ~ 0, 
                                 is.na(pident) == FALSE & pident > 50  ~ 1, )) %>% 
  filter(BLAST_P_HIT == 1) %>% 
  group_by(gene_cluster) %>% 
  unique() %>% 
  View()



merged_genes_with_blastp %>% 
  mutate(BLAST_P_HIT = case_when(is.na(pident) == TRUE ~ 0, 
                                 is.na(pident) == FALSE ~ 1, )) %>% 
  filter(yes_no_identical_GO_term == 1, BLAST_P_HIT == 1) %>% 
  group_by(gene_cluster) %>% 
  unique() %>% 
  View()
100 - 1.4

128 + 62 + 87

277 + 4

277/281

bed_col_names <-  c("crhom", "start", "stop" ,"name", "none", "strand", "other", "other_1")
quick_med_gene_len <- readr::read_delim("/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure3_5/00.data/bed_files/Zea_mays.AGPv4.38.chr.all.bed", delim = "\t", col_names = bed_col_names)

quick_med_gene_len %>% 
  mutate(total_len = stop - start) %>% 
  group_by(other_1) %>% 
  summarise(gene_len_med = median(total_len))
                                 