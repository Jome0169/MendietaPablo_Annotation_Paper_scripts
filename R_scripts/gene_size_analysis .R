library("tidyverse")
library(ggridges)
library(refGenome)


H3K36me3_colors <- "#CB2026"
H3K4me1_colors <-  "#642165"
H3K4me3_colors <- "#B89BC9"
H3K56ac_colors <- "#99CA3C"
ATAC_colors <- "#FFCD05"
H2A_colors <-  "#66c2a5"
H3K27me3_colors <- "#234bb8"


color_list <- c(H3K36me3_colors, H3K4me3_colors, H3K27me3_colors, H3K56ac_colors, ATAC_colors, H2A_colors,H3K4me1_colors)


setwd("/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure_7/")

read_bed <- function(updated_bed_anno, type){
  bed_cols = c("chrom", "start", "stop", "gene_name", "score", "strand")
  
  updated_bed_group <- read_delim(updated_bed_anno, '\t', col_names = bed_cols, col_types = "cddcccc")
  
  updated_bed <- updated_bed_group %>% 
    mutate(group_ID = as.factor(type)) %>% 
    mutate(gene_length = stop - start) %>% 
    mutate(group_name = "Original Annotation")
  
  
  
  return(updated_bed)
  
}



read_before_after_bed_take_mean <- function(ISO_seq_bed, WGS_bed, type){
  bed_cols = c("chrom", "start", "stop", "gene_name", "score", "strand")

  
  iso_seq_bed <- read_delim(ISO_seq_bed, '\t', col_names = bed_cols, col_types = "dddcccc")
  iso_updated_bed <- iso_seq_bed %>% 
    mutate(group_ID = as.factor(type)) %>% 
    mutate(group_name = "Updated Annotation") %>% 
    mutate(gene_length = stop - start) %>% 
    mutate(data_type = "Iso_seq") %>% 
    select(gene_name, group_ID, group_name, gene_length)
  
  wgs_seq_bed <- read_delim(WGS_bed, '\t', col_names = bed_cols, col_types = "dddcccc")
  wgs_updated_bed <- wgs_seq_bed %>% 
    mutate(group_ID = as.factor(type)) %>% 
    mutate(group_name = "Updated Annotation") %>% 
    mutate(gene_length = stop - start) %>% 
    mutate(data_type = "WGS") %>% 
    select(gene_name, group_ID, group_name, gene_length)
  
  
  merge_iso_wgs_assemblies <-  bind_rows(wgs_updated_bed, iso_updated_bed) %>% 
    group_by(gene_name, group_ID, group_name) %>% 
    summarise(new_len = mean(gene_length)) %>% 
    rename(gene_length = new_len)
  
  
  
  combined_bed_file <- merge_iso_wgs_assemblies
  return(combined_bed_file)
  
}



Zea_mays_un_modified_genes <- read_bed("00.data/Annotations_passing/Zea_mays.AGPv4.38.un_modified_genes.bed", "un_modified_genes")
hyper_large_updated <- read_before_after_bed_take_mean("00.data/Annotations_passing/06.iso_seq_passing_assemblies/hyper_large_genes_updated_passing.bed", "00.data/Annotations_passing/07.WGS_passing_assemblies/hyper_large_genes_updated_passing.bed", "hyper_large")
major_extension_updated <- read_before_after_bed_take_mean("00.data/Annotations_passing/06.iso_seq_passing_assemblies/major_extension_genes_updated_passing.bed", "00.data/Annotations_passing/07.WGS_passing_assemblies/major_extension_genes_updated_passing.bed", "major_extension")
merged_updated <- read_before_after_bed_take_mean("00.data/Annotations_passing/06.iso_seq_passing_assemblies/merged_genes_updated_passing.bed", "00.data/Annotations_passing/07.WGS_passing_assemblies/merged_genes_updated_passing.bed", "merged")
minor_extension_updated <- read_before_after_bed_take_mean("00.data/Annotations_passing/06.iso_seq_passing_assemblies/minor_extension_genes_updated_passing.bed", "00.data/Annotations_passing/07.WGS_passing_assemblies/minor_extension_genes_updated_passing.bed", "minor_extension")
novel_updated <- read_before_after_bed_take_mean("00.data/Annotations_passing/06.iso_seq_passing_assemblies/novel_genes_updated_passing.bed", "00.data/Annotations_passing/07.WGS_passing_assemblies/novel_genes_updated_passing.bed", "Novel")






all_genes <- bind_rows(Zea_mays_un_modified_genes,
hyper_large_updated,
major_extension_updated,
merged_updated,
minor_extension_updated,novel_updated)


bed_gene_lens <- all_genes %>% 
  ggplot(., aes(x = group_ID , y = gene_length)) + 
  geom_violin(trim=TRUE) + geom_jitter(shape=2, position=position_jitter(0.2), size = .1) 


bed_gene_lens

ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure_7/imgs", filename="disiribution_of_annotation_lengths.png", plot=bed_gene_lens, width = 8, height = 8, units = "in")


top_1_percent <- all_genes %>% 
  top_frac(., .01, gene_length) %>% 
  mutate(base = 'none') %>% 
  mutate(base_group = "top_1_percent_genes")



top_1_percent_graph <- top_1_percent %>% 
  ggplot(.,aes(x = base, fill = group_ID)) + 
  geom_bar(position = "fill")

top_1_percent_graph
ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure_7/imgs", filename="proportion_of_top_1_percent.png", plot=top_1_percent_graph, width = 8, height = 8, units = "in")


count_number_rows <- nrow(top_1_percent)

set.seed(123)
equal_sample <- sample_n(all_genes, count_number_rows, replace = FALSE) %>% 
  mutate(base_group = "random_sample_normal_genes")

bound_top_1_percetn_random_sample <- bind_rows(top_1_percent,equal_sample)

bound_top_1_percetn_random_sample %>% 
  ggplot(., aes(x=gene_length, y=base_group)) +  
  geom_density_ridges(scale = .8, alpha = .8)  +
  labs(title = "Novel Class Length", x = "Log2(Gene Length)", y = "test") +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() + pablo_lncRNA_theme() +
  theme(legend.justification = c(1, 1), legend.position = c(.99, .99)) + 
  theme(axis.title.y =element_blank()) 
  


top_001_percent_genes <- top_1_percent %>% 
  select(chrom, start,stop,gene_name,score,strand)



write.table(top_001_percent_genes, file="01.pull_GTF_files/top_001_percent_genes.bed", sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)  
