library("tidyverse")
library(ggridges)
library(extrafont)


pablo_lncRNA_theme <- function(base_size = 22,
                               base_family = "Arial",
                               base_line_size = base_size / 170,
                               base_rect_size = base_size / 170){
  theme_bw(base_size = base_size, 
           base_family = base_family,
           base_line_size = base_line_size) %+replace%
    theme(
      plot.title = element_text(
        color = "black", 
        hjust = .5, 
        ),
      axis.title = element_text(
        color = "black",
        ),
      panel.border = element_blank(),
      axis.text = element_text(
        color = "black",
        ),
      legend.text = element_text(
        color = "black",
        ),
      panel.grid.major = element_line(colour = "white",linetype = "solid"),   
      panel.grid.minor = element_line(colour = "white",
                                      linetype = "solid", 
                                      size = rel(1)),  
      text = element_text(family = "Arial",, face = "plain",
                          color = "black", size = base_size,
                          hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                          margin = margin(), debug = FALSE),
      
      complete = TRUE
    )
}





H3K36me3_colors <- "#CB2026"
H3K4me1_colors <-  "#642165"
H3K4me3_colors <- "#B89BC9"
H3K56ac_colors <- "#99CA3C"
ATAC_colors <- "#FFCD05"
H2A_colors <-  "#66c2a5"
H3K27me3_colors <- "#234bb8"


color_list <- c(H3K36me3_colors, H3K4me3_colors, H3K27me3_colors, H3K56ac_colors, ATAC_colors, H2A_colors,H3K4me1_colors)


pablo_lncRNA_theme_2 <- function(base_size = 22,
                               base_family = "sans",
                               base_line_size = base_size / 170,
                               base_rect_size = base_size / 170){
  theme_bw(base_size = base_size, 
           base_family = base_family,
           base_line_size = base_line_size) %+replace%
    theme(
      plot.title = element_text(
        color = "black", 
        hjust = .5, 
        size = rel(1.3)),
      axis.title = element_text(
        color = "black",
        size = rel(1)),
      axis.text = element_text(
        color = "black",
        size = rel(1)),
      text = element_text(family = "sans",, face = "plain",
                          color = "black", size = base_size,
                          hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                          margin = margin(), debug = FALSE),
      
      complete = TRUE
    )
}

setwd("/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure_6/00.data/gathered_species_files")

read_bed <- function(updated_bed_anno, type, species){
  bed_cols = c("chrom", "start", "stop", "gene_name", "score", "strand")
  
  updated_bed_group <- read_delim(updated_bed_anno, '\t', col_names = bed_cols, col_types = "cddcccc")
  
  updated_bed <- updated_bed_group %>% 
    mutate(group_ID = as.factor(type)) %>% 
    mutate(gene_length = stop - start) %>% 
    mutate(species_name = species)
  
  
  
  return(updated_bed)
  
}


read_file_set <- function(file_base, species){
  
  
  hyper_large_genes <- str_c(file_base, "_tis_leaf_annotation_hyper_large_genes.bed", sep="")
  major_extension_genes <- str_c(file_base, "_tis_leaf_annotation_major_extension_genes.bed", sep="")
  merged_genes <- str_c(file_base, "_tis_leaf_annotation_merged_genes.bed", sep="")
  minor_extension_genes <- str_c(file_base, "_tis_leaf_annotation_minor_extension_genes.bed", sep="")
  novel_genes <- str_c(file_base, "_tis_leaf_annotation_novel_genes.bed", sep="")
  un_altered_genes <- str_c(file_base, "_tis_leaf_annotation_un_altered_genes.bed", sep="")
  
  

  
  hyper_large_genes_read <- read_bed(hyper_large_genes,"Hyper large", species)
  major_extension_genes_read <- read_bed(major_extension_genes,"Major extension", species)
  merged_genes_read <- read_bed(merged_genes,"Merged genes", species)
  minor_extension_genes_read <- read_bed(minor_extension_genes,"Minor extension", species)
  novel_genes_read <- read_bed(novel_genes,"Novel", species)
  un_altered_genes_read <- read_bed(un_altered_genes,"Unalterd genes", species)
  
  
  final_list <- bind_rows(hyper_large_genes_read,major_extension_genes_read,merged_genes_read,minor_extension_genes_read,novel_genes_read,un_altered_genes_read)
  
  
  
  return(final_list)
  
}
read_file_set_mays <- function(species){
  
  
  hyper_large_genes_read <- read_bed("hyper_large_genes_updated_passing.bed","Hyper large", species)
  major_extension_genes_read <- read_bed("major_extension_genes_updated_passing.bed","Major extension", species)
  merged_genes_read <- read_bed("merged_genes_updated_passing.bed","Merged genes", species)
  minor_extension_genes_read <- read_bed("minor_extension_genes_updated_passing.bed","Minor extension", species)
  novel_genes_read <- read_bed("novel_genes_updated_passing.bed","Novel", species)
  un_altered_genes_read <- read_bed("Zea_mays.AGPv4.38.chr.unaltered.bed", "Unalterd genes", species)
  
  
  final_list <- bind_rows(hyper_large_genes_read,major_extension_genes_read,merged_genes_read,minor_extension_genes_read,novel_genes_read, un_altered_genes_read)
  
  
  
  return(final_list)
  
}


group_annotation_classes <- function(species_level_tribble) {
  
  new_tribble <- species_level_tribble %>% 
    group_by(group_ID, species_name) %>% 
    summarise(counts = n())
  
  return(new_tribble)

  }





asapargus <- read_file_set("Asparagus_10days", "A. officianlis")
phaseolus <- read_file_set("Phaseolus_10days", "P. vulgaris")
setaria <- read_file_set("Setaria_7days", "S. viridis")
sorghum <- read_file_set("Sorghum_7days", "S. bicolor")
soybean <- read_file_set("Soybean_10days", "G. max")
zea_mays <- read_file_set_mays('Z. mays')



asapargus_collapsed <- group_annotation_classes(asapargus)
#brachy_collapsed <- group_annotation_classes(brachy)
phaseolus_collapsed <- group_annotation_classes(phaseolus)
setaria_collapsed <- group_annotation_classes(setaria)
sorghum_collapsed <- group_annotation_classes(sorghum)
soybean_collapsed <- group_annotation_classes(soybean)
zea_mays_collapsed <- group_annotation_classes(zea_mays)



collapsed_all_species <- bind_rows(asapargus_collapsed,
#brachy_collapsed,
phaseolus_collapsed,
setaria_collapsed,
sorghum_collapsed, 
soybean_collapsed, 
zea_mays_collapsed) 

collapsed_all_species_no_unaletered <- collapsed_all_species %>% 
  filter(group_ID != "Unalterd genes") %>% 
  select(species_name, group_ID, counts)



collapsed_all_species_no_unaletered

write.table(collapsed_all_species_no_unaletered, file = "/Users/feilab/Desktop/test.csv", sep = ",", quote = FALSE, row.names = F)

collapsed_all_species_no_unaletered %>% 
  ungroup() %>% 
  group_by(species_name) %>% 
  summarise(sum(counts))

View(collapsed_all_species)

ggplot(collapsed_all_species, aes(fill=group_ID, y=counts, x=species_name)) + 
  geom_bar(position="dodge", stat="identity")

quick_plot <- ggplot(collapsed_all_species_no_unaletered, aes(fill=group_ID, y=counts, x=species_name))  + pablo_lncRNA_theme_2() + 
  geom_bar(position="dodge", stat="identity") + scale_fill_manual("Legend", values = color_list)  + theme(legend.justification = c(1, 1), legend.position = c(.99, .99)) +
  xlab("Species") + ylab("Count") + ggtitle("Counts of annotation classes found in each species") + theme(axis.text.x = element_text(face = "italic"))

#unique(collapsed_all_species_no_unaletered$species_name)

quick_plot

ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure_6/imgs", filename="multi_species_counts.svg", plot=quick_plot, width = 4, height = 4, units = "in")





gene_count_file_name <- "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure_6/00.data/number_genes_species.csv"
gene_count_read <- read_csv(gene_count_file_name, col_names = TRUE)
gene_count_read





species_number_annotation_errors_gene_number <- collapsed_all_species_no_unaletered %>% 
  group_by(species_name) %>% 
  summarise(number_annot_errors = sum(counts)) %>% 
  left_join(.,  gene_count_read, by = "species_name")

gene_number_vs_annotation_error <- species_number_annotation_errors_gene_number %>% 
  ggplot(., aes(y = number_annot_errors, x = number_genes, label=species_name)) + geom_point() + pablo_lncRNA_theme() + 
  geom_text(aes(label=species_name),vjust="inward",hjust="inward") + 
  xlab("Number of Genes in Genome") + ylab("Counts")



genome_size_vs_annotation_error_all_counts <- species_number_annotation_errors_gene_number %>% 
  ggplot(., aes(y = number_annot_errors, x = (genome_size_bp/1000000000), label=species_name)) + geom_point() + pablo_lncRNA_theme() + 
  geom_text(aes(label=species_name),vjust="inward",hjust="inward") + 
  xlab("Genome Size (Gb)") + ylab("Counts") 


genome_size_vs_annotation_error_all_counts


species_number_annotation_errors_gene_number_2 <- collapsed_all_species_no_unaletered %>% 
  left_join(.,  gene_count_read, by = "species_name")

gene_number_vs_annotation_error_by_type <- species_number_annotation_errors_gene_number_2 %>% 
  ggplot(., aes(y = counts, x = number_genes, label=species_name)) + geom_point() + pablo_lncRNA_theme() + 
  geom_text(aes(label=species_name),vjust="inward",hjust="inward") + facet_grid(group_ID~.,scales = "free") +
  xlab("Number of Genes in Genome") + ylab("Counts")

gene_number_vs_annotation_error_by_type

library(ggrepel)
genome_size_vs_annotation_error_by_type <- species_number_annotation_errors_gene_number_2 %>% 
  ggplot(., aes(y = counts, x = (genome_size_bp/1000000000), label=species_name)) + geom_point(size = 3) + pablo_lncRNA_theme_2() +
  #geom_text(aes(label=species_name, fontface = "italic"),vjust="inward",hjust="inward") + 
  facet_grid(group_ID~.,  scales = "free") +
  geom_text_repel(aes(label=species_name, fontface = "italic"), size = rel(22 * .3), box.padding = .45) + 
  xlab("Genome Size (Gb)") + ylab("Counts") + ggtitle("Genome size versus annotation counts") 
  

genome_size_vs_annotation_error_by_type


ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure_6/imgs", filename="gene_number_vs_annotation_error.pdf", plot=gene_number_vs_annotation_error, width = 8, height = 8, units = "in")
ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure_6/imgs", filename="genome_size_vs_annotation_error_all_counts.pdf", plot=genome_size_vs_annotation_error_all_counts, width = 10, height = 8, units = "in")
ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure_6/imgs", filename="gene_number_vs_annotation_error_by_type.pdf", plot=gene_number_vs_annotation_error_by_type, width = 10, height = 12, units = "in")
ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure_6/imgs", filename="genome_size_vs_annotation_error_by_type.svg", plot=genome_size_vs_annotation_error_by_type, width = 4, height = 4, units = "in")





combine_all_bed_files <- bind_rows(asapargus,phaseolus,setaria,sorghum,soybean,zea_mays)


gene_density_chart <- combine_all_bed_files %>% 
  ggplot(., aes(x=log10(gene_length), y=group_ID, fill=group_ID)) +  
  geom_density_ridges() +
  labs(title = "Length of Each Gene Class Found in Each Species", x = "Log2(Gene Length)", y = "test") +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  #theme(legend.justification = c(1, 1), legend.position = c(.99, .99)) + 
  theme(axis.title.y =element_blank()) + 
  facet_grid(species_name~.) + 
  scale_fill_manual("legend", values = color_list)

gene_density_chart

ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure_6/imgs", filename="gene_density_chart_by_species.pdf", plot=gene_density_chart, width = 10, height = 12, units = "in")


combine_all_bed_files_2 <- combine_all_bed_files %>% 
  mutate(annotation_class_broad = case_when(group_ID  == "Hyper large" ~ "Updated annotations",
                                            group_ID == "Major extension" ~ "Updated annotations",
                                            group_ID == "Merged genes" ~ "Updated annotations",
                                            group_ID == "Minor extension" ~ "Updated annotations",
                                            group_ID == "Novel" ~ "Novel",
                                            group_ID == "Unalterd genes" ~ "Unalterd annotations"
                                            ))


gene_density_chart_2 <- combine_all_bed_files_2 %>% 
  ggplot(., aes(x=log10(gene_length), y=annotation_class_broad, fill=annotation_class_broad)) +  
  geom_density_ridges() +
  labs(title = "Length of Each Gene Class Found in Each Species", x = "Log2(Gene Length)", y = "test") +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  #theme(legend.justification = c(1, 1), legend.position = c(.99, .99)) + 
  theme(axis.title.y =element_blank()) + 
  facet_grid(species_name~.) + 
  scale_fill_manual("legend", values = color_list)

gene_density_chart_2

ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure_6/imgs", filename="gene_density_chart_by_species_broaded_types.pdf", plot=gene_density_chart_2, width = 8, height = 12, units = "in")











combine_all_bed_files %>% 
  ggplot(., aes(x=(gene_length), y=group_ID, fill=group_ID)) +  
  geom_density_ridges() +
  labs(title = "Novel Class Length", x = "Log2(Gene Length)", y = "test") +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  theme(legend.justification = c(1, 1), legend.position = c(.99, .99)) + 
  theme(axis.title.y =element_blank()) + 
  facet_grid(species_name~.) + 
  scale_fill_manual("legend", values = color_list)


combine_all_bed_files %>% 
  filter(group_ID == "novel_genes") %>% 
  ggplot(., aes(x=(gene_length), y=group_ID, fill=group_ID)) +  
  geom_density_ridges() +
  labs(title = "Novel Class Length", x = "Log2(Gene Length)", y = "test") +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  theme(legend.justification = c(1, 1), legend.position = c(.99, .99)) + 
  theme(axis.title.y =element_blank()) + 
  facet_grid(species_name~.) + 
  scale_fill_manual("legend", values = color_list)

combine_all_bed_files %>% 
  filter(group_ID == "hyper_large_genes" |group_ID == "major_extension_genes" | group_ID == "minor_extension_genes") %>% 
  ggplot(., aes(x=log(gene_length), y=group_ID, fill=group_ID)) +  
  geom_density_ridges() +
  labs(title = "Novel Class Length", x = "Log2(Gene Length)", y = "test") +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  theme(legend.justification = c(1, 1), legend.position = c(.99, .99)) + 
  theme(axis.title.y =element_blank()) + 
  facet_grid(species_name~.) + 
  scale_fill_manual("legend", values = color_list)

combine_all_bed_files %>% 
  filter(group_ID == "merged_genes") %>% 
  ggplot(., aes(x=log(gene_length), y=group_ID, fill=group_ID)) +  
  geom_density_ridges() +
  labs(title = "Novel Class Length", x = "Log2(Gene Length)", y = "test") +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  theme(legend.justification = c(1, 1), legend.position = c(.99, .99)) + 
  theme(axis.title.y =element_blank()) + 
  facet_grid(species_name~.) + 
  scale_fill_manual("legend", values = color_list)

