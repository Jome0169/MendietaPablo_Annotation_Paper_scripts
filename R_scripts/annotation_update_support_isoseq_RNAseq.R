#################
#Purpose of this script is to generate the bar charts of the number of updated assembleis that are supported by either the
#iso-seq data, or the RNA-seq data. This script is currently for V4 of the maize genome annotaiton. 


#Analysis on Server Location:
#/scratch/jpm73279/04.lncRNA/02.Analysis/34.generate_updated_annotations

#Input File Location:
#/Users/feilab/Projects/03.ncRNA_project/02.Analysis/lncRNA_copy_files/2020-03-10_ven_diagram_correct_assemblies
#################


library(tidyverse)
library(ggforce)
library(cowplot)

#args <- commandArgs(trailingOnly = TRUE)
#cat(args, sep = "\n")


setwd("/Users/feilab/Projects/03.ncRNA_project/02.Analysis/lncRNA_copy_files/2020-03-10_ven_diagram_correct_assemblies")


H3K36me3_colors <- "#CB2026"
H3K4me1_colors <-  "#642165"
H3K4me3_colors <- "#B89BC9"
H3K56ac_colors <- "#99CA3C"
ATAC_colors <- "#FFCD05"
H2A_colors <-  "#66c2a5"
H3K27me3_colors <- "#234bb8"

colors_to_use <- list(H3K27me3_colors,ATAC_colors,H3K56ac_colors,H3K36me3_colors,H3K4me1_colors,H3K4me3_colors,H3K56ac_colors)




# Pablo Theme -------------------------------------------------------------
pablo_lncRNA_theme <- function(base_size = 25,
                               base_family = "Arial",
                               base_line_size = base_size / 170,
                               base_rect_size = base_size / 170){
  theme_bw(base_size = base_size, 
           base_family = base_family,
           base_line_size = base_line_size) %+replace%
    theme(
      plot.title = element_text(
        color = "black", 
        face = "bold",
        hjust = .5, 
        size = rel(1.3)),
      axis.title = element_text(
        color = "black",
        size = rel(1)),
      axis.text = element_text(
        color = "black",
        size = rel(1)),
      legend.text = element_text(
        color = "black",
        size = rel(1)),
      #Facets
      strip.text = element_text(
        color = "black",
        size=rel(1)),
      
      
      panel.grid.major = element_line(colour = "grey", linetype = "solid"),   
      panel.grid.minor = element_line(colour = "grey",
                                      linetype = "solid", 
                                      size = rel(1)),  
      text = element_text(family = "Arial",, face = "plain",
                          color = "black", size = base_size,
                          hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                          margin = margin(), debug = FALSE),
      
      complete = TRUE
    )
}






# Read in Files -----------------------------------------------------------




read_in_file <- function(iso_pass, iso_fail, WGS_pass, WGS_fail, type_name){
  
 `%notin%` <- Negate(`%in%`)
  
  iso_fail_read <- read_delim(iso_fail, '\t', col_names = c("gene_name"))# %>% mutate(class_name = "iso_fail") %>% mutate(overall_class = "fail", type_seq = "iso_seq") 
  iso_pass_read <- read_delim(iso_pass, '\t', col_names = c("gene_name"))# %>% mutate(class_name = "iso_pass") %>% mutate(overall_class = "pass", type_seq = "iso_seq") 
  WGS_pass_read <- read_delim(WGS_pass, '\t', col_names = c("gene_name"))# %>% mutate(class_name = "WGS_pass") %>% mutate(overall_class = "pass", type_seq = "WGS") 
  WGS_fail_read <- read_delim(WGS_fail, '\t', col_names = c("gene_name"))# %>% mutate(class_name = "WGS_fail") %>% mutate(overall_class = "fail", type_seq = "WGS") 
  
  
  iso_WGS_backed <- sum(iso_pass_read$gene_name %in% WGS_pass_read$gene_name)
  ISO_only <- sum(!iso_pass_read$gene_name %in% WGS_pass_read$gene_name)
  WGS_only <- sum(WGS_pass_read$gene_name %notin% iso_pass_read$gene_name)
  unsupported_annotaiton <- sum(iso_fail_read$gene_name %in% WGS_fail_read$gene_name)
  
  
  #all_combined <- bind_rows(iso_fail_read,
  #iso_pass_read,
  #WGS_pass_read,
  #WGS_fail_read)
  
  counts <- c(iso_WGS_backed,
  ISO_only,
  WGS_only,
  unsupported_annotaiton)
  
  
  type <- c(rep(type_name, 4))
  condition <- c("Iso-seq & RNA-seq", "Iso-seq", 'RNA-seq', "Unsupported")
  
  data <- data.frame(type, condition, counts)
  
  return(data)
}




hyper_large_genes <- read_in_file( "06.iso_seq_passing_assemblies/hyper_large_genes_updated_passing.txt",
                          "06.iso_seq_passing_assemblies/hyper_large_genes_updated_failing.txt",
                          "07.WGS_passing_assemblies/hyper_large_genes_updated_passing.txt",
                          "07.WGS_passing_assemblies/hyper_large_genes_updated_failing.txt",
                          "Hyper Large")

major_extension_genes <- read_in_file("06.iso_seq_passing_assemblies/major_extension_genes_updated_passing.txt",
                                      "06.iso_seq_passing_assemblies/major_extension_genes_updated_failing.txt",
                                      "07.WGS_passing_assemblies/major_extension_genes_updated_passing.txt",
                                      "07.WGS_passing_assemblies/major_extension_genes_updated_failing.txt",
                                      "Major Extension")

minor_extension_genes <- read_in_file("06.iso_seq_passing_assemblies/minor_extension_genes_updated_passing.txt",
                                      "06.iso_seq_passing_assemblies/minor_extension_genes_updated_failing.txt",
                                      "07.WGS_passing_assemblies/minor_extension_genes_updated_passing.txt",
                                      "07.WGS_passing_assemblies/minor_extension_genes_updated_failing.txt",
                                      "Minor Extension")

merged_genes <- read_in_file("06.iso_seq_passing_assemblies/merged_genes_updated_passing.txt",
                             "06.iso_seq_passing_assemblies/merged_genes_updated_failing.txt",
                             "07.WGS_passing_assemblies/merged_genes_updated_passing.txt",
                             "07.WGS_passing_assemblies/merged_genes_updated_failing.txt",
                             "Merged")

all_groups <- bind_rows(hyper_large_genes,
major_extension_genes,
minor_extension_genes,
merged_genes)

all_groups

final_plotting_group <- all_groups %>% 
  dplyr::group_by(type) %>% 
  mutate(total_counts = sum(counts)) %>% 
  mutate(plot_proportion = counts / total_counts) %>% 
  mutate(write_total_n = str_c("N=", as.character(total_counts)))

all_groups %>% 
  dplyr::group_by(type) %>% 
  mutate(total_counts = sum(counts)) %>% 
  mutate(Prop_class = counts/total_counts * 100)
  

final_plotting_group$condition

annotaion_support <- ggplot(final_plotting_group, aes(fill=factor(condition, levels=c("Iso-seq & RNA-seq", "Iso-seq","RNA-seq","Unsupported" )) , y=counts, x=type, label=counts)) + 
  geom_bar(position="fill", stat="identity") + pablo_lncRNA_theme() +
  geom_text(aes(x = type, y = plot_proportion, label = counts), position=position_fill(vjust=0.5)) + 
  geom_text(aes(x = type, y = 1.02 , label = write_total_n)) + theme_cowplot() + scale_fill_manual(values = alpha(colors_to_use, .8)) +
  xlab("Annotation Type") + ylab("Proportion") + labs(fill = "Annotation Support") + ggtitle("Proportion of Annotations Supported") + 
  theme(legend.position="top")

annotaion_support

ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure3_5/imgs/annotation_support/", filename="UpdateAnnotation_support.svg", plot=annotaion_support, width = 10, height = 8, units = "in")


View(final_plotting_group)

final_plotting_group

test <- final_plotting_group %>% 
  mutate(passing_failing_group = case_when(condition == "Iso-seq & RNA-seq" ~ 1,
            condition =="Iso-seq" ~ 1, 
            condition =='RNA-seq' ~ 1, 
            condition == "Unsupported" ~ 0)) %>% 
  group_by(passing_failing_group) %>% 
  summarise(total_number = sum(counts))

test
  
6385 + 2828

6385/9213
