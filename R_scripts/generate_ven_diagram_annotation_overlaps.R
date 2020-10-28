library(VennDiagram)
library("extrafont")
library(patchwork)
library(dplyr)
library(cowplot)
library(eulerr)
library(ggplot2)
library(devtools)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

setwd("/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure3_5")

read_in_annot_counts <- function(file_name, annotation_type){
  
  
  anno_file <- read_delim(file_name,'\t', col_names = TRUE) %>% 
    mutate(anno_type = annotation_type) 
  
  return(anno_file)
  
}
ven_diagram_drawer <- function(read_in_counts){
  
  
  p_venn <- draw.triple.venn(area1 = read_in_counts$leaf, 
                             area2 = read_in_counts$ear, 
                             area3 = read_in_counts$root, 
                             n12 = read_in_counts$ear_leaf, 
                             n23 = read_in_counts$ear_root, 
                             n13 = read_in_counts$leaf_root, 
                             n123 = read_in_counts$all,
                             fontfamily = "sans",
                             category = c("Leaf", "Ear", "Root"),
                             #Unusre how this works below, but allows me to operate on Gtree object 
                             #From https://stackoverflow.com/questions/21202876/adding-title-and-sub-title-to-venn-diagram 
                             ind = FALSE,
                             cex = 2,
                             cat.cex = 2)
  
  return(p_venn)
  
}



hyper_large_counts <- "01.generate_annotation_ven_diagrams/hyper_large_counts.txt"
major_extension_counts <- "01.generate_annotation_ven_diagrams/major_extension_counts.txt"
merged_genes_counts <- "01.generate_annotation_ven_diagrams/merged_genes_counts.txt"
minor_extension_counts <- "01.generate_annotation_ven_diagrams/minor_extension_counts.txt"
novel_counts <- "01.generate_annotation_ven_diagrams/novel_counts.txt"


hyper_large <- read_in_annot_counts(hyper_large_counts, "Hyper Large Gene Extensions")
major_extension <- read_in_annot_counts(major_extension_counts, "Major Gene Extensions")
merged_genes <- read_in_annot_counts(merged_genes_counts, "Merged Genes")
minor_extension <- read_in_annot_counts(minor_extension_counts, "Minor Gene Extension")
novel <- read_in_annot_counts(novel_counts, "Novel Genes")

minor_extension



hyper_large_ven <- ven_diagram_drawer(hyper_large)
minor_extension_ven <- ven_diagram_drawer(minor_extension)
major_extension_ven <- ven_diagram_drawer(major_extension)
merged_genes_ven <- ven_diagram_drawer(merged_genes)
novel_ven <- ven_diagram_drawer(novel)



hyper_large_ven_title <- grid.arrange(gTree(children=hyper_large_ven),top=textGrob(hyper_large$anno_type, gp=gpar(fontsize=25,font=1)))
minor_extension_ven_title <- grid.arrange(gTree(children=minor_extension_ven),top=textGrob(minor_extension$anno_type, gp=gpar(fontsize=25,font=1)))
major_extension_ven_title <- grid.arrange(gTree(children=major_extension_ven),top=textGrob(major_extension$anno_type, gp=gpar(fontsize=25,font=1)))
merged_genes_ven_title <- grid.arrange(gTree(children=merged_genes_ven),top=textGrob(merged_genes$anno_type, gp=gpar(fontsize=25,font=1)))
novel_ven_title <- grid.arrange(gTree(children=novel_ven),top=textGrob(novel$anno_type, gp=gpar(fontsize=25,font=1)))
minor_extension_ven_plot

hyper_large_ven_plot <- ggdraw(hyper_large_ven_title) +theme(plot.background = element_rect(fill = NA),plot.margin = margin(12, 12, 12, 12)) 
minor_extension_ven_plot <- ggdraw(minor_extension_ven_title) +theme(plot.background = element_rect(fill = NA),plot.margin = margin(12, 12, 12, 12))
major_extension_ven_plot <- ggdraw(major_extension_ven_title) +theme(plot.background = element_rect(fill = NA),plot.margin = margin(12, 12, 12, 12))
merged_genes_ven_plot <- ggdraw(merged_genes_ven_title) +theme(plot.background = element_rect(fill = NA),plot.margin = margin(12, 12, 12, 12))
novel_ven_plot <- ggdraw(novel_ven_title) +theme(plot.background = element_rect(fill = NA),plot.margin = margin(12, 12, 12, 12))

ggsave(plot= hyper_large_ven_plot, filename = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure3_5/imgs/hyper_large_ven.svg", width = 7, height =6 )
ggsave(plot= minor_extension_ven_plot, filename = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure3_5/imgs/minor_extension_ven.svg", width = 7, height =6 )
ggsave(plot= major_extension_ven_plot, filename = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure3_5/imgs/major_extension_ven.svg", width = 7, height =6 )
ggsave(plot= merged_genes_ven_plot, filename = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure3_5/imgs/merged_genes_ven.svg", width = 7, height =6 )
ggsave(plot= novel_ven_plot, filename = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure3_5/imgs/novel_genes_ven.svg", width = 7, height =6 )
