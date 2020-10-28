# 5/20/2020 
#This script generates upset plots which look at the expressed genes, and how they overlap
#various histone modoifications. The point these upset plots are trying to make is that 
#genes which overlap both of these modification types are more expressed - AND most genes 
#in the genome when controlling for mappability have both of these modifications.



library(tidylog)
library(tidyr)
library(dplyr)
library(tidyverse)
library(UpSetR)
library(ggplot2)
library(grid)
library(plyr)
library(ComplexHeatmap)


assign_overlap_yes_no <- function(histone_mod_tribble) {
  
  histone_mod_tribble_passing <- histone_mod_tribble %>% 
    mutate(modification_pass_fail = case_when(score == '.' ~ 0,
                                              score != '.' ~ 1))
  
  return(histone_mod_tribble_passing)
}
alternative_TPM_file <- function(TPM_file) {
  col_names_TPM_file <- c("gene_ID", "Chr", "Start", "End", "Length", "Reads", "TPM")
  read_TPM_gene <- read_delim(TPM_file, '\t', col_names = col_names_TPM_file) %>% 
    mutate(updated_name = str_c("gene", gene_ID , sep =':')) %>% 
    select(-gene_ID) %>% 
    dplyr::rename(gene_ID = updated_name) %>% 
    select(gene_ID, TPM)
  
  
  read_TPM_gene = read_TPM_gene[-1,]
  
  pound_sign_carriers <- read_TPM_gene %>% 
    dplyr::filter(stringr::str_detect(gene_ID, '#')) %>% 
    dplyr::mutate(sample = gsub("#.*", "", gene_ID)) %>% 
    dplyr::group_by(sample) %>%
    dplyr::summarize(new_TPM = mean(TPM)) %>% 
    dplyr::rename(gene_ID = sample) %>% 
    dplyr::rename(TPM = new_TPM)
  
  
  read_TPM_gene_2 <- bind_rows(read_TPM_gene,pound_sign_carriers)
  
  read_TPM_gene_2$TPM <- as.numeric(read_TPM_gene_2$TPM)
  
  
  return(read_TPM_gene_2)
}








combine_overlap_TPM <- function(TPM_file, K36_file, K4me1_file, K4me3_file, K56ac_file, mappability_file){
  

  col_names_TPM_file <- c("gene_ID", "TPM")
  read_TPM_gene <- alternative_TPM_file(TPM_file)
  
  mapp_hit_cols <- c("chrom", "start", "stop", "gene_ID", "NA", "strand", "origin", "bio_ID", "other","other_2", "mapp")
  mappability <- read_delim(mappability_file, '\t', col_names = mapp_hit_cols, col_types = "ccccccccccd" ) %>% select(gene_ID, mapp)
  
  bed_hit_col_names <- c("chrom", "start", "stop", "gene_ID", "NA", "strand", "origin", "bio_ID","none", "ID_string", "score", "hit_chrom", "hit_start", "hit_stop")
  gene_H3K36me3 <- read_delim(K36_file, '\t', col_names = bed_hit_col_names) 
  gene_H3K4me1 <- read_delim(K4me1_file, '\t', col_names = bed_hit_col_names) 
  gene_H3K4me3 <- read_delim(K4me3_file, '\t', col_names = bed_hit_col_names) 
  gene_H3K56ac <- read_delim(K56ac_file, '\t', col_names = bed_hit_col_names) 
  
  
  
  gene_H3K36me3_pass_fail <- assign_overlap_yes_no(gene_H3K36me3) %>% select(gene_ID, modification_pass_fail) %>%  dplyr::rename("H3K36me3" = modification_pass_fail) 
  gene_H3K4me1_pass_fail <- assign_overlap_yes_no(gene_H3K4me1) %>% select(gene_ID, modification_pass_fail) %>%  dplyr::rename("H3K4me1" = modification_pass_fail) 
  gene_H3K4me3_pass_fail <- assign_overlap_yes_no(gene_H3K4me3) %>% select(gene_ID, modification_pass_fail) %>%  dplyr::rename("H3K4me3" = modification_pass_fail) 
  gene_H3K56ac_pass_fail <- assign_overlap_yes_no(gene_H3K56ac) %>% select(gene_ID, modification_pass_fail) %>%  dplyr::rename("H3K56ac" = modification_pass_fail)
  
  
  gene_all_combined <- left_join(read_TPM_gene, gene_H3K36me3_pass_fail,  by = c("gene_ID")) %>% 
    left_join(., gene_H3K4me1_pass_fail,  by = c("gene_ID")) %>% 
    left_join(., gene_H3K4me3_pass_fail,  by = c("gene_ID")) %>% 
    left_join(., gene_H3K56ac_pass_fail,  by = c("gene_ID")) %>% 
    filter(is.na(H3K36me3) == FALSE & is.na(H3K4me1) == FALSE & is.na(H3K4me3) == FALSE &  is.na(H3K56ac) == FALSE) %>% 
    left_join(., mappability, by = c("gene_ID")) %>% 
    filter(mapp > .65)
  
  return(gene_all_combined)
  
}

setwd("/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2")





leaf_overlap_expression <- combine_overlap_TPM(
"00.data/TPM_vals/RNA_B73_leaf_TPM.surviving.txt",
"00.data/histone_mods_intersecting_genes/all_genes_leaf_overlapping_H3K36me3_broad.bed",
"00.data/histone_mods_intersecting_genes/all_genes_leaf_overlapping_H3K4me1_broad.bed",
"00.data/histone_mods_intersecting_genes/all_genes_leaf_overlapping_H3K4me3_narrow.bed",
"00.data/histone_mods_intersecting_genes/all_genes_leaf_overlapping_H3K56ac_narrow.bed",
"00.data/genome_annotation_mappability.values.bed")





# Testing section -- Looking at expressed genes with NO marks -------------
#bed_hit_col_names <- c("chrom", "start", "stop", "gene_ID", "NA", "strand", "origin", "bio_ID","none", "ID_string", "score", "hit_chrom", "hit_start", "hit_stop")
#test_this <- read_delim("00.data/histone_mods_intersecting_genes/genes_leaf_overlapping_H3K36me3_broad.bed", '\t', col_names = bed_hit_col_names)
#test_this2 <- assign_overlap_yes_no(test_this)#

#test_this2 %>% 
#  dplyr::group_by(modification_pass_fail) %>% 
#  dplyr::summarise(counts = n())#

#leaf_overlap_expression %>% 
#  select(H3K36me3) %>% 
#  dplyr::group_by(H3K36me3) %>% 
#  dplyr::summarise(counts = n())
#  
#leaf_overlap_expression %>% 
#  filter(is.na(H3K36me3) == TRUE) %>% 
#  View()#
#

#col_names_TPM_file <- c("gene_ID", "TPM")
#read_TPM_gene <- read_delim("00.data/TPM_vals/RNA_B73_leaf_TPM.surviving.txt", '\t', col_names = col_names_TPM_file) %>% 
#  mutate(updated_name = str_c("gene", gene_ID , sep =':')) %>% 
#  select(-gene_ID) %>% 
#  dplyr::rename(gene_ID = updated_name)#

#final$modification_pass_fail
#final <- full_join(read_TPM_gene, test_this2, by = "gene_ID") 
#idk <- final %>% 
#  filter(is.na(modification_pass_fail) == TRUE)


# Load Other Datasets ---------------------------------------------------------------


root_overlap_expression <- combine_overlap_TPM(
  "00.data/TPM_vals/RNA_B73_root_TPM.surviving.txt",
  "00.data/histone_mods_intersecting_genes/all_genes_root_overlapping_H3K36me3_broad.bed",
  "00.data/histone_mods_intersecting_genes/all_genes_root_overlapping_H3K4me1_broad.bed",
  "00.data/histone_mods_intersecting_genes/all_genes_root_overlapping_H3K4me3_narrow.bed",
  "00.data/histone_mods_intersecting_genes/all_genes_root_overlapping_H3K56ac_narrow.bed",
  "00.data/genome_annotation_mappability.values.bed")


ear_overlap_expression <- combine_overlap_TPM(
  "00.data/TPM_vals/RNA_B73_ear_TPM.surviving.txt",
  "00.data/histone_mods_intersecting_genes/all_genes_ear_overlapping_H3K36me3_broad.bed",
  "00.data/histone_mods_intersecting_genes/all_genes_ear_overlapping_H3K4me1_broad.bed",
  "00.data/histone_mods_intersecting_genes/all_genes_ear_overlapping_H3K4me3_narrow.bed",
  "00.data/histone_mods_intersecting_genes/all_genes_ear_overlapping_H3K56ac_narrow.bed",
  "00.data/genome_annotation_mappability.values.bed")



generate_K36me3_K4me1ation_K56ac_K4me3iaion_classes <- function(TPM_histone_mod_intersect){
  
  final <- TPM_histone_mod_intersect %>% 
  mutate(Grouping = case_when( TPM >= 1 ~ "Expressed",
                               TPM < 1 ~ "Not Expressed")) %>% 
    mutate(K36me3_K4me1ation_mods = H3K36me3 + H3K4me1) %>%
    mutate(iniation_mods =  H3K56ac + H3K4me3) %>% 
    mutate(K56ac_K4me3 = case_when( iniation_mods == 2 ~ 1,
                             iniation_mods < 2 ~ 0)) %>% 
    mutate(K36me3_K4me1 = case_when(K36me3_K4me1ation_mods == 2 ~ 1,
                             K36me3_K4me1ation_mods < 2 ~ 0)) %>% 
    filter(Grouping == "Expressed") %>% 
    #Count either inclusion of 1 mark as "K56ac_K4me3iation" or "K36me3_K4me1ation" 
    mutate(K56ac_K4me3iation = case_when( iniation_mods >= 1 ~ 1,
                                   iniation_mods < 1 ~ 0)) %>% 
    mutate(K36me3_K4me1ation = case_when(K36me3_K4me1ation_mods >= 1 ~ 1,
                                  K36me3_K4me1ation_mods < 1 ~ 0)) %>% 
    unique()
  
  return(final)
  
}


leaf_gene_all_combined_gr_1 <- generate_K36me3_K4me1ation_K56ac_K4me3iaion_classes(leaf_overlap_expression) 
root_gene_all_combined_gr_1 <- generate_K36me3_K4me1ation_K56ac_K4me3iaion_classes(root_overlap_expression)
ear_gene_all_combined_gr_1 <- generate_K36me3_K4me1ation_K56ac_K4me3iaion_classes(ear_overlap_expression)



# Old Upset plot generation using default Upset plot  --------
#Old Upset plot branch is no longer under active developemnt, and the author is
#not activly developing the package
#things <- c("K56ac_K4me3iation","K36me3_K4me1ation")

#generate_upset_R2 <- function(tissue_tribble) {
#  
#  transformed_df <- tissue_tribble %>% 
#    mutate(Log2_TPM = log2(TPM))
#  
#  data_frame_for_plotting <- as.data.frame(transformed_df)
#  generated_graph <- upset(data_frame_for_plotting, boxplot.summary = c("Log2_TPM"), sets = c("K36me3_K4me1ation", "K56ac_K4me3iation"), mb.ratio = c(0.55, 0.45), order.by = "freq",
#                           text.scale = c(2.5,2.5,2,1.5,2.5,2.5),  point.size = 4, line.size = 2, empty.intersections = "on")
#  
#  return(generated_graph)
#}
#
#pdf(file="imgs/upset_plots/leaf_upset.pdf", onefile=FALSE, width = 12, height = 9) # or other device
#leaf_upset_final <- generate_upset_R2(leaf_gene_all_combined_gr_1)
#leaf_upset_final
#dev.off()
#
#
#pdf(file="imgs/upset_plots/root_upset.pdf", onefile=FALSE, width = 12, height = 9) # or other device
#root_upset_final <- generate_upset_R2(root_gene_all_combined_gr_1)
#root_upset_final
#dev.off()
#
#pdf(file="imgs/upset_plots/ear_upset.pdf", onefile=FALSE, width = 12, height = 9) # or other device
#ear_upset_final <- generate_upset_R2(ear_gene_all_combined_gr_1)
#ear_upset_final
#dev.off()


# Try ComplexHeatMap Upset Plot -------------------------------------------

#Generate the df to be used
generate_complex_heatmap_matrix <- function(tissue_name_tribble){
  
  #Take the log2 of the tpm calye
  transformed_df <- tissue_name_tribble %>% 
    mutate(Log2_TPM = log2(TPM + 1))
  
  upset_R_df <- transformed_df %>% 
    select(gene_ID, Log2_TPM, K36me3_K4me1ation, K56ac_K4me3iation) %>% 
    unique()
  
  samp2 <- as.data.frame(upset_R_df)[,-1]
  
  #On-the-fly indicator function for use in formulae
  #leaf_comb = make_comb_mat(samp2, remove_complement_set =TRUE)
  return(samp2)
  
}
#Run the complex heatmap matrix calculation, as well as generate upset plot - return
generate_complex_heatmap_UpSetR_plot <- function(data_frame_filtered) { 
 
  comb_value <- make_comb_mat(data_frame_filtered, remove_complement_set = FALSE)
  comb_elements_leaf <- lapply(comb_name(comb_value), function(nm) extract_comb(comb_value, nm))

  
  #Generate  the list to used in the Annotation Density plots on either the top or bottom of the values
  TPM_val = lapply(comb_elements_leaf, function(ind) data_frame_filtered$Log2_TPM[ind])
  print(TPM_val)
  
  #Putting the Expression information on the BOTTOM of the Plot
  final <- UpSet(comb_value, #gp = gpar(fontfamily = "sans"),
                 top_annotation = upset_top_annotation(comb_value, 
                                                       annotation_name_rot = 90,
                                                       annotation_name_side = "left",
                                                       axis_param = list(side = "left")),
                 right_annotation = NULL)  %v% 
  HeatmapAnnotation("Log2(TPM +1)" = anno_density(TPM_val, type = "violin"), annotation_name_rot = 90, annotation_name_side = "left")

  return(final)
}
#Add the annotations in forms of counts
append_annotations <- function(heatmap, data_frame_filtered){
  
  comb_val <- make_comb_mat(data_frame_filtered, remove_complement_set =FALSE)
  od = column_order(heatmap)
  cs = comb_size(comb_val)
  decorate_annotation("Intersection\nsize", {
    grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
              default.units = "native", just = "bottom", gp = gpar(fontsize = 8))
  })
  
}

library(cowplot)


#Generate Leaf Upset plot
leaf_com_matrix <- generate_complex_heatmap_matrix(leaf_gene_all_combined_gr_1)
leaf_upset_plot <- generate_complex_heatmap_UpSetR_plot(leaf_com_matrix)
pdf(file="imgs/upset_plots/leaf_upset.pdf", onefile=FALSE, width = 4.5, height = 4) # or other device
leaf_upset_plot
append_annotations(leaf_upset_plot, leaf_com_matrix)
dev.off()

root_com_matrix <- generate_complex_heatmap_matrix(root_gene_all_combined_gr_1)
root_upset_plot <- generate_complex_heatmap_UpSetR_plot(root_com_matrix)
pdf(file="imgs/upset_plots/root_upset.pdf", onefile=FALSE, width = 4.5, height = 4) # or other device
root_upset_plot
append_annotations(root_upset_plot, root_com_matrix)
dev.off()


ear_com_matrix <- generate_complex_heatmap_matrix(ear_gene_all_combined_gr_1)
ear_upset_plot <- generate_complex_heatmap_UpSetR_plot(ear_com_matrix)
pdf(file="imgs/upset_plots/ear_upset.pdf", onefile=FALSE, width = 4.5, height = 4) # or other device
ear_upset_plot
append_annotations(ear_upset_plot, ear_com_matrix)
dev.off()


#Based off these upset plots, I need to look at the classes which have only K56ac_K4me3iation, or K36me3_K4me1ation marks
#Need to see if this is biological, or just an oddity that I'm picking up.
#Below is a function to output the bed files of those regions which only have K56ac_K4me3iation/K36me3_K4me1ation.




###########################################s###############
###########################################s###############
###########################################s###############
###########################################s###############
# Same Analysis, but split by all marks -----------------------------------

#Generate the df to be used
generate_complex_heatmap_matrix_split_mods<- function(tissue_name_tribble){
  
  #Take the log2 of the tpm calye
  transformed_df <- tissue_name_tribble %>% 
    mutate(Log2_TPM = log2(TPM + 1))
  
  upset_R_df <- transformed_df %>% 
    select(gene_ID, Log2_TPM, H3K36me3, H3K4me1, H3K56ac, H3K4me3) %>% 
    unique()
  
  samp2 <- as.data.frame(upset_R_df)[,-1]
  
  #On-the-fly indicator function for use in formulae
  #leaf_comb = make_comb_mat(samp2, remove_complement_set =TRUE)
  return(samp2)
  
}
generate_complex_heatmap_UpSetR_plot_split_mods <- function(data_frame_filtered) { 
  
  comb_value <- make_comb_mat(data_frame_filtered, remove_complement_set = FALSE)
  comb_elements_leaf <- lapply(comb_name(comb_value), function(nm) extract_comb(comb_value, nm))
  
  
  #Generate  the list to used in the Annotation Density plots on either the top or bottom of the values
  TPM_val = lapply(comb_elements_leaf, function(ind) data_frame_filtered$Log2_TPM[ind])
  print(TPM_val)
  
  #Putting the Expression information on the BOTTOM of the Plot
  final <- UpSet(comb_value, #gp = gpar(fontfamily = "sans"),
                 top_annotation = upset_top_annotation(comb_value, 
                                                       annotation_name_rot = 90,
                                                       annotation_name_side = "left",
                                                       axis_param = list(side = "left")),
                 right_annotation = NULL)  %v% 
    HeatmapAnnotation("Log2(TPM +1)" = anno_density(TPM_val, type = "violin"), annotation_name_rot = 90, annotation_name_side = "left")
  
  return(final)
}
#Run the complex heatmap matrix calculation, as well as generate upset plot - return

#Generate Leaf Upset plot
leaf_com_matrix <- generate_complex_heatmap_matrix_split_mods(leaf_gene_all_combined_gr_1)
leaf_upset_plot <- generate_complex_heatmap_UpSetR_plot(leaf_com_matrix)


pdf(file="imgs/upset_plots/leaf_upset.split_mods.pdf", onefile=FALSE, width = 4.5, height = 4) # or other device
leaf_upset_plot
append_annotations(leaf_upset_plot, leaf_com_matrix)
dev.off()



root_com_matrix <- generate_complex_heatmap_matrix_split_mods(root_gene_all_combined_gr_1)
root_upset_plot <- generate_complex_heatmap_UpSetR_plot(root_com_matrix)
pdf(file="imgs/upset_plots/root_upset.split_mods.pdf", onefile=FALSE, width = 4.5, height = 4) # or other device
root_upset_plot
append_annotations(root_upset_plot, root_com_matrix)
dev.off()


ear_com_matrix <- generate_complex_heatmap_matrix_split_mods(ear_gene_all_combined_gr_1)
ear_upset_plot <- generate_complex_heatmap_UpSetR_plot(ear_com_matrix)
pdf(file="imgs/upset_plots/ear_upset.split_mods.pdf", onefile=FALSE, width = 4.5, height = 4) # or other device
ear_upset_plot
append_annotations(ear_upset_plot, ear_com_matrix)
dev.off()



#################################/########################3/#?3################/########
#################################/########################3/#?3################/########
#################################/########################3/#?3################/########
#################################/########################3/#?3################/########

generate_complex_heatmap_matrix_independent_mod <- function(tissue_name_tribble){
  
  #Take the log2 of the tpm calye
  transformed_df <- tissue_name_tribble %>% 
    mutate(Log2_TPM = log2(TPM + 1))
  
  upset_R_df <- transformed_df %>% 
    select(gene_ID, Log2_TPM, K56ac_K4me3, K36me3_K4me1) %>% 
    unique()
  
  samp2 <- as.data.frame(upset_R_df)[,-1]
  
  #On-the-fly indicator function for use in formulae
  #leaf_comb = make_comb_mat(samp2, remove_complement_set =TRUE)
  return(samp2)
}


#Run the complex heatmap matrix calculation, as well as generate upset plot - return
#Generate Leaf Upset plot
leaf_com_matrix <- generate_complex_heatmap_matrix_independent_mod(leaf_gene_all_combined_gr_1)
leaf_upset_plot <- generate_complex_heatmap_UpSetR_plot_split_mods(leaf_com_matrix)
pdf(file="imgs/upset_plots/leaf_upset.both_mods_req.pdf", onefile=FALSE, width = 4.5, height = 4) # or other device
leaf_upset_plot
append_annotations(leaf_upset_plot, leaf_com_matrix)
dev.off()

root_com_matrix <- generate_complex_heatmap_matrix_independent_mod(root_gene_all_combined_gr_1)
root_upset_plot <- generate_complex_heatmap_UpSetR_plot_split_mods(root_com_matrix)
pdf(file="imgs/upset_plots/root_upset.both_mods_req.pdf", onefile=FALSE, width = 4.5, height = 4) # or other device
root_upset_plot
append_annotations(root_upset_plot, root_com_matrix)
dev.off()

ear_com_matrix <- generate_complex_heatmap_matrix_independent_mod(ear_gene_all_combined_gr_1)
ear_upset_plot <- generate_complex_heatmap_UpSetR_plot_split_mods(ear_com_matrix)
pdf(file="imgs/upset_plots/ear_upset.both_mods_req.pdf", onefile=FALSE, width = 4.5, height = 4) # or other device
ear_upset_plot
append_annotations(ear_upset_plot, ear_com_matrix)
dev.off()





#Write bed files of those genes which are INTIATION ONLY or K36me3_K4me1ATION ONLY based off the 
#the inclusiong of BOTH histone mods of each type, and NO other overlap with otehr Histone mods

#Silly, but based off the way I did the intial read in, have to re-generate a list of all genes in each tissue
gather_gene_list <- function(K36_file, K4me1_file, K4me3_file, K56ac_file){
  
  
  
  bed_hit_col_names <- c("chrom", "start", "stop", "gene_ID", "score", "strand", "origin", "bio_ID","none", "ID_string", "score", "hit_chrom", "hit_start", "hit_stop")
  gene_H3K36me3 <- read_delim(K36_file, '\t', col_names = bed_hit_col_names) 
  gene_H3K4me1 <- read_delim(K4me1_file, '\t', col_names = bed_hit_col_names) 
  gene_H3K4me3 <- read_delim(K4me3_file, '\t', col_names = bed_hit_col_names) 
  gene_H3K56ac <- read_delim(K56ac_file, '\t', col_names = bed_hit_col_names) 
  
  
  
  gene_H3K36me3_pass_fail <- assign_overlap_yes_no(gene_H3K36me3) %>% select(chrom:strand, modification_pass_fail) %>%  dplyr::rename("H3K36me3" = modification_pass_fail) 
  gene_H3K4me1_pass_fail <- assign_overlap_yes_no(gene_H3K4me1) %>% select(chrom:strand, modification_pass_fail) %>%  dplyr::rename("H3K4me1" = modification_pass_fail) 
  gene_H3K4me3_pass_fail <- assign_overlap_yes_no(gene_H3K4me3) %>% select(chrom:strand, modification_pass_fail) %>%  dplyr::rename("H3K4me3" = modification_pass_fail) 
  gene_H3K56ac_pass_fail <- assign_overlap_yes_no(gene_H3K56ac) %>% select(chrom:strand, modification_pass_fail) %>%  dplyr::rename("H3K56ac" = modification_pass_fail)
  
  
  gene_all_combined <- full_join(gene_H3K36me3_pass_fail, gene_H3K4me1_pass_fail,  by = c("chrom", "start", "stop", "gene_ID", "score", "strand")) %>% 
    full_join(., gene_H3K4me3_pass_fail,  by = c("chrom", "start", "stop", "gene_ID", "score", "strand")) %>% 
    full_join(., gene_H3K56ac_pass_fail,  by = c("chrom", "start", "stop", "gene_ID", "score", "strand")) %>% 
    distinct()
  
  return(gene_all_combined)
}

#pull the correct gene names from the finalized gene list
write_bed_table_with_names <- function(origina_file_name, combined_genes, out_file_name_base){
  
  
  K36me3_K4me1ation_bed_name = str_c(out_file_name_base, "_elongation_only.bed")
  K56ac_K4me3iation_bed_name = str_c(out_file_name_base, "_initiation_only.bed")
  

  
  K36me3_K4me1ation_only <- combined_genes %>% 
    dplyr::filter((H3K36me3 + H3K4me1) == 2)  %>% 
    dplyr::filter((H3K4me3 + H3K56ac) == 0 ) %>% 
    distinct()
  
  K56ac_K4me3iation_only <- combined_genes %>% 
    dplyr::filter((H3K36me3 + H3K4me1) == 0)  %>% 
    dplyr::filter((H3K4me3 + H3K56ac) == 2 ) %>% 
    distinct()
  
  print(dim(K36me3_K4me1ation_only))
  print(dim(K56ac_K4me3iation_only))
  
  K36me3_K4me1ation_only_bed <- origina_file_name %>% 
    filter(gene_ID %in% K36me3_K4me1ation_only$gene_ID)
  
  K56ac_K4me3iation_only_bed <- origina_file_name %>% 
    filter(gene_ID %in% K56ac_K4me3iation_only$gene_ID)


  write.table(K36me3_K4me1ation_only_bed, file=K36me3_K4me1ation_bed_name, sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)  
  write.table(K56ac_K4me3iation_only_bed, file=K56ac_K4me3iation_bed_name, sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)  
  
}





#Load Tissue Specific bed Files
leaf_total_gene_list <- gather_gene_list("00.data/histone_mods_intersecting_genes/genes_leaf_overlapping_H3K36me3_broad.bed",
  "00.data/histone_mods_intersecting_genes/genes_leaf_overlapping_H3K4me1_broad.bed",
  "00.data/histone_mods_intersecting_genes/genes_leaf_overlapping_H3K4me3_narrow.bed",
  "00.data/histone_mods_intersecting_genes/genes_leaf_overlapping_H3K56ac_narrow.bed")

root_total_gene_list <- gather_gene_list("00.data/histone_mods_intersecting_genes/genes_root_overlapping_H3K36me3_broad.bed",
  "00.data/histone_mods_intersecting_genes/genes_root_overlapping_H3K4me1_broad.bed",
  "00.data/histone_mods_intersecting_genes/genes_root_overlapping_H3K4me3_narrow.bed",
  "00.data/histone_mods_intersecting_genes/genes_root_overlapping_H3K56ac_narrow.bed")


ear_total_gene_list <- gather_gene_list("00.data/histone_mods_intersecting_genes/genes_ear_overlapping_H3K36me3_broad.bed",
  "00.data/histone_mods_intersecting_genes/genes_ear_overlapping_H3K4me1_broad.bed",
  "00.data/histone_mods_intersecting_genes/genes_ear_overlapping_H3K4me3_narrow.bed",
  "00.data/histone_mods_intersecting_genes/genes_ear_overlapping_H3K56ac_narrow.bed")


#I generally move these files over to '02.isolated_markes', a directory within ~/Projects/03.ncRNA_project/03.Figures/Figure2
write_bed_table_with_names(leaf_total_gene_list, leaf_gene_all_combined_gr_1, "leaf")
write_bed_table_with_names(root_total_gene_list, root_gene_all_combined_gr_1, "root")
write_bed_table_with_names(ear_total_gene_list, ear_gene_all_combined_gr_1, "ear")


# Notes -------------------------------------------------------------------

#NOTES ON GENERATING AN UPSET PLOT IN COMPLEXHEATMAP
#

#Take the log2 of the tpm calye
#transformed_df <- leaf_gene_all_combined_gr_1 %>% 
#  mutate(Log2_TPM = log2(TPM + 1))
#
#upset_R_df <- transformed_df %>% 
#  select(gene_ID, Log2_TPM, K36me3_K4me1ation, K56ac_K4me3iation) %>% 
#  unique()
#
#samp2 <- as.data.frame(upset_R_df)[,-1]
#
#
##On-the-fly indicator function for use in formulae
#leaf_comb = make_comb_mat(samp2, remove_complement_set =TRUE)
#comb_elements_leaf = lapply(comb_name(leaf_comb), function(nm) extract_comb(leaf_comb, nm))
#View(comb_elements_leaf)
#
##Generate  the list to used in the Annotation Density plots on either the top or bottom of the values
#TPM_val = lapply(comb_elements_leaf, function(ind) samp2$Log2_TPM[ind])
#
#
##Putting the Expression information on the BOTTOM of the Plot
#final <- UpSet(leaf_comb, 
#               top_annotation = upset_top_annotation(leaf_comb, 
#                                                     annotation_name_rot = 90,
#                                                     annotation_name_side = "left",
#                                                     axis_param = list(side = "left")),
#               right_annotation = NULL)  %v% 
#  HeatmapAnnotation("Log2(TPM +1)" = anno_density(TPM_val, type = "violin"), annotation_name_rot = 90, annotation_name_side = "left")
#
#final
#dev.off()
#
#od = column_order(final)
#cs = comb_size(leaf_comb)
#decorate_annotation("Intersection\nsize", {
#  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
#            default.units = "native", just = "bottom", gp = gpar(fontsize = 8))
#})
#
#
#
#
#
##Putting the Expression information on the top of the Plot - possibly prettier in the long run, but is a bit nasty to look at imo
#
#top_ha = HeatmapAnnotation(
#  Log2TPM = anno_density(TPM_val, type = "violin"),
#  "Intersection\nsize" = anno_barplot(comb_size(leaf_comb), 
#                                      border = FALSE, 
#                                      gp = gpar(fill = "black"), 
#                                      height = unit(2, "cm")),
#  gap = unit(2, "mm"), annotation_name_side = "left", annotation_name_rot = 0)
## the same for using m2 or m3
#
#UpSet(leaf_comb, top_annotation = top_ha)