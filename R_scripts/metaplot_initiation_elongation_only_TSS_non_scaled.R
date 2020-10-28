library(tidyverse)
library(EnrichedHeatmap)
library(cowplot)
library("grid")
library("ggplotify")
library("here")
library(circlize)

setwd("/Users/feilab/Projects/03.ncRNA_project/02.Analysis/lncRNA_copy_files/2020-01-06_metaplots_lncRNA_RNAOnly_vs_Myset")




H3K36me3_colors <- "#CB2026"
H3K4me1_colors <-  "#642165"
H3K4me3_colors <- "#B89BC9"
H3K56ac_colors <- "#99CA3C"
ATAC_colors <- "#FFCD05"
H2A_colors <-  "#66c2a5"
H3K27me3_colors <- "#234bb8"



read_gunzipped_file <- function(zipped_file){
  
  bw_zipped_file <- read_delim(zipped_file, delim='\t', col_names = FALSE, skip = 1)
  
  final_return <- bw_zipped_file %>% 
    unite(combined_name, X1,X2,X3,X4,X5,X6, sep='_') %>% 
    column_to_rownames(var="combined_name")
  
  return_matrix <- as.matrix(final_return)
  
  return(return_matrix)
  
  
}

generate_tss_heatmap_scaled_region <- function(TSS_file, mark_name) {
  
  #TSS_file <- TSS_file[, -c(151:200)] # delete columns 5 through 7
  
  TSS <- as.normalizedMatrix(TSS_file, 
                             k_upstream = 100, 
                             k_downstream = 100, 
                             k_target = 0,
                             extend = c(2000, 2000), 
                             signal_name = mark_name, 
                             target_name = c("TTS"),
                             keep = c(0,.99), smooth = FALSE)
  
  
  
  return(TSS)
}

#Generate the plots for scaled region matricies 
generate_complex_heatmaps_chip_scaled_region <- function(TSS_matrix, mark_name, row_split_array, color_hex) {
  
  #If row title is NOT missing - pass on it and don't add it 
  
  TSS_matrix <- TSS_matrix
  
  
  #Scale colors across matricies the same
  #common_min <-  min(TSS_matrix)
  common_max <-  max(TSS_matrix)
  col_fun <- circlize::colorRamp2(c(0, common_max), c("white", color_hex))
  
  #Get the same max value for each 
  #y_max_val <- round(max(colMeans(TES_matrix, na.rm = TRUE), colMeans(TSS_matrix, na.rm = TRUE))) + 2 
  
  axis_name = c("-2kb","TSS", "2kb")
  final_graph <-  EnrichedHeatmap(TSS_matrix,  column_title =  mark_name, show_heatmap_legend = FALSE, use_raster = TRUE,
                                  row_split = row_split_array, cluster_rows = FALSE, 
                                  column_title_gp = gpar(fontsize = 25, fontface = "bold"), 
                                  col = col_fun, axis_name_rot = 90, row_title_rot = 90, axis_name_gp = gpar(fontsize = 18), axis_name = axis_name,pos_line_gp = gpar(lty = 3),
                                  top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = c(H3K27me3_colors, "#888A8D"))))) 
  
  
  
  
  
  return(final_graph)
  
  
}
normalize_to_input <- function(mark_name_1, input_mark){
  normalized_to_input <- (mark_name_1 + 1)-(input_mark +1)
  return(normalized_to_input)
}
combine_identical_marks_on_off <- function(genes_on, genes_off) {
  
  #genes_1 <- map_df(genes) 
  #lncRNA_1 <- map_df(lncRNA)
  #unannotated_1 <- map_df(unannotated)
  
  bind_TSS_rows <- bind_rows(
    as.data.frame(genes_on),
    as.data.frame(genes_off))
  
  final_return <- as.matrix(bind_TSS_rows)
  
  return(final_return)
}

take_subsample_multi_data_frame <- function(number_desired, file_1, file_2, file_3, file_4, file_5, file_6, file_7) {
  file_1_sub_sample <- file_1[sample(nrow(file_1),size=number_desired,replace=FALSE),]
  
  take_subsampled_named <- rownames(file_1_sub_sample)
  file_2_subsampled <- subset(file_2, rownames(file_2) %in% take_subsampled_named)
  file_3_subsampled <- subset(file_3, rownames(file_3) %in% take_subsampled_named)
  file_4_subsampled <- subset(file_4, rownames(file_4) %in% take_subsampled_named)
  file_5_subsampled <- subset(file_5, rownames(file_5) %in% take_subsampled_named)
  file_6_subsampled <- subset(file_6, rownames(file_6) %in% take_subsampled_named)
  file_7_subsampled <- subset(file_7, rownames(file_7) %in% take_subsampled_named)
  
  
  
  return(list(file_1_sub_sample,file_2_subsampled,file_3_subsampled,file_4_subsampled,file_5_subsampled,file_6_subsampled,file_7_subsampled))
}



#Call sub functions depending on mark_type
generate_final_plot_scaled <- function(file_1, mark_name, color_hex, row_split_array) {
  Chip_heatmap <- generate_tss_heatmap_scaled_region(file_1, mark_name)
  Chip_plot_heatmap <- generate_complex_heatmaps_chip_scaled_region(Chip_heatmap, mark_name, row_split_array, color_hex)
  return(Chip_plot_heatmap)
  }







# leaf --------------------------------------------------------------------
leaf_initiation_only_leaf_H2AZ_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_initiation_only_leaf_H2AZ_TSS.gz")
leaf_initiation_only_leaf_H3K27me3_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_initiation_only_leaf_H3K27me3_TSS.gz")
leaf_initiation_only_leaf_H3K36me3_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_initiation_only_leaf_H3K36me3_TSS.gz")
leaf_initiation_only_leaf_H3K4me1_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_initiation_only_leaf_H3K4me1_TSS.gz")
leaf_initiation_only_leaf_H3K4me3_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_initiation_only_leaf_H3K4me3_TSS.gz")
leaf_initiation_only_leaf_H3K56ac_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_initiation_only_leaf_H3K56ac_TSS.gz")
leaf_initiation_only_leaf_input_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_initiation_only_leaf_input_TSS.gz")





bw_zipped_file <- read_delim("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_initiation_only_leaf_input_TSS.gz", delim='\t', col_names = FALSE, skip = 1)



# Load Control Class ------------------------------------------------------
leaf_on_genes_leaf_H2AZ <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_on_genes_leaf_H2AZ_TSS.gz")
leaf_on_genes_leaf_H3K27me3 <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_on_genes_leaf_H3K27me3_TSS.gz")
leaf_on_genes_leaf_H3K36me3 <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_on_genes_leaf_H3K36me3_TSS.gz")
leaf_on_genes_leaf_H3K4me1 <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_on_genes_leaf_H3K4me1_TSS.gz")
leaf_on_genes_leaf_H3K4me3 <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_on_genes_leaf_H3K4me3_TSS.gz")
leaf_on_genes_leaf_H3K56ac <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_on_genes_leaf_H3K56ac_TSS.gz")
leaf_on_genes_leaf_input <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_on_genes_leaf_input_TSS.gz")



take_number <- nrow(leaf_initiation_only_leaf_H3K36me3_TSS)
take_sub_samples <- take_subsample_multi_data_frame(take_number, leaf_on_genes_leaf_H2AZ, leaf_on_genes_leaf_H3K27me3, leaf_on_genes_leaf_H3K36me3, leaf_on_genes_leaf_H3K4me1, leaf_on_genes_leaf_H3K4me3, leaf_on_genes_leaf_H3K56ac, leaf_on_genes_leaf_input)

leaf_on_genes_leaf_H2AZ_sub_sampled <- take_sub_samples[1]
leaf_on_genes_leaf_H3K27me3_sub_sampled <- take_sub_samples[2]
leaf_on_genes_leaf_H3K36me3_sub_sampled <- take_sub_samples[3]
leaf_on_genes_leaf_H3K4me1_sub_sampled <- take_sub_samples[4]
leaf_on_genes_leaf_H3K4me3_sub_sampled <- take_sub_samples[5]
leaf_on_genes_leaf_H3K56ac_sub_sampled <- take_sub_samples[6]
leaf_on_genes_leaf_inp_sub_sampled <- take_sub_samples[7]



leaf_initiation_only <- nrow(leaf_initiation_only_leaf_H2AZ_TSS)
values <- c( "Leaf Initiation Only", "Control Expressed Genes")
final_split <- rep(values, times = c(leaf_initiation_only, leaf_initiation_only))


leaf_initiation_only_control_merged_H2AZ <- combine_identical_marks_on_off(leaf_initiation_only_leaf_H2AZ_TSS,leaf_on_genes_leaf_H2AZ_sub_sampled)
leaf_initiation_only_control_merged_H3K27me3 <- combine_identical_marks_on_off(leaf_initiation_only_leaf_H3K27me3_TSS,leaf_on_genes_leaf_H3K27me3_sub_sampled)
leaf_initiation_only_control_merged_H3K36me3 <- combine_identical_marks_on_off(leaf_initiation_only_leaf_H3K36me3_TSS,leaf_on_genes_leaf_H3K36me3_sub_sampled)
leaf_initiation_only_control_merged_H3K4me1 <- combine_identical_marks_on_off(leaf_initiation_only_leaf_H3K4me1_TSS,leaf_on_genes_leaf_H3K4me1_sub_sampled)
leaf_initiation_only_control_merged_H3K4me3 <- combine_identical_marks_on_off(leaf_initiation_only_leaf_H3K4me3_TSS,leaf_on_genes_leaf_H3K4me3_sub_sampled)
leaf_initiation_only_control_merged_H3K56ac <- combine_identical_marks_on_off(leaf_initiation_only_leaf_H3K56ac_TSS,leaf_on_genes_leaf_H3K56ac_sub_sampled)
leaf_initiation_only_control_merged_input <- combine_identical_marks_on_off(leaf_initiation_only_leaf_input_TSS,leaf_on_genes_leaf_inp_sub_sampled)



leaf_initiation_only_leaf_H2AZ_matrix_normalized <- normalize_to_input(leaf_initiation_only_control_merged_H2AZ, leaf_initiation_only_control_merged_input)
leaf_initiation_only_leaf_H3K27me3_matrix_normalized <- normalize_to_input(leaf_initiation_only_control_merged_H3K27me3, leaf_initiation_only_control_merged_input)
leaf_initiation_only_leaf_H3K36me3_matrix_normalized <- normalize_to_input(leaf_initiation_only_control_merged_H3K36me3, leaf_initiation_only_control_merged_input)
leaf_initiation_only_leaf_H3K4me1_matrix_normalized <- normalize_to_input(leaf_initiation_only_control_merged_H3K4me1, leaf_initiation_only_control_merged_input)
leaf_initiation_only_leaf_H3K4me3_matrix_normalized <- normalize_to_input(leaf_initiation_only_control_merged_H3K4me3, leaf_initiation_only_control_merged_input)
leaf_initiation_only_leaf_H3K56ac_matrix_normalized <- normalize_to_input(leaf_initiation_only_control_merged_H3K56ac, leaf_initiation_only_control_merged_input)


dim(leaf_initiation_only_leaf_H2AZ_matrix_normalized)
length(final_split)

leaf_initiation_only_leaf_H2AZ_plots <- generate_final_plot_scaled(leaf_initiation_only_leaf_H2AZ_matrix_normalized, "H2AZ", H2A_colors, final_split)
leaf_initiation_only_leaf_H3K27me3_plots <- generate_final_plot_scaled(leaf_initiation_only_leaf_H3K27me3_matrix_normalized, "H3K27me3", H3K27me3_colors, final_split)
leaf_initiation_only_leaf_H3K36me3_plots <- generate_final_plot_scaled(leaf_initiation_only_leaf_H3K36me3_matrix_normalized, "H3K36me3", H3K36me3_colors, final_split)
leaf_initiation_only_leaf_H3K4me1_plots <- generate_final_plot_scaled(leaf_initiation_only_leaf_H3K4me1_matrix_normalized, "H3K4me1", H3K4me1_colors, final_split)
leaf_initiation_only_leaf_H3K4me3_plots <- generate_final_plot_scaled(leaf_initiation_only_leaf_H3K4me3_matrix_normalized, "H3K4me3", H3K4me3_colors, final_split)
leaf_initiation_only_leaf_H3K56ac_plots <- generate_final_plot_scaled(leaf_initiation_only_leaf_H3K56ac_matrix_normalized, "H3K56ac", H3K56ac_colors, final_split)



leaf_initiation_only_leaf_plot <- (leaf_initiation_only_leaf_H3K4me3_plots + 
                                     leaf_initiation_only_leaf_H3K56ac_plots + 
                                     leaf_initiation_only_leaf_H3K36me3_plots + 
                                     leaf_initiation_only_leaf_H3K4me1_plots + 
                                     leaf_initiation_only_leaf_H2AZ_plots + 
                                     leaf_initiation_only_leaf_H3K27me3_plots)


units_apart <- c(.7)
draw_units <- rep(units_apart, 6)

leaf_initiation_only_leaf_plot_draw = grid.grabExpr(draw(leaf_initiation_only_leaf_plot, ht_gap = unit(draw_units, "cm")))
ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/", filename="leaf_initiation_only_leaf_TSS_non_scaled.pdf", plot=leaf_initiation_only_leaf_plot_draw, width = 6, height = 4, units = "in")
system("open /Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots")







# root --------------------------------------------------------------------
root_initiation_only_root_H2AZ_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_initiation_only_root_H2AZ_TSS.gz")
root_initiation_only_root_H3K27me3_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_initiation_only_root_H3K27me3_TSS.gz")
root_initiation_only_root_H3K36me3_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_initiation_only_root_H3K36me3_TSS.gz")
root_initiation_only_root_H3K4me1_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_initiation_only_root_H3K4me1_TSS.gz")
root_initiation_only_root_H3K4me3_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_initiation_only_root_H3K4me3_TSS.gz")
root_initiation_only_root_H3K56ac_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_initiation_only_root_H3K56ac_TSS.gz")
root_initiation_only_root_input_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_initiation_only_root_input_TSS.gz")


nrow(root_on_genes_root_H2AZ)


# Load Control Class ------------------------------------------------------
root_on_genes_root_H2AZ <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_on_genes_root_H2AZ_TSS.gz")
root_on_genes_root_H3K27me3 <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_on_genes_root_H3K27me3_TSS.gz")
root_on_genes_root_H3K36me3 <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_on_genes_root_H3K36me3_TSS.gz")
root_on_genes_root_H3K4me1 <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_on_genes_root_H3K4me1_TSS.gz")
root_on_genes_root_H3K4me3 <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_on_genes_root_H3K4me3_TSS.gz")
root_on_genes_root_H3K56ac <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_on_genes_root_H3K56ac_TSS.gz")
root_on_genes_root_input <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_on_genes_root_input_TSS.gz")



take_number <- nrow(root_initiation_only_root_H3K36me3_TSS)
take_sub_samples <- take_subsample_multi_data_frame(take_number, root_on_genes_root_H2AZ, root_on_genes_root_H3K27me3, root_on_genes_root_H3K36me3, root_on_genes_root_H3K4me1, root_on_genes_root_H3K4me3, root_on_genes_root_H3K56ac, root_on_genes_root_input)

root_on_genes_root_H2AZ_sub_sampled <- take_sub_samples[1]
root_on_genes_root_H3K27me3_sub_sampled <- take_sub_samples[2]
root_on_genes_root_H3K36me3_sub_sampled <- take_sub_samples[3]
root_on_genes_root_H3K4me1_sub_sampled <- take_sub_samples[4]
root_on_genes_root_H3K4me3_sub_sampled <- take_sub_samples[5]
root_on_genes_root_H3K56ac_sub_sampled <- take_sub_samples[6]
root_on_genes_root_inp_sub_sampled <- take_sub_samples[7]



root_initiation_only <- nrow(root_initiation_only_root_H2AZ_TSS)
values <- c( "root Initiation Only", "Control Expressed Genes")
final_split <- rep(values, times = c(root_initiation_only, root_initiation_only))


root_initiation_only_control_merged_H2AZ <- combine_identical_marks_on_off(root_initiation_only_root_H2AZ_TSS,root_on_genes_root_H2AZ_sub_sampled)
root_initiation_only_control_merged_H3K27me3 <- combine_identical_marks_on_off(root_initiation_only_root_H3K27me3_TSS,root_on_genes_root_H3K27me3_sub_sampled)
root_initiation_only_control_merged_H3K36me3 <- combine_identical_marks_on_off(root_initiation_only_root_H3K36me3_TSS,root_on_genes_root_H3K36me3_sub_sampled)
root_initiation_only_control_merged_H3K4me1 <- combine_identical_marks_on_off(root_initiation_only_root_H3K4me1_TSS,root_on_genes_root_H3K4me1_sub_sampled)
root_initiation_only_control_merged_H3K4me3 <- combine_identical_marks_on_off(root_initiation_only_root_H3K4me3_TSS,root_on_genes_root_H3K4me3_sub_sampled)
root_initiation_only_control_merged_H3K56ac <- combine_identical_marks_on_off(root_initiation_only_root_H3K56ac_TSS,root_on_genes_root_H3K56ac_sub_sampled)
root_initiation_only_control_merged_input <- combine_identical_marks_on_off(root_initiation_only_root_input_TSS,root_on_genes_root_inp_sub_sampled)



root_initiation_only_root_H2AZ_matrix_normalized <- normalize_to_input(root_initiation_only_control_merged_H2AZ, root_initiation_only_control_merged_input)
root_initiation_only_root_H3K27me3_matrix_normalized <- normalize_to_input(root_initiation_only_control_merged_H3K27me3, root_initiation_only_control_merged_input)
root_initiation_only_root_H3K36me3_matrix_normalized <- normalize_to_input(root_initiation_only_control_merged_H3K36me3, root_initiation_only_control_merged_input)
root_initiation_only_root_H3K4me1_matrix_normalized <- normalize_to_input(root_initiation_only_control_merged_H3K4me1, root_initiation_only_control_merged_input)
root_initiation_only_root_H3K4me3_matrix_normalized <- normalize_to_input(root_initiation_only_control_merged_H3K4me3, root_initiation_only_control_merged_input)
root_initiation_only_root_H3K56ac_matrix_normalized <- normalize_to_input(root_initiation_only_control_merged_H3K56ac, root_initiation_only_control_merged_input)



root_initiation_only_root_H2AZ_plots <- generate_final_plot_scaled(root_initiation_only_root_H2AZ_matrix_normalized, "H2AZ", H2A_colors, final_split)
root_initiation_only_root_H3K27me3_plots <- generate_final_plot_scaled(root_initiation_only_root_H3K27me3_matrix_normalized, "H3K27me3", H3K27me3_colors, final_split)
root_initiation_only_root_H3K36me3_plots <- generate_final_plot_scaled(root_initiation_only_root_H3K36me3_matrix_normalized, "H3K36me3", H3K36me3_colors, final_split)
root_initiation_only_root_H3K4me1_plots <- generate_final_plot_scaled(root_initiation_only_root_H3K4me1_matrix_normalized, "H3K4me1", H3K4me1_colors, final_split)
root_initiation_only_root_H3K4me3_plots <- generate_final_plot_scaled(root_initiation_only_root_H3K4me3_matrix_normalized, "H3K4me3", H3K4me3_colors, final_split)
root_initiation_only_root_H3K56ac_plots <- generate_final_plot_scaled(root_initiation_only_root_H3K56ac_matrix_normalized, "H3K56ac", H3K56ac_colors, final_split)



root_initiation_only_root_plot <- (root_initiation_only_root_H3K4me3_plots + 
                                     root_initiation_only_root_H3K56ac_plots + 
                                     root_initiation_only_root_H3K36me3_plots + 
                                     root_initiation_only_root_H3K4me1_plots + 
                                     root_initiation_only_root_H2AZ_plots + 
                                     root_initiation_only_root_H3K27me3_plots)

root_initiation_only_root_plot
units_apart <- c(.7)
draw_units <- rep(units_apart, 6)

root_initiation_only_root_plot_draw = grid.grabExpr(draw(root_initiation_only_root_plot, ht_gap = unit(draw_units, "cm")))
ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/", filename="root_initiation_only_root_TSS_non_scaled.png", plot=root_initiation_only_root_plot_draw, width = 15, height = 18, units = "in")
system("open /Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/root_initiation_only_root_TSS_non_scaled.png")




# ear --------------------------------------------------------------------
ear_initiation_only_ear_H2AZ_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_initiation_only_ear_H2AZ_TSS.gz")
ear_initiation_only_ear_H3K27me3_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_initiation_only_ear_H3K27me3_TSS.gz")
ear_initiation_only_ear_H3K36me3_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_initiation_only_ear_H3K36me3_TSS.gz")
ear_initiation_only_ear_H3K4me1_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_initiation_only_ear_H3K4me1_TSS.gz")
ear_initiation_only_ear_H3K4me3_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_initiation_only_ear_H3K4me3_TSS.gz")
ear_initiation_only_ear_H3K56ac_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_initiation_only_ear_H3K56ac_TSS.gz")
ear_initiation_only_ear_input_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_initiation_only_ear_input_TSS.gz")


nrow(ear_on_genes_ear_H2AZ)


# Load Control Class ------------------------------------------------------
ear_on_genes_ear_H2AZ <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_on_genes_ear_H2AZ_TSS.gz")
ear_on_genes_ear_H3K27me3 <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_on_genes_ear_H3K27me3_TSS.gz")
ear_on_genes_ear_H3K36me3 <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_on_genes_ear_H3K36me3_TSS.gz")
ear_on_genes_ear_H3K4me1 <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_on_genes_ear_H3K4me1_TSS.gz")
ear_on_genes_ear_H3K4me3 <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_on_genes_ear_H3K4me3_TSS.gz")
ear_on_genes_ear_H3K56ac <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_on_genes_ear_H3K56ac_TSS.gz")
ear_on_genes_ear_input <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_on_genes_ear_input_TSS.gz")



take_number <- nrow(ear_initiation_only_ear_H3K36me3_TSS)
take_sub_samples <- take_subsample_multi_data_frame(take_number, ear_on_genes_ear_H2AZ, ear_on_genes_ear_H3K27me3, ear_on_genes_ear_H3K36me3, ear_on_genes_ear_H3K4me1, ear_on_genes_ear_H3K4me3, ear_on_genes_ear_H3K56ac, ear_on_genes_ear_input)

ear_on_genes_ear_H2AZ_sub_sampled <- take_sub_samples[1]
ear_on_genes_ear_H3K27me3_sub_sampled <- take_sub_samples[2]
ear_on_genes_ear_H3K36me3_sub_sampled <- take_sub_samples[3]
ear_on_genes_ear_H3K4me1_sub_sampled <- take_sub_samples[4]
ear_on_genes_ear_H3K4me3_sub_sampled <- take_sub_samples[5]
ear_on_genes_ear_H3K56ac_sub_sampled <- take_sub_samples[6]
ear_on_genes_ear_inp_sub_sampled <- take_sub_samples[7]



ear_initiation_only <- nrow(ear_initiation_only_ear_H2AZ_TSS)

values <- c( "ear Initiation Only", "Control Expressed Genes")
final_split <- rep(values, times = c(ear_initiation_only, ear_initiation_only))


ear_initiation_only_control_merged_H2AZ <- combine_identical_marks_on_off(ear_initiation_only_ear_H2AZ_TSS,ear_on_genes_ear_H2AZ_sub_sampled)
ear_initiation_only_control_merged_H3K27me3 <- combine_identical_marks_on_off(ear_initiation_only_ear_H3K27me3_TSS,ear_on_genes_ear_H3K27me3_sub_sampled)
ear_initiation_only_control_merged_H3K36me3 <- combine_identical_marks_on_off(ear_initiation_only_ear_H3K36me3_TSS,ear_on_genes_ear_H3K36me3_sub_sampled)
ear_initiation_only_control_merged_H3K4me1 <- combine_identical_marks_on_off(ear_initiation_only_ear_H3K4me1_TSS,ear_on_genes_ear_H3K4me1_sub_sampled)
ear_initiation_only_control_merged_H3K4me3 <- combine_identical_marks_on_off(ear_initiation_only_ear_H3K4me3_TSS,ear_on_genes_ear_H3K4me3_sub_sampled)
ear_initiation_only_control_merged_H3K56ac <- combine_identical_marks_on_off(ear_initiation_only_ear_H3K56ac_TSS,ear_on_genes_ear_H3K56ac_sub_sampled)
ear_initiation_only_control_merged_input <- combine_identical_marks_on_off(ear_initiation_only_ear_input_TSS,ear_on_genes_ear_inp_sub_sampled)



ear_initiation_only_ear_H2AZ_matrix_normalized <- normalize_to_input(ear_initiation_only_control_merged_H2AZ, ear_initiation_only_control_merged_input)
ear_initiation_only_ear_H3K27me3_matrix_normalized <- normalize_to_input(ear_initiation_only_control_merged_H3K27me3, ear_initiation_only_control_merged_input)
ear_initiation_only_ear_H3K36me3_matrix_normalized <- normalize_to_input(ear_initiation_only_control_merged_H3K36me3, ear_initiation_only_control_merged_input)
ear_initiation_only_ear_H3K4me1_matrix_normalized <- normalize_to_input(ear_initiation_only_control_merged_H3K4me1, ear_initiation_only_control_merged_input)
ear_initiation_only_ear_H3K4me3_matrix_normalized <- normalize_to_input(ear_initiation_only_control_merged_H3K4me3, ear_initiation_only_control_merged_input)
ear_initiation_only_ear_H3K56ac_matrix_normalized <- normalize_to_input(ear_initiation_only_control_merged_H3K56ac, ear_initiation_only_control_merged_input)



ear_initiation_only_ear_H2AZ_plots <- generate_final_plot_scaled(ear_initiation_only_ear_H2AZ_matrix_normalized, "H2AZ", H2A_colors, final_split)
ear_initiation_only_ear_H3K27me3_plots <- generate_final_plot_scaled(ear_initiation_only_ear_H3K27me3_matrix_normalized, "H3K27me3", H3K27me3_colors, final_split)
ear_initiation_only_ear_H3K36me3_plots <- generate_final_plot_scaled(ear_initiation_only_ear_H3K36me3_matrix_normalized, "H3K36me3", H3K36me3_colors, final_split)
ear_initiation_only_ear_H3K4me1_plots <- generate_final_plot_scaled(ear_initiation_only_ear_H3K4me1_matrix_normalized, "H3K4me1", H3K4me1_colors, final_split)
ear_initiation_only_ear_H3K4me3_plots <- generate_final_plot_scaled(ear_initiation_only_ear_H3K4me3_matrix_normalized, "H3K4me3", H3K4me3_colors, final_split)
ear_initiation_only_ear_H3K56ac_plots <- generate_final_plot_scaled(ear_initiation_only_ear_H3K56ac_matrix_normalized, "H3K56ac", H3K56ac_colors, final_split)



ear_initiation_only_ear_plot <- (ear_initiation_only_ear_H3K4me3_plots + 
                                   ear_initiation_only_ear_H3K56ac_plots + 
                                   ear_initiation_only_ear_H3K36me3_plots + 
                                   ear_initiation_only_ear_H3K4me1_plots + 
                                   ear_initiation_only_ear_H2AZ_plots + 
                                   ear_initiation_only_ear_H3K27me3_plots)


units_apart <- c(.7)
draw_units <- rep(units_apart, 6)

ear_initiation_only_ear_plot_draw = grid.grabExpr(draw(ear_initiation_only_ear_plot, ht_gap = unit(draw_units, "cm")))
ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/", filename="ear_initiation_only_ear_TSS_non_scaled.png", plot=ear_initiation_only_ear_plot_draw, width = 15, height = 18, units = "in")
system("open /Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/ear_initiation_only_ear_TSS_non_scaled.png")


######################s######################s###################################################s#######################################################S
######################s######################s###################################################s#######################################################S
######################s######################s###################################################s#######################################################S
######################s######################s###################################################s#######################################################S
######################s######################s###################################################s#######################################################S
######################s######################s###################################################s#######################################################S
#ELONGATION ONLY
######################s######################s###################################################s#######################################################S
######################s######################s###################################################s#######################################################S
######################s######################s###################################################s#######################################################S
######################s######################s###################################################s#######################################################S
######################s######################s###################################################s#######################################################S
######################s######################s###################################################s#######################################################S
######################s######################s###################################################s#######################################################S
######################s######################s###################################################s#######################################################S
######################s######################s###################################################s#######################################################S
######################s######################s###################################################s#######################################################S
######################s######################s###################################################s#######################################################S
######################s######################s###################################################s#######################################################S






# leaf --------------------------------------------------------------------
leaf_elongation_only_leaf_H2AZ_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_elongation_only_leaf_H2AZ_TSS.gz")
leaf_elongation_only_leaf_H3K27me3_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_elongation_only_leaf_H3K27me3_TSS.gz")
leaf_elongation_only_leaf_H3K36me3_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_elongation_only_leaf_H3K36me3_TSS.gz")
leaf_elongation_only_leaf_H3K4me1_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_elongation_only_leaf_H3K4me1_TSS.gz")
leaf_elongation_only_leaf_H3K4me3_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_elongation_only_leaf_H3K4me3_TSS.gz")
leaf_elongation_only_leaf_H3K56ac_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_elongation_only_leaf_H3K56ac_TSS.gz")
leaf_elongation_only_leaf_input_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_elongation_only_leaf_input_TSS.gz")





# Load Control Class ------------------------------------------------------
leaf_on_genes_leaf_H2AZ <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_on_genes_leaf_H2AZ_TSS.gz")
leaf_on_genes_leaf_H3K27me3 <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_on_genes_leaf_H3K27me3_TSS.gz")
leaf_on_genes_leaf_H3K36me3 <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_on_genes_leaf_H3K36me3_TSS.gz")
leaf_on_genes_leaf_H3K4me1 <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_on_genes_leaf_H3K4me1_TSS.gz")
leaf_on_genes_leaf_H3K4me3 <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_on_genes_leaf_H3K4me3_TSS.gz")
leaf_on_genes_leaf_H3K56ac <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_on_genes_leaf_H3K56ac_TSS.gz")
leaf_on_genes_leaf_input <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_on_genes_leaf_input_TSS.gz")



take_number <- nrow(leaf_elongation_only_leaf_H3K36me3_TSS)
take_sub_samples <- take_subsample_multi_data_frame(take_number, leaf_on_genes_leaf_H2AZ, leaf_on_genes_leaf_H3K27me3, leaf_on_genes_leaf_H3K36me3, leaf_on_genes_leaf_H3K4me1, leaf_on_genes_leaf_H3K4me3, leaf_on_genes_leaf_H3K56ac, leaf_on_genes_leaf_input)

leaf_on_genes_leaf_H2AZ_sub_sampled <- take_sub_samples[1]
leaf_on_genes_leaf_H3K27me3_sub_sampled <- take_sub_samples[2]
leaf_on_genes_leaf_H3K36me3_sub_sampled <- take_sub_samples[3]
leaf_on_genes_leaf_H3K4me1_sub_sampled <- take_sub_samples[4]
leaf_on_genes_leaf_H3K4me3_sub_sampled <- take_sub_samples[5]
leaf_on_genes_leaf_H3K56ac_sub_sampled <- take_sub_samples[6]
leaf_on_genes_leaf_inp_sub_sampled <- take_sub_samples[7]



leaf_elongation_only <- nrow(leaf_elongation_only_leaf_H2AZ_TSS)

values <- c( "Leaf elongation Only", "Control Expressed Genes")
final_split <- rep(values, times = c(leaf_elongation_only, leaf_elongation_only))


leaf_elongation_only_control_merged_H2AZ <- combine_identical_marks_on_off(leaf_elongation_only_leaf_H2AZ_TSS,leaf_on_genes_leaf_H2AZ_sub_sampled)
leaf_elongation_only_control_merged_H3K27me3 <- combine_identical_marks_on_off(leaf_elongation_only_leaf_H3K27me3_TSS,leaf_on_genes_leaf_H3K27me3_sub_sampled)
leaf_elongation_only_control_merged_H3K36me3 <- combine_identical_marks_on_off(leaf_elongation_only_leaf_H3K36me3_TSS,leaf_on_genes_leaf_H3K36me3_sub_sampled)
leaf_elongation_only_control_merged_H3K4me1 <- combine_identical_marks_on_off(leaf_elongation_only_leaf_H3K4me1_TSS,leaf_on_genes_leaf_H3K4me1_sub_sampled)
leaf_elongation_only_control_merged_H3K4me3 <- combine_identical_marks_on_off(leaf_elongation_only_leaf_H3K4me3_TSS,leaf_on_genes_leaf_H3K4me3_sub_sampled)
leaf_elongation_only_control_merged_H3K56ac <- combine_identical_marks_on_off(leaf_elongation_only_leaf_H3K56ac_TSS,leaf_on_genes_leaf_H3K56ac_sub_sampled)
leaf_elongation_only_control_merged_input <- combine_identical_marks_on_off(leaf_elongation_only_leaf_input_TSS,leaf_on_genes_leaf_inp_sub_sampled)



leaf_elongation_only_leaf_H2AZ_matrix_normalized <- normalize_to_input(leaf_elongation_only_control_merged_H2AZ, leaf_elongation_only_control_merged_input)
leaf_elongation_only_leaf_H3K27me3_matrix_normalized <- normalize_to_input(leaf_elongation_only_control_merged_H3K27me3, leaf_elongation_only_control_merged_input)
leaf_elongation_only_leaf_H3K36me3_matrix_normalized <- normalize_to_input(leaf_elongation_only_control_merged_H3K36me3, leaf_elongation_only_control_merged_input)
leaf_elongation_only_leaf_H3K4me1_matrix_normalized <- normalize_to_input(leaf_elongation_only_control_merged_H3K4me1, leaf_elongation_only_control_merged_input)
leaf_elongation_only_leaf_H3K4me3_matrix_normalized <- normalize_to_input(leaf_elongation_only_control_merged_H3K4me3, leaf_elongation_only_control_merged_input)
leaf_elongation_only_leaf_H3K56ac_matrix_normalized <- normalize_to_input(leaf_elongation_only_control_merged_H3K56ac, leaf_elongation_only_control_merged_input)




leaf_elongation_only_leaf_H2AZ_plots <- generate_final_plot_scaled(leaf_elongation_only_leaf_H2AZ_matrix_normalized, "H2AZ", H2A_colors, final_split)
leaf_elongation_only_leaf_H3K27me3_plots <- generate_final_plot_scaled(leaf_elongation_only_leaf_H3K27me3_matrix_normalized, "H3K27me3", H3K27me3_colors, final_split)
leaf_elongation_only_leaf_H3K36me3_plots <- generate_final_plot_scaled(leaf_elongation_only_leaf_H3K36me3_matrix_normalized, "H3K36me3", H3K36me3_colors, final_split)
leaf_elongation_only_leaf_H3K4me1_plots <- generate_final_plot_scaled(leaf_elongation_only_leaf_H3K4me1_matrix_normalized, "H3K4me1", H3K4me1_colors, final_split)
leaf_elongation_only_leaf_H3K4me3_plots <- generate_final_plot_scaled(leaf_elongation_only_leaf_H3K4me3_matrix_normalized, "H3K4me3", H3K4me3_colors, final_split)
leaf_elongation_only_leaf_H3K56ac_plots <- generate_final_plot_scaled(leaf_elongation_only_leaf_H3K56ac_matrix_normalized, "H3K56ac", H3K56ac_colors, final_split)



leaf_elongation_only_leaf_plot <- (leaf_elongation_only_leaf_H3K4me3_plots + 
                                     leaf_elongation_only_leaf_H3K56ac_plots + 
                                     leaf_elongation_only_leaf_H3K36me3_plots + 
                                     leaf_elongation_only_leaf_H3K4me1_plots + 
                                     leaf_elongation_only_leaf_H2AZ_plots + 
                                     leaf_elongation_only_leaf_H3K27me3_plots)


units_apart <- c(.7)
draw_units <- rep(units_apart, 6)

leaf_elongation_only_leaf_plot_draw = grid.grabExpr(draw(leaf_elongation_only_leaf_plot, ht_gap = unit(draw_units, "cm")))
ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/", filename="leaf_elongation_only_leaf_TSS_non_scaled.pdf", plot=leaf_elongation_only_leaf_plot_draw, width = 6, height = 6, units = "in")
system("open /Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/leaf_elongation_only_leaf_TSS_non_scaled.pdf")







# root --------------------------------------------------------------------
root_elongation_only_root_H2AZ_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_elongation_only_root_H2AZ_TSS.gz")
root_elongation_only_root_H3K27me3_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_elongation_only_root_H3K27me3_TSS.gz")
root_elongation_only_root_H3K36me3_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_elongation_only_root_H3K36me3_TSS.gz")
root_elongation_only_root_H3K4me1_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_elongation_only_root_H3K4me1_TSS.gz")
root_elongation_only_root_H3K4me3_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_elongation_only_root_H3K4me3_TSS.gz")
root_elongation_only_root_H3K56ac_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_elongation_only_root_H3K56ac_TSS.gz")
root_elongation_only_root_input_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_elongation_only_root_input_TSS.gz")


nrow(root_on_genes_root_H2AZ)


# Load Control Class ------------------------------------------------------
root_on_genes_root_H2AZ <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_on_genes_root_H2AZ_TSS.gz")
root_on_genes_root_H3K27me3 <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_on_genes_root_H3K27me3_TSS.gz")
root_on_genes_root_H3K36me3 <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_on_genes_root_H3K36me3_TSS.gz")
root_on_genes_root_H3K4me1 <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_on_genes_root_H3K4me1_TSS.gz")
root_on_genes_root_H3K4me3 <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_on_genes_root_H3K4me3_TSS.gz")
root_on_genes_root_H3K56ac <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_on_genes_root_H3K56ac_TSS.gz")
root_on_genes_root_input <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_on_genes_root_input_TSS.gz")



take_number <- nrow(root_elongation_only_root_H3K36me3_TSS)
take_sub_samples <- take_subsample_multi_data_frame(take_number, root_on_genes_root_H2AZ, root_on_genes_root_H3K27me3, root_on_genes_root_H3K36me3, root_on_genes_root_H3K4me1, root_on_genes_root_H3K4me3, root_on_genes_root_H3K56ac, root_on_genes_root_input)

root_on_genes_root_H2AZ_sub_sampled <- take_sub_samples[1]
root_on_genes_root_H3K27me3_sub_sampled <- take_sub_samples[2]
root_on_genes_root_H3K36me3_sub_sampled <- take_sub_samples[3]
root_on_genes_root_H3K4me1_sub_sampled <- take_sub_samples[4]
root_on_genes_root_H3K4me3_sub_sampled <- take_sub_samples[5]
root_on_genes_root_H3K56ac_sub_sampled <- take_sub_samples[6]
root_on_genes_root_inp_sub_sampled <- take_sub_samples[7]



root_elongation_only <- nrow(root_elongation_only_root_H2AZ_TSS)
values <- c( "root elongation Only", "Control Expressed Genes")
final_split <- rep(values, times = c(root_elongation_only, root_elongation_only))


root_elongation_only_control_merged_H2AZ <- combine_identical_marks_on_off(root_elongation_only_root_H2AZ_TSS,root_on_genes_root_H2AZ_sub_sampled)f
root_elongation_only_control_merged_H3K27me3 <- combine_identical_marks_on_off(root_elongation_only_root_H3K27me3_TSS,root_on_genes_root_H3K27me3_sub_sampled)
root_elongation_only_control_merged_H3K36me3 <- combine_identical_marks_on_off(root_elongation_only_root_H3K36me3_TSS,root_on_genes_root_H3K36me3_sub_sampled)
root_elongation_only_control_merged_H3K4me1 <- combine_identical_marks_on_off(root_elongation_only_root_H3K4me1_TSS,root_on_genes_root_H3K4me1_sub_sampled)
root_elongation_only_control_merged_H3K4me3 <- combine_identical_marks_on_off(root_elongation_only_root_H3K4me3_TSS,root_on_genes_root_H3K4me3_sub_sampled)
root_elongation_only_control_merged_H3K56ac <- combine_identical_marks_on_off(root_elongation_only_root_H3K56ac_TSS,root_on_genes_root_H3K56ac_sub_sampled)
root_elongation_only_control_merged_input <- combine_identical_marks_on_off(root_elongation_only_root_input_TSS,root_on_genes_root_inp_sub_sampled)



root_elongation_only_root_H2AZ_matrix_normalized <- normalize_to_input(root_elongation_only_control_merged_H2AZ, root_elongation_only_control_merged_input)
root_elongation_only_root_H3K27me3_matrix_normalized <- normalize_to_input(root_elongation_only_control_merged_H3K27me3, root_elongation_only_control_merged_input)
root_elongation_only_root_H3K36me3_matrix_normalized <- normalize_to_input(root_elongation_only_control_merged_H3K36me3, root_elongation_only_control_merged_input)
root_elongation_only_root_H3K4me1_matrix_normalized <- normalize_to_input(root_elongation_only_control_merged_H3K4me1, root_elongation_only_control_merged_input)
root_elongation_only_root_H3K4me3_matrix_normalized <- normalize_to_input(root_elongation_only_control_merged_H3K4me3, root_elongation_only_control_merged_input)
root_elongation_only_root_H3K56ac_matrix_normalized <- normalize_to_input(root_elongation_only_control_merged_H3K56ac, root_elongation_only_control_merged_input)



root_elongation_only_root_H2AZ_plots <- generate_final_plot_scaled(root_elongation_only_root_H2AZ_matrix_normalized, "H2AZ", H2A_colors, final_split)
root_elongation_only_root_H3K27me3_plots <- generate_final_plot_scaled(root_elongation_only_root_H3K27me3_matrix_normalized, "H3K27me3", H3K27me3_colors, final_split)
root_elongation_only_root_H3K36me3_plots <- generate_final_plot_scaled(root_elongation_only_root_H3K36me3_matrix_normalized, "H3K36me3", H3K36me3_colors, final_split)
root_elongation_only_root_H3K4me1_plots <- generate_final_plot_scaled(root_elongation_only_root_H3K4me1_matrix_normalized, "H3K4me1", H3K4me1_colors, final_split)
root_elongation_only_root_H3K4me3_plots <- generate_final_plot_scaled(root_elongation_only_root_H3K4me3_matrix_normalized, "H3K4me3", H3K4me3_colors, final_split)
root_elongation_only_root_H3K56ac_plots <- generate_final_plot_scaled(root_elongation_only_root_H3K56ac_matrix_normalized, "H3K56ac", H3K56ac_colors, final_split)



root_elongation_only_root_plot <- (root_elongation_only_root_H3K4me3_plots + 
                                     root_elongation_only_root_H3K56ac_plots + 
                                     root_elongation_only_root_H3K36me3_plots + 
                                     root_elongation_only_root_H3K4me1_plots + 
                                     root_elongation_only_root_H2AZ_plots + 
                                     root_elongation_only_root_H3K27me3_plots)

root_elongation_only_root_plot
units_apart <- c(.7)
draw_units <- rep(units_apart, 6)

root_elongation_only_root_plot_draw = grid.grabExpr(draw(root_elongation_only_root_plot, ht_gap = unit(draw_units, "cm")))
ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/", filename="root_elongation_only_root_TSS_non_scaled.png", plot=root_elongation_only_root_plot_draw, width = 15, height = 18, units = "in")
system("open /Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/root_elongation_only_root_TSS_non_scaled.png")




# ear --------------------------------------------------------------------
ear_elongation_only_ear_H2AZ_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_elongation_only_ear_H2AZ_TSS.gz")
ear_elongation_only_ear_H3K27me3_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_elongation_only_ear_H3K27me3_TSS.gz")
ear_elongation_only_ear_H3K36me3_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_elongation_only_ear_H3K36me3_TSS.gz")
ear_elongation_only_ear_H3K4me1_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_elongation_only_ear_H3K4me1_TSS.gz")
ear_elongation_only_ear_H3K4me3_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_elongation_only_ear_H3K4me3_TSS.gz")
ear_elongation_only_ear_H3K56ac_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_elongation_only_ear_H3K56ac_TSS.gz")
ear_elongation_only_ear_input_TSS <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_elongation_only_ear_input_TSS.gz")


nrow(ear_on_genes_ear_H2AZ)


# Load Control Class ------------------------------------------------------
ear_on_genes_ear_H2AZ <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_on_genes_ear_H2AZ_TSS.gz")
ear_on_genes_ear_H3K27me3 <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_on_genes_ear_H3K27me3_TSS.gz")
ear_on_genes_ear_H3K36me3 <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_on_genes_ear_H3K36me3_TSS.gz")
ear_on_genes_ear_H3K4me1 <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_on_genes_ear_H3K4me1_TSS.gz")
ear_on_genes_ear_H3K4me3 <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_on_genes_ear_H3K4me3_TSS.gz")
ear_on_genes_ear_H3K56ac <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_on_genes_ear_H3K56ac_TSS.gz")
ear_on_genes_ear_input <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_on_genes_ear_input_TSS.gz")



take_number <- nrow(ear_elongation_only_ear_H3K36me3_TSS)
take_sub_samples <- take_subsample_multi_data_frame(take_number, ear_on_genes_ear_H2AZ, ear_on_genes_ear_H3K27me3, ear_on_genes_ear_H3K36me3, ear_on_genes_ear_H3K4me1, ear_on_genes_ear_H3K4me3, ear_on_genes_ear_H3K56ac, ear_on_genes_ear_input)

ear_on_genes_ear_H2AZ_sub_sampled <- take_sub_samples[1]
ear_on_genes_ear_H3K27me3_sub_sampled <- take_sub_samples[2]
ear_on_genes_ear_H3K36me3_sub_sampled <- take_sub_samples[3]
ear_on_genes_ear_H3K4me1_sub_sampled <- take_sub_samples[4]
ear_on_genes_ear_H3K4me3_sub_sampled <- take_sub_samples[5]
ear_on_genes_ear_H3K56ac_sub_sampled <- take_sub_samples[6]
ear_on_genes_ear_inp_sub_sampled <- take_sub_samples[7]



ear_elongation_only <- nrow(ear_elongation_only_ear_H2AZ_TSS)

values <- c( "ear elongation Only", "Control Expressed Genes")
final_split <- rep(values, times = c(ear_elongation_only, ear_elongation_only))


ear_elongation_only_control_merged_H2AZ <- combine_identical_marks_on_off(ear_elongation_only_ear_H2AZ_TSS,ear_on_genes_ear_H2AZ_sub_sampled)
ear_elongation_only_control_merged_H3K27me3 <- combine_identical_marks_on_off(ear_elongation_only_ear_H3K27me3_TSS,ear_on_genes_ear_H3K27me3_sub_sampled)
ear_elongation_only_control_merged_H3K36me3 <- combine_identical_marks_on_off(ear_elongation_only_ear_H3K36me3_TSS,ear_on_genes_ear_H3K36me3_sub_sampled)
ear_elongation_only_control_merged_H3K4me1 <- combine_identical_marks_on_off(ear_elongation_only_ear_H3K4me1_TSS,ear_on_genes_ear_H3K4me1_sub_sampled)
ear_elongation_only_control_merged_H3K4me3 <- combine_identical_marks_on_off(ear_elongation_only_ear_H3K4me3_TSS,ear_on_genes_ear_H3K4me3_sub_sampled)
ear_elongation_only_control_merged_H3K56ac <- combine_identical_marks_on_off(ear_elongation_only_ear_H3K56ac_TSS,ear_on_genes_ear_H3K56ac_sub_sampled)
ear_elongation_only_control_merged_input <- combine_identical_marks_on_off(ear_elongation_only_ear_input_TSS,ear_on_genes_ear_inp_sub_sampled)


ear_elongation_only_ear_H2AZ_matrix_normalized <- normalize_to_input(ear_elongation_only_control_merged_H2AZ, ear_elongation_only_control_merged_input)
ear_elongation_only_ear_H3K27me3_matrix_normalized <- normalize_to_input(ear_elongation_only_control_merged_H3K27me3, ear_elongation_only_control_merged_input)
ear_elongation_only_ear_H3K36me3_matrix_normalized <- normalize_to_input(ear_elongation_only_control_merged_H3K36me3, ear_elongation_only_control_merged_input)
ear_elongation_only_ear_H3K4me1_matrix_normalized <- normalize_to_input(ear_elongation_only_control_merged_H3K4me1, ear_elongation_only_control_merged_input)
ear_elongation_only_ear_H3K4me3_matrix_normalized <- normalize_to_input(ear_elongation_only_control_merged_H3K4me3, ear_elongation_only_control_merged_input)
ear_elongation_only_ear_H3K56ac_matrix_normalized <- normalize_to_input(ear_elongation_only_control_merged_H3K56ac, ear_elongation_only_control_merged_input)



ear_elongation_only_ear_H2AZ_plots <- generate_final_plot_scaled(ear_elongation_only_ear_H2AZ_matrix_normalized, "H2AZ", H2A_colors, final_split)
ear_elongation_only_ear_H3K27me3_plots <- generate_final_plot_scaled(ear_elongation_only_ear_H3K27me3_matrix_normalized, "H3K27me3", H3K27me3_colors, final_split)
ear_elongation_only_ear_H3K36me3_plots <- generate_final_plot_scaled(ear_elongation_only_ear_H3K36me3_matrix_normalized, "H3K36me3", H3K36me3_colors, final_split)
ear_elongation_only_ear_H3K4me1_plots <- generate_final_plot_scaled(ear_elongation_only_ear_H3K4me1_matrix_normalized, "H3K4me1", H3K4me1_colors, final_split)
ear_elongation_only_ear_H3K4me3_plots <- generate_final_plot_scaled(ear_elongation_only_ear_H3K4me3_matrix_normalized, "H3K4me3", H3K4me3_colors, final_split)
ear_elongation_only_ear_H3K56ac_plots <- generate_final_plot_scaled(ear_elongation_only_ear_H3K56ac_matrix_normalized, "H3K56ac", H3K56ac_colors, final_split)



ear_elongation_only_ear_plot <- (ear_elongation_only_ear_H3K4me3_plots + 
                                   ear_elongation_only_ear_H3K56ac_plots + 
                                   ear_elongation_only_ear_H3K36me3_plots + 
                                   ear_elongation_only_ear_H3K4me1_plots + 
                                   ear_elongation_only_ear_H2AZ_plots + 
                                   ear_elongation_only_ear_H3K27me3_plots)


units_apart <- c(.7)
draw_units <- rep(units_apart, 6)

ear_elongation_only_ear_plot_draw = grid.grabExpr(draw(ear_elongation_only_ear_plot, ht_gap = unit(draw_units, "cm")))
ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/", filename="ear_elongation_only_ear_TSS_non_scaled.png", plot=ear_elongation_only_ear_plot_draw, width = 15, height = 18, units = "in")
system("open /Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/ear_elongation_only_ear_TSS_non_scaled.png")




