library(tidyverse)
library(EnrichedHeatmap)
library(cowplot)
library("grid")
library("ggplotify")
library("here")
library(circlize)

# load required packages
library(factoextra)
library(NbClust)


setwd("/Users/feilab/Projects/03.ncRNA_project/02.Analysis/lncRNA_copy_files/2020-06-10_intiation_intragenic_class")



H3K36me3_colors <- "#CB2026"
H3K4me1_colors <-  "#642165"
H3K4me3_colors <- "#B89BC9"
H3K56ac_colors <- "#99CA3C"
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
                             extend = c(100, 100), 
                             signal_name = mark_name, 
                             target_name = c("TSS"),
                             keep = c(0,.99), smooth = TRUE)
  
  
  
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
  
  axis_name = c("TSS")
  final_graph <-  EnrichedHeatmap(TSS_matrix,  column_title =  mark_name, show_heatmap_legend = FALSE, use_raster = TRUE, 
                                  row_split = row_split_array, cluster_rows = TRUE, 
                                  column_title_gp = gpar(fontsize = 25, fontface = "bold"), 
                                  col = col_fun, axis_name_rot = 90, row_title_rot = 90, axis_name_gp = gpar(fontsize = 18), axis_name = axis_name,pos_line_gp = gpar(lty = 3),
                                  top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 1:15)))) 
  
  
  
  
  
  return(final_graph)
  
  
}
normalize_to_input <- function(mark_name_1, input_mark){
  normalized_to_input <- (mark_name_1 + 1) - (input_mark +1)
  return(normalized_to_input)
}
combine_identical_marks_on_off <- function(genes_on, genes_off) {
  
  #genes_1 <- map_df(genes) 
  #lncRNA_1 <- map_df(lncRNA)
  #unannotated_1 <- map_df(unannotated)
  
  bind_all_rows <- bind_rows(
    as.data.frame(genes_on),
    as.data.frame(genes_off))
  
  final_return <- as.matrix(bind_all_rows)
  
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
generate_final_plot_scaled <- function(file_1, mark_name, color_hex, row_split_array, mark_type) {
  
  if (mark_type == "R2C2"){
    TSS_heatmap <- generate_tss_heatmap_r2c2_scaled_region(file_1, mark_name)
    R2C2_graphed_heatmap <- generate_complex_heatmaps_R2C2_scaled_region(TSS_heatmap, mark_name, row_split_array, color_hex)
    return(R2C2_graphed_heatmap)
    
  } else if (mark_type == "ATAC") {
    ATAC_heatmap <- generate_tss_heatmap_scaled_region(file_1, mark_name)
    ATAC_graphed_heatmap <- generate_complex_heatmaps_ATAC_scaled_region(ATAC_heatmap, mark_name, row_split_array, color_hex)
    return(ATAC_graphed_heatmap)
    
  } else if (mark_type == "CHIP") {
    Chip_heatmap <- generate_tss_heatmap_scaled_region(file_1, mark_name)
    Chip_plot_heatmap <- generate_complex_heatmaps_chip_scaled_region(Chip_heatmap, mark_name, row_split_array, color_hex)
    return(Chip_plot_heatmap)
  } 
}



dim(leaf_initiation_only_leaf_H3K36me3_all)



# leaf --------------------------------------------------------------------
leaf_initiation_only_leaf_H2AZ_all <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_initiation_only_no_gene_leaf_H2AZ_center.gz")
leaf_initiation_only_leaf_H3K27me3_all <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_initiation_only_no_gene_leaf_H3K27me3_center.gz")
leaf_initiation_only_leaf_H3K36me3_all <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_initiation_only_no_gene_leaf_H3K36me3_center.gz")
leaf_initiation_only_leaf_H3K4me1_all <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_initiation_only_no_gene_leaf_H3K4me1_center.gz")
leaf_initiation_only_leaf_H3K4me3_all <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_initiation_only_no_gene_leaf_H3K4me3_center.gz")
leaf_initiation_only_leaf_H3K56ac_all <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_initiation_only_no_gene_leaf_H3K56ac_center.gz")
leaf_initiation_only_leaf_input_all <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/leaf_initiation_only_no_gene_leaf_input_center.gz")



leaf_initiation_only_control_merged_H2AZ <-leaf_initiation_only_leaf_H2AZ_all
leaf_initiation_only_control_merged_H3K27me3 <-leaf_initiation_only_leaf_H3K27me3_all
leaf_initiation_only_control_merged_H3K36me3 <-leaf_initiation_only_leaf_H3K36me3_all
leaf_initiation_only_control_merged_H3K4me1 <-leaf_initiation_only_leaf_H3K4me1_all
leaf_initiation_only_control_merged_H3K4me3 <-leaf_initiation_only_leaf_H3K4me3_all
leaf_initiation_only_control_merged_H3K56ac <-leaf_initiation_only_leaf_H3K56ac_all
leaf_initiation_only_control_merged_input <-leaf_initiation_only_leaf_input_all



leaf_initiation_only_leaf_H2AZ_matrix_normalized <- normalize_to_input(leaf_initiation_only_control_merged_H2AZ, leaf_initiation_only_control_merged_input)
leaf_initiation_only_leaf_H3K27me3_matrix_normalized <- normalize_to_input(leaf_initiation_only_control_merged_H3K27me3, leaf_initiation_only_control_merged_input)
leaf_initiation_only_leaf_H3K36me3_matrix_normalized <- normalize_to_input(leaf_initiation_only_control_merged_H3K36me3, leaf_initiation_only_control_merged_input)
leaf_initiation_only_leaf_H3K4me1_matrix_normalized <- normalize_to_input(leaf_initiation_only_control_merged_H3K4me1, leaf_initiation_only_control_merged_input)
leaf_initiation_only_leaf_H3K4me3_matrix_normalized <- normalize_to_input(leaf_initiation_only_control_merged_H3K4me3, leaf_initiation_only_control_merged_input)
leaf_initiation_only_leaf_H3K56ac_matrix_normalized <- normalize_to_input(leaf_initiation_only_control_merged_H3K56ac, leaf_initiation_only_control_merged_input)



leaf_combined_all_matricies <- cbind(leaf_initiation_only_leaf_H3K36me3_matrix_normalized,leaf_initiation_only_leaf_H3K4me1_matrix_normalized,
            leaf_initiation_only_leaf_H3K4me3_matrix_normalized,leaf_initiation_only_leaf_H3K56ac_matrix_normalized)



wss <- function(k) {
  kmeans(leaf_combined_all_matricies, k, nstart = 30 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

png("/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/Initiaion_only_non_gene_class/leaf_squares_test.png")
plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")


dev.off()






x <- c(2:8)
for (val in x) {
  
  
  
  
  leaf_partition = paste0("cluster", kmeans(leaf_combined_all_matricies, centers = val)$cluster)
  leaf_partition_hm <- Heatmap(leaf_partition, col = structure(1:val, names = paste0("cluster", 1:val)), name = "partition",
                               show_row_names = FALSE, width = unit(3, "mm"))
  
  
  leaf_initiation_only_leaf_H2AZ_plots <- generate_final_plot_scaled(leaf_initiation_only_leaf_H2AZ_matrix_normalized, "H2AZ", H2A_colors, leaf_partition, "CHIP")
  leaf_initiation_only_leaf_H3K27me3_plots <- generate_final_plot_scaled(leaf_initiation_only_leaf_H3K27me3_matrix_normalized, "H3K27me3", H3K27me3_colors, leaf_partition, "CHIP")
  leaf_initiation_only_leaf_H3K36me3_plots <- generate_final_plot_scaled(leaf_initiation_only_leaf_H3K36me3_matrix_normalized, "H3K36me3", H3K36me3_colors, leaf_partition, "CHIP")
  leaf_initiation_only_leaf_H3K4me1_plots <- generate_final_plot_scaled(leaf_initiation_only_leaf_H3K4me1_matrix_normalized, "H3K4me1", H3K4me1_colors, leaf_partition, "CHIP")
  leaf_initiation_only_leaf_H3K4me3_plots <- generate_final_plot_scaled(leaf_initiation_only_leaf_H3K4me3_matrix_normalized, "H3K4me3", H3K4me3_colors, leaf_partition, "CHIP")
  leaf_initiation_only_leaf_H3K56ac_plots <- generate_final_plot_scaled(leaf_initiation_only_leaf_H3K56ac_matrix_normalized, "H3K56ac", H3K56ac_colors, leaf_partition, "CHIP")
  
  
  
  leaf_initiation_only_leaf_plot <- (leaf_partition_hm + leaf_initiation_only_leaf_H3K4me3_plots + 
                                       leaf_initiation_only_leaf_H3K56ac_plots + 
                                       leaf_initiation_only_leaf_H3K36me3_plots + 
                                       leaf_initiation_only_leaf_H3K4me1_plots + 
                                       leaf_initiation_only_leaf_H2AZ_plots + 
                                       leaf_initiation_only_leaf_H3K27me3_plots)
  
  
  units_apart <- c(.7)
  draw_units <- rep(units_apart, 6)
  
  k_means_file_name_base = str_c("leaf_initiation_only_leaf_plot_kmeans_cluster", as.character(val), sep = '_')
  k_means_file_name = str_c(k_means_file_name_base, ".pdf", sep="_")
  leaf_initiation_only_leaf_plot_draw = grid.grabExpr(draw(leaf_initiation_only_leaf_plot, split = leaf_partition, ht_gap = unit(draw_units, "cm")))
  ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/Initiaion_only_non_gene_class", filename=k_means_file_name, plot=leaf_initiation_only_leaf_plot_draw, width = 15, height = 18, units = "in")
  
}



#leaf_partition = paste0("cluster", kmeans(leaf_combined_all_matricies, centers = 4)$cluster)
#leaf_partition_hm <- Heatmap(leaf_partition, col = structure(1:4, names = paste0("cluster", 1:4)), name = "partition",
#                        show_row_names = FALSE, width = unit(3, "mm"))
#
#
#
#
#leaf_initiation_only_leaf_H2AZ_plots <- generate_final_plot_scaled(leaf_initiation_only_leaf_H2AZ_matrix_normalized, "H2AZ", H2A_colors, leaf_partition, "CHIP")
#leaf_initiation_only_leaf_H3K27me3_plots <- generate_final_plot_scaled(leaf_initiation_only_leaf_H3K27me3_matrix_normalized, "H3K27me3", H3K27me3_colors, leaf_partition, "CHIP")
#leaf_initiation_only_leaf_H3K36me3_plots <- generate_final_plot_scaled(leaf_initiation_only_leaf_H3K36me3_matrix_normalized, "H3K36me3", H3K36me3_colors, leaf_partition, "CHIP")
#leaf_initiation_only_leaf_H3K4me1_plots <- generate_final_plot_scaled(leaf_initiation_only_leaf_H3K4me1_matrix_normalized, "H3K4me1", H3K4me1_colors, leaf_partition, "CHIP")
#leaf_initiation_only_leaf_H3K4me3_plots <- generate_final_plot_scaled(leaf_initiation_only_leaf_H3K4me3_matrix_normalized, "H3K4me3", H3K4me3_colors, leaf_partition, "CHIP")
#leaf_initiation_only_leaf_H3K56ac_plots <- generate_final_plot_scaled(leaf_initiation_only_leaf_H3K56ac_matrix_normalized, "H3K56ac", H3K56ac_colors, leaf_partition, "CHIP")
#
#
#
#leaf_initiation_only_leaf_plot <- (leaf_partition_hm + leaf_initiation_only_leaf_H3K4me3_plots + 
#                                     leaf_initiation_only_leaf_H3K56ac_plots + 
#                                     leaf_initiation_only_leaf_H3K36me3_plots + 
#                                     leaf_initiation_only_leaf_H3K4me1_plots + 
#                                     leaf_initiation_only_leaf_H2AZ_plots + 
#                                     leaf_initiation_only_leaf_H3K27me3_plots)
#
#
#units_apart <- c(.7)
#draw_units <- rep(units_apart, 6)
#
#leaf_initiation_only_leaf_plot_draw = grid.grabExpr(draw(leaf_initiation_only_leaf_plot, split = leaf_partition, ht_gap = unit(draw_units, "cm")))
#ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/04.initiation_overlapping_elongation", filename="leaf_initiation_only_leaf_plot_draw_2.pdf", plot=leaf_initiation_only_leaf_plot_draw, width = 15, height = 18, units = "in")
#system("open /Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/04.initiation_overlapping_elongation/leaf_initiation_only_leaf_plot_draw_2.pdf")
#






# root --------------------------------------------------------------------
root_initiation_only_root_H2AZ_all <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_initiation_only_no_gene_root_H2AZ_center.gz")
root_initiation_only_root_H3K27me3_all <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_initiation_only_no_gene_root_H3K27me3_center.gz")
root_initiation_only_root_H3K36me3_all <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_initiation_only_no_gene_root_H3K36me3_center.gz")
root_initiation_only_root_H3K4me1_all <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_initiation_only_no_gene_root_H3K4me1_center.gz")
root_initiation_only_root_H3K4me3_all <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_initiation_only_no_gene_root_H3K4me3_center.gz")
root_initiation_only_root_H3K56ac_all <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_initiation_only_no_gene_root_H3K56ac_center.gz")
root_initiation_only_root_input_all <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/root_initiation_only_no_gene_root_input_center.gz")



root_initiation_only_control_merged_H2AZ <-root_initiation_only_root_H2AZ_all
root_initiation_only_control_merged_H3K27me3 <-root_initiation_only_root_H3K27me3_all
root_initiation_only_control_merged_H3K36me3 <-root_initiation_only_root_H3K36me3_all
root_initiation_only_control_merged_H3K4me1 <-root_initiation_only_root_H3K4me1_all
root_initiation_only_control_merged_H3K4me3 <-root_initiation_only_root_H3K4me3_all
root_initiation_only_control_merged_H3K56ac <-root_initiation_only_root_H3K56ac_all
root_initiation_only_control_merged_input <-root_initiation_only_root_input_all



root_initiation_only_root_H2AZ_matrix_normalized <- normalize_to_input(root_initiation_only_control_merged_H2AZ, root_initiation_only_control_merged_input)
root_initiation_only_root_H3K27me3_matrix_normalized <- normalize_to_input(root_initiation_only_control_merged_H3K27me3, root_initiation_only_control_merged_input)
root_initiation_only_root_H3K36me3_matrix_normalized <- normalize_to_input(root_initiation_only_control_merged_H3K36me3, root_initiation_only_control_merged_input)
root_initiation_only_root_H3K4me1_matrix_normalized <- normalize_to_input(root_initiation_only_control_merged_H3K4me1, root_initiation_only_control_merged_input)
root_initiation_only_root_H3K4me3_matrix_normalized <- normalize_to_input(root_initiation_only_control_merged_H3K4me3, root_initiation_only_control_merged_input)
root_initiation_only_root_H3K56ac_matrix_normalized <- normalize_to_input(root_initiation_only_control_merged_H3K56ac, root_initiation_only_control_merged_input)



root_combined_all_matricies <- cbind(root_initiation_only_root_H3K36me3_matrix_normalized,root_initiation_only_root_H3K4me1_matrix_normalized,
                                     root_initiation_only_root_H3K4me3_matrix_normalized,root_initiation_only_root_H3K56ac_matrix_normalized)



wss <- function(k) {
  kmeans(root_combined_all_matricies, k, nstart = 30 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

root_test_squares <- plot(k.values, wss_values,
                          type="b", pch = 19, frame = FALSE, 
                          xlab="Number of clusters K",
                          ylab="Total within-clusters sum of squares")


ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/Initiaion_only_non_gene_class", filename="root_squares_test.png", plot=root_test_squares)

x <- c(2:8)
for (val in x) {
  
  
  
  
  root_partition = paste0("cluster", kmeans(root_combined_all_matricies, centers = val)$cluster)
  root_partition_hm <- Heatmap(root_partition, col = structure(1:val, names = paste0("cluster", 1:val)), name = "partition",
                               show_row_names = FALSE, width = unit(3, "mm"))
  
  
  root_initiation_only_root_H2AZ_plots <- generate_final_plot_scaled(root_initiation_only_root_H2AZ_matrix_normalized, "H2AZ", H2A_colors, root_partition, "CHIP")
  root_initiation_only_root_H3K27me3_plots <- generate_final_plot_scaled(root_initiation_only_root_H3K27me3_matrix_normalized, "H3K27me3", H3K27me3_colors, root_partition, "CHIP")
  root_initiation_only_root_H3K36me3_plots <- generate_final_plot_scaled(root_initiation_only_root_H3K36me3_matrix_normalized, "H3K36me3", H3K36me3_colors, root_partition, "CHIP")
  root_initiation_only_root_H3K4me1_plots <- generate_final_plot_scaled(root_initiation_only_root_H3K4me1_matrix_normalized, "H3K4me1", H3K4me1_colors, root_partition, "CHIP")
  root_initiation_only_root_H3K4me3_plots <- generate_final_plot_scaled(root_initiation_only_root_H3K4me3_matrix_normalized, "H3K4me3", H3K4me3_colors, root_partition, "CHIP")
  root_initiation_only_root_H3K56ac_plots <- generate_final_plot_scaled(root_initiation_only_root_H3K56ac_matrix_normalized, "H3K56ac", H3K56ac_colors, root_partition, "CHIP")
  
  
  
  root_initiation_only_root_plot <- (root_partition_hm + root_initiation_only_root_H3K4me3_plots + 
                                       root_initiation_only_root_H3K56ac_plots + 
                                       root_initiation_only_root_H3K36me3_plots + 
                                       root_initiation_only_root_H3K4me1_plots + 
                                       root_initiation_only_root_H2AZ_plots + 
                                       root_initiation_only_root_H3K27me3_plots)
  
  
  units_apart <- c(.7)
  draw_units <- rep(units_apart, 6)
  
  k_means_file_name_base = str_c("root_initiation_only_root_plot_kmeans_cluster", as.character(val), sep = '_')
  k_means_file_name = str_c(k_means_file_name_base, ".pdf", sep="_")
  root_initiation_only_root_plot_draw = grid.grabExpr(draw(root_initiation_only_root_plot, split = root_partition, ht_gap = unit(draw_units, "cm")))
  ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/Initiaion_only_non_gene_class", filename=k_means_file_name, plot=root_initiation_only_root_plot_draw, width = 15, height = 18, units = "in")
  
}



# ear --------------------------------------------------------------------
ear_initiation_only_ear_H2AZ_all <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_initiation_only_no_gene_ear_H2AZ_center.gz")
ear_initiation_only_ear_H3K27me3_all <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_initiation_only_no_gene_ear_H3K27me3_center.gz")
ear_initiation_only_ear_H3K36me3_all <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_initiation_only_no_gene_ear_H3K36me3_center.gz")
ear_initiation_only_ear_H3K4me1_all <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_initiation_only_no_gene_ear_H3K4me1_center.gz")
ear_initiation_only_ear_H3K4me3_all <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_initiation_only_no_gene_ear_H3K4me3_center.gz")
ear_initiation_only_ear_H3K56ac_all <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_initiation_only_no_gene_ear_H3K56ac_center.gz")
ear_initiation_only_ear_input_all <- read_gunzipped_file("05.Generate_matrix_TSS_TES_overlapping_H3K36me3/ear_initiation_only_no_gene_ear_input_center.gz")



ear_initiation_only_control_merged_H2AZ <-ear_initiation_only_ear_H2AZ_all
ear_initiation_only_control_merged_H3K27me3 <-ear_initiation_only_ear_H3K27me3_all
ear_initiation_only_control_merged_H3K36me3 <-ear_initiation_only_ear_H3K36me3_all
ear_initiation_only_control_merged_H3K4me1 <-ear_initiation_only_ear_H3K4me1_all
ear_initiation_only_control_merged_H3K4me3 <-ear_initiation_only_ear_H3K4me3_all
ear_initiation_only_control_merged_H3K56ac <-ear_initiation_only_ear_H3K56ac_all
ear_initiation_only_control_merged_input <-ear_initiation_only_ear_input_all



ear_initiation_only_ear_H2AZ_matrix_normalized <- normalize_to_input(ear_initiation_only_control_merged_H2AZ, ear_initiation_only_control_merged_input)
ear_initiation_only_ear_H3K27me3_matrix_normalized <- normalize_to_input(ear_initiation_only_control_merged_H3K27me3, ear_initiation_only_control_merged_input)
ear_initiation_only_ear_H3K36me3_matrix_normalized <- normalize_to_input(ear_initiation_only_control_merged_H3K36me3, ear_initiation_only_control_merged_input)
ear_initiation_only_ear_H3K4me1_matrix_normalized <- normalize_to_input(ear_initiation_only_control_merged_H3K4me1, ear_initiation_only_control_merged_input)
ear_initiation_only_ear_H3K4me3_matrix_normalized <- normalize_to_input(ear_initiation_only_control_merged_H3K4me3, ear_initiation_only_control_merged_input)
ear_initiation_only_ear_H3K56ac_matrix_normalized <- normalize_to_input(ear_initiation_only_control_merged_H3K56ac, ear_initiation_only_control_merged_input)



ear_combined_all_matricies <- cbind(ear_initiation_only_ear_H3K36me3_matrix_normalized,ear_initiation_only_ear_H3K4me1_matrix_normalized,
                                     ear_initiation_only_ear_H3K4me3_matrix_normalized,ear_initiation_only_ear_H3K56ac_matrix_normalized)



wss <- function(k) {
  kmeans(ear_combined_all_matricies, k, nstart = 30 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

png("/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/Initiaion_only_non_gene_class/ear_squares_test.png")
plot(k.values, wss_values,
                          type="b", pch = 19, frame = FALSE, 
                          xlab="Number of clusters K",
                          ylab="Total within-clusters sum of squares")


dev.off()


x <- c(2:8)
for (val in x) {
  
  
  
  
  ear_partition = paste0("cluster", kmeans(ear_combined_all_matricies, centers = val)$cluster)
  ear_partition_hm <- Heatmap(ear_partition, col = structure(1:val, names = paste0("cluster", 1:val)), name = "partition",
                               show_row_names = FALSE, width = unit(3, "mm"))
  
  
  ear_initiation_only_ear_H2AZ_plots <- generate_final_plot_scaled(ear_initiation_only_ear_H2AZ_matrix_normalized, "H2AZ", H2A_colors, ear_partition, "CHIP")
  ear_initiation_only_ear_H3K27me3_plots <- generate_final_plot_scaled(ear_initiation_only_ear_H3K27me3_matrix_normalized, "H3K27me3", H3K27me3_colors, ear_partition, "CHIP")
  ear_initiation_only_ear_H3K36me3_plots <- generate_final_plot_scaled(ear_initiation_only_ear_H3K36me3_matrix_normalized, "H3K36me3", H3K36me3_colors, ear_partition, "CHIP")
  ear_initiation_only_ear_H3K4me1_plots <- generate_final_plot_scaled(ear_initiation_only_ear_H3K4me1_matrix_normalized, "H3K4me1", H3K4me1_colors, ear_partition, "CHIP")
  ear_initiation_only_ear_H3K4me3_plots <- generate_final_plot_scaled(ear_initiation_only_ear_H3K4me3_matrix_normalized, "H3K4me3", H3K4me3_colors, ear_partition, "CHIP")
  ear_initiation_only_ear_H3K56ac_plots <- generate_final_plot_scaled(ear_initiation_only_ear_H3K56ac_matrix_normalized, "H3K56ac", H3K56ac_colors, ear_partition, "CHIP")
  
  
  
  ear_initiation_only_ear_plot <- (ear_partition_hm + ear_initiation_only_ear_H3K4me3_plots + 
                                       ear_initiation_only_ear_H3K56ac_plots + 
                                       ear_initiation_only_ear_H3K36me3_plots + 
                                       ear_initiation_only_ear_H3K4me1_plots + 
                                       ear_initiation_only_ear_H2AZ_plots + 
                                       ear_initiation_only_ear_H3K27me3_plots)
  
  
  units_apart <- c(.7)
  draw_units <- rep(units_apart, 6)
  
  k_means_file_name_base = str_c("ear_initiation_only_ear_plot_kmeans_cluster", as.character(val), sep = '_')
  k_means_file_name = str_c(k_means_file_name_base, ".pdf", sep="_")
  ear_initiation_only_ear_plot_draw = grid.grabExpr(draw(ear_initiation_only_ear_plot, split = ear_partition, ht_gap = unit(draw_units, "cm")))
  ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/Initiaion_only_non_gene_class", filename=k_means_file_name, plot=ear_initiation_only_ear_plot_draw, width = 15, height = 18, units = "in")
  
}
