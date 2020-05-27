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
H2A_colors <-  "#66c2a5"
H3K27me3_colors <- "#234bb8"

generate_tss_heatmap_scaled_region <- function(TSS_file, mark_name) {
  
  #TSS_file <- TSS_file[, -c(151:200)] # delete columns 5 through 7
  
  TSS <- as.normalizedMatrix(TSS_file, 
                             k_upstream = 50, 
                             k_downstream = 50, 
                             k_target = 250,
                             extend = c(50, 50), 
                             signal_name = mark_name, 
                             target_name = c("TSS","TTS"),
                             keep = c(0,.95), smooth = FALSE)
  
  
  
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
  
  axis_name = c("-1000bp","TSS","TTS","1000bp")
  final_graph <-  EnrichedHeatmap(TSS_matrix,  column_title =  mark_name, show_heatmap_legend = TRUE, use_raster = TRUE, 
                                  row_split = row_split_array, cluster_rows = TRUE, 
                                  column_title_gp = gpar(fontsize = 25, fontface = "bold"), 
                                  col = col_fun, axis_name_rot = 90, row_title_rot = 0, axis_name_gp = gpar(fontsize = 18), axis_name = axis_name,pos_line_gp = gpar(lty = 3),
                                  top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4)))) 
  
  
  
  
  
  return(final_graph)
  
  
}
generate_complex_heatmaps_ATAC_scaled_region <- function(final_matrix_list, mark_name, row_split_array, color_hex) {
  
  #If row title is NOT missing - pass on it and don't add it 
  
  TSS_matrix <- final_matrix_list
  common_max <-  max(TSS_matrix)
  #col_fun <- circlize::colorRamp2(c(0, common_max), c("white", color_hex))
  col_fun = colorRamp2(quantile(TSS_matrix, c(0, 0.95)), c("white", color_hex))
  
  axis_name = c("-1000bp","TSS","TTS","1000bp")
  #This is the only difference between ATAC and ChIP-seq metaplots. Basically we can't scale the same way we did 
  #previously. This allows us to scale the top part of the metaplot appropriatly.
  
  final_graph <-  EnrichedHeatmap(TSS_matrix,  column_title =  mark_name, show_heatmap_legend = TRUE, use_raster = TRUE,
                                  row_split = row_split_array, cluster_rows = TRUE, 
                                  column_title_gp = gpar(fontsize = 25, fontface = "bold"), axis_name = axis_name, pos_line_gp = gpar(lty = 3),
                                  col = col_fun, axis_name_rot = 90, row_title_rot = 0, axis_name_gp = gpar(fontsize = 18))
  return(final_graph)
}
normalize_to_input <- function(mark_name_1, input_mark){
  normalized_to_input <- (mark_name_1 + 1)/(input_mark +1)
  return(normalized_to_input)
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


# leaf --------------------------------------------------------------------
leaf_initiation_only_leaf_H2AZ_all <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_initiation_only_leaf_H2AZ_all.gz")
leaf_initiation_only_leaf_H3K27me3_all <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_initiation_only_leaf_H3K27me3_all.gz")
leaf_initiation_only_leaf_H3K36me3_all <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_initiation_only_leaf_H3K36me3_all.gz")
leaf_initiation_only_leaf_H3K4me1_all <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_initiation_only_leaf_H3K4me1_all.gz")
leaf_initiation_only_leaf_H3K4me3_all <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_initiation_only_leaf_H3K4me3_all.gz")
leaf_initiation_only_leaf_H3K56ac_all <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_initiation_only_leaf_H3K56ac_all.gz")
leaf_initiation_only_leaf_input_all <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_initiation_only_leaf_input_all.gz")


leaf_initiation_only_leaf_H2AZ_matrix <- generate_tss_heatmap_scaled_region(leaf_initiation_only_leaf_H2AZ_all, "H2AZ")
leaf_initiation_only_leaf_H3K27me3_matrix <- generate_tss_heatmap_scaled_region(leaf_initiation_only_leaf_H3K27me3_all, "H3K27me")
leaf_initiation_only_leaf_H3K36me3_matrix <- generate_tss_heatmap_scaled_region(leaf_initiation_only_leaf_H3K36me3_all, "H3K36me")
leaf_initiation_only_leaf_H3K4me1_matrix <- generate_tss_heatmap_scaled_region(leaf_initiation_only_leaf_H3K4me1_all, "H3K4me1")
leaf_initiation_only_leaf_H3K4me3_matrix <- generate_tss_heatmap_scaled_region(leaf_initiation_only_leaf_H3K4me3_all, "H3K4me3")
leaf_initiation_only_leaf_H3K56ac_matrix <- generate_tss_heatmap_scaled_region(leaf_initiation_only_leaf_H3K56ac_all, "H3K56ac")
leaf_initiation_only_leaf_input_matrix <- generate_tss_heatmap_scaled_region(leaf_initiation_only_leaf_input_all, "input")







leaf_initiation_only_leaf_H2AZ_matrix_normalized <- normalize_to_input(leaf_initiation_only_leaf_H2AZ_matrix, leaf_initiation_only_leaf_input_matrix)
leaf_initiation_only_leaf_H3K27me3_matrix_normalized <- normalize_to_input(leaf_initiation_only_leaf_H3K27me3_matrix, leaf_initiation_only_leaf_input_matrix)
leaf_initiation_only_leaf_H3K36me3_matrix_normalized <- normalize_to_input(leaf_initiation_only_leaf_H3K36me3_matrix, leaf_initiation_only_leaf_input_matrix)
leaf_initiation_only_leaf_H3K4me1_matrix_normalized <- normalize_to_input(leaf_initiation_only_leaf_H3K4me1_matrix, leaf_initiation_only_leaf_input_matrix)
leaf_initiation_only_leaf_H3K4me3_matrix_normalized <- normalize_to_input(leaf_initiation_only_leaf_H3K4me3_matrix, leaf_initiation_only_leaf_input_matrix)
leaf_initiation_only_leaf_H3K56ac_matrix_normalized <- normalize_to_input(leaf_initiation_only_leaf_H3K56ac_matrix, leaf_initiation_only_leaf_input_matrix)


number_on <- nrow(leaf_initiation_only_leaf_H2AZ_matrix_normalized)
values <- c( "Leaf Init Only")
leaf_initiation_only_leaf_split <- rep(values, times = c(number_on))


leaf_initiation_only_leaf_H2AZ_plot <- generate_complex_heatmaps_chip_scaled_region(leaf_initiation_only_leaf_H2AZ_matrix_normalized, "H2NA", leaf_initiation_only_leaf_split, H2A_colors)
leaf_initiation_only_leaf_H3K27me3_plot <- generate_complex_heatmaps_chip_scaled_region(leaf_initiation_only_leaf_H3K27me3_matrix_normalized, "H3K27me3", leaf_initiation_only_leaf_split, H3K27me3_colors)
leaf_initiation_only_leaf_H3K36me3_plot <- generate_complex_heatmaps_chip_scaled_region(leaf_initiation_only_leaf_H3K36me3_matrix_normalized, "H3K36me3", leaf_initiation_only_leaf_split, H3K36me3_colors)
leaf_initiation_only_leaf_H3K4me1_plot <- generate_complex_heatmaps_chip_scaled_region(leaf_initiation_only_leaf_H3K4me1_matrix_normalized, "H3K4me1",leaf_initiation_only_leaf_split,  H3K4me1_colors)
leaf_initiation_only_leaf_H3K4me3_plot <- generate_complex_heatmaps_chip_scaled_region(leaf_initiation_only_leaf_H3K4me3_matrix_normalized, "H3K4me3",leaf_initiation_only_leaf_split,  H3K4me3_colors)
leaf_initiation_only_leaf_H3K56ac_plot <- generate_complex_heatmaps_chip_scaled_region(leaf_initiation_only_leaf_H3K56ac_matrix_normalized, "H3K56ac",leaf_initiation_only_leaf_split,  H3K56ac_colors)


leaf_initiation_only_leaf_plot <- (leaf_initiation_only_leaf_H3K36me3_plot + 
leaf_initiation_only_leaf_H3K4me1_plot + 
leaf_initiation_only_leaf_H3K4me3_plot + 
leaf_initiation_only_leaf_H3K56ac_plot + 
leaf_initiation_only_leaf_H2AZ_plot + 
leaf_initiation_only_leaf_H3K27me3_plot)


units_apart <- c(.7)
draw_units <- rep(units_apart, 6)

leaf_initiation_only_leaf_plot_draw = grid.grabExpr(draw(leaf_initiation_only_leaf_plot, ht_gap = unit(draw_units, "cm")))
ggsave(path = "~/Desktop/", filename="leaf_initiation_only_leaf_plot_draw.pdf", plot=leaf_initiation_only_leaf_plot_draw, width = 15, height = 18, units = "in")
system("open ~/Desktop/leaf_final_on_off.pdf")



# Root Init Only  ---------------------------------------------------------

root_initiation_only_root_H2AZ_all <- read_gunzipped_file("07.scaled_to_region/CHIP/root_initiation_only_root_H2AZ_all.gz")
root_initiation_only_root_H3K27me3_all <- read_gunzipped_file("07.scaled_to_region/CHIP/root_initiation_only_root_H3K27me3_all.gz")
root_initiation_only_root_H3K36me3_all <- read_gunzipped_file("07.scaled_to_region/CHIP/root_initiation_only_root_H3K36me3_all.gz")
root_initiation_only_root_H3K4me1_all <- read_gunzipped_file("07.scaled_to_region/CHIP/root_initiation_only_root_H3K4me1_all.gz")
root_initiation_only_root_H3K4me3_all <- read_gunzipped_file("07.scaled_to_region/CHIP/root_initiation_only_root_H3K4me3_all.gz")
root_initiation_only_root_H3K56ac_all <- read_gunzipped_file("07.scaled_to_region/CHIP/root_initiation_only_root_H3K56ac_all.gz")
root_initiation_only_root_input_all <- read_gunzipped_file("07.scaled_to_region/CHIP/root_initiation_only_root_input_all.gz")


root_initiation_only_root_H2AZ_matrix <- generate_tss_heatmap_scaled_region(root_initiation_only_root_H2AZ_all, "H2AZ")
root_initiation_only_root_H3K27me3_matrix <- generate_tss_heatmap_scaled_region(root_initiation_only_root_H3K27me3_all, "H3K27me")
root_initiation_only_root_H3K36me3_matrix <- generate_tss_heatmap_scaled_region(root_initiation_only_root_H3K36me3_all, "H3K36me")
root_initiation_only_root_H3K4me1_matrix <- generate_tss_heatmap_scaled_region(root_initiation_only_root_H3K4me1_all, "H3K4me1")
root_initiation_only_root_H3K4me3_matrix <- generate_tss_heatmap_scaled_region(root_initiation_only_root_H3K4me3_all, "H3K4me3")
root_initiation_only_root_H3K56ac_matrix <- generate_tss_heatmap_scaled_region(root_initiation_only_root_H3K56ac_all, "H3K56ac")
root_initiation_only_root_input_matrix <- generate_tss_heatmap_scaled_region(root_initiation_only_root_input_all, "input")


root_initiation_only_root_H2AZ_matrix_normalized <- normalize_to_input(root_initiation_only_root_H2AZ_matrix, root_initiation_only_root_input_matrix)
root_initiation_only_root_H3K27me3_matrix_normalized <- normalize_to_input(root_initiation_only_root_H3K27me3_matrix, root_initiation_only_root_input_matrix)
root_initiation_only_root_H3K36me3_matrix_normalized <- normalize_to_input(root_initiation_only_root_H3K36me3_matrix, root_initiation_only_root_input_matrix)
root_initiation_only_root_H3K4me1_matrix_normalized <- normalize_to_input(root_initiation_only_root_H3K4me1_matrix, root_initiation_only_root_input_matrix)
root_initiation_only_root_H3K4me3_matrix_normalized <- normalize_to_input(root_initiation_only_root_H3K4me3_matrix, root_initiation_only_root_input_matrix)
root_initiation_only_root_H3K56ac_matrix_normalized <- normalize_to_input(root_initiation_only_root_H3K56ac_matrix, root_initiation_only_root_input_matrix)


number_on <- nrow(root_initiation_only_root_H2AZ_matrix_normalized)
values <- c( "root Init Only")
root_initiation_only_root_split <- rep(values, times = c(number_on))


root_initiation_only_root_H2AZ_plot <- generate_complex_heatmaps_chip_scaled_region(root_initiation_only_root_H2AZ_matrix_normalized, "H2NA", root_initiation_only_root_split, H2A_colors)
root_initiation_only_root_H3K27me3_plot <- generate_complex_heatmaps_chip_scaled_region(root_initiation_only_root_H3K27me3_matrix_normalized, "H3K27me3", root_initiation_only_root_split, H3K27me3_colors)
root_initiation_only_root_H3K36me3_plot <- generate_complex_heatmaps_chip_scaled_region(root_initiation_only_root_H3K36me3_matrix_normalized, "H3K36me3", root_initiation_only_root_split, H3K36me3_colors)
root_initiation_only_root_H3K4me1_plot <- generate_complex_heatmaps_chip_scaled_region(root_initiation_only_root_H3K4me1_matrix_normalized, "H3K4me1",root_initiation_only_root_split,  H3K4me1_colors)
root_initiation_only_root_H3K4me3_plot <- generate_complex_heatmaps_chip_scaled_region(root_initiation_only_root_H3K4me3_matrix_normalized, "H3K4me3",root_initiation_only_root_split,  H3K4me3_colors)
root_initiation_only_root_H3K56ac_plot <- generate_complex_heatmaps_chip_scaled_region(root_initiation_only_root_H3K56ac_matrix_normalized, "H3K56ac",root_initiation_only_root_split,  H3K56ac_colors)


root_initiation_only_root_plot <- (root_initiation_only_root_H3K36me3_plot + 
                                         root_initiation_only_root_H3K4me1_plot + 
                                         root_initiation_only_root_H3K4me3_plot + 
                                         root_initiation_only_root_H3K56ac_plot + 
                                         root_initiation_only_root_H2AZ_plot + 
                                         root_initiation_only_root_H3K27me3_plot)


units_apart <- c(.7)
draw_units <- rep(units_apart, 6)

root_initiation_only_root_plot_draw = grid.grabExpr(draw(root_initiation_only_root_plot, ht_gap = unit(draw_units, "cm")))
ggsave(path = "~/Desktop/", filename="root_initiation_only_root_plot_draw.pdf", plot=root_initiation_only_root_plot_draw, width = 15, height = 18, units = "in")
system("open ~/Desktop/root_initiation_only_root_plot_draw.pdf")




# ear Init ONly---------------------------------------------------------------------

ear_initiation_only_ear_H2AZ_all <- read_gunzipped_file("07.scaled_to_region/CHIP/ear_initiation_only_ear_H2AZ_all.gz")
ear_initiation_only_ear_H3K27me3_all <- read_gunzipped_file("07.scaled_to_region/CHIP/ear_initiation_only_ear_H3K27me3_all.gz")
ear_initiation_only_ear_H3K36me3_all <- read_gunzipped_file("07.scaled_to_region/CHIP/ear_initiation_only_ear_H3K36me3_all.gz")
ear_initiation_only_ear_H3K4me1_all <- read_gunzipped_file("07.scaled_to_region/CHIP/ear_initiation_only_ear_H3K4me1_all.gz")
ear_initiation_only_ear_H3K4me3_all <- read_gunzipped_file("07.scaled_to_region/CHIP/ear_initiation_only_ear_H3K4me3_all.gz")
ear_initiation_only_ear_H3K56ac_all <- read_gunzipped_file("07.scaled_to_region/CHIP/ear_initiation_only_ear_H3K56ac_all.gz")
ear_initiation_only_ear_input_all <- read_gunzipped_file("07.scaled_to_region/CHIP/ear_initiation_only_ear_input_all.gz")


ear_initiation_only_ear_H2AZ_matrix <- generate_tss_heatmap_scaled_region(ear_initiation_only_ear_H2AZ_all, "H2AZ")
ear_initiation_only_ear_H3K27me3_matrix <- generate_tss_heatmap_scaled_region(ear_initiation_only_ear_H3K27me3_all, "H3K27me")
ear_initiation_only_ear_H3K36me3_matrix <- generate_tss_heatmap_scaled_region(ear_initiation_only_ear_H3K36me3_all, "H3K36me")
ear_initiation_only_ear_H3K4me1_matrix <- generate_tss_heatmap_scaled_region(ear_initiation_only_ear_H3K4me1_all, "H3K4me1")
ear_initiation_only_ear_H3K4me3_matrix <- generate_tss_heatmap_scaled_region(ear_initiation_only_ear_H3K4me3_all, "H3K4me3")
ear_initiation_only_ear_H3K56ac_matrix <- generate_tss_heatmap_scaled_region(ear_initiation_only_ear_H3K56ac_all, "H3K56ac")
ear_initiation_only_ear_input_matrix <- generate_tss_heatmap_scaled_region(ear_initiation_only_ear_input_all, "input")


ear_initiation_only_ear_H2AZ_matrix_normalized <- normalize_to_input(ear_initiation_only_ear_H2AZ_matrix, ear_initiation_only_ear_input_matrix)
ear_initiation_only_ear_H3K27me3_matrix_normalized <- normalize_to_input(ear_initiation_only_ear_H3K27me3_matrix, ear_initiation_only_ear_input_matrix)
ear_initiation_only_ear_H3K36me3_matrix_normalized <- normalize_to_input(ear_initiation_only_ear_H3K36me3_matrix, ear_initiation_only_ear_input_matrix)
ear_initiation_only_ear_H3K4me1_matrix_normalized <- normalize_to_input(ear_initiation_only_ear_H3K4me1_matrix, ear_initiation_only_ear_input_matrix)
ear_initiation_only_ear_H3K4me3_matrix_normalized <- normalize_to_input(ear_initiation_only_ear_H3K4me3_matrix, ear_initiation_only_ear_input_matrix)
ear_initiation_only_ear_H3K56ac_matrix_normalized <- normalize_to_input(ear_initiation_only_ear_H3K56ac_matrix, ear_initiation_only_ear_input_matrix)


number_on <- nrow(ear_initiation_only_ear_H2AZ_matrix_normalized)
values <- c( "ear Init Only")
ear_initiation_only_ear_split <- rep(values, times = c(number_on))


ear_initiation_only_ear_H2AZ_plot <- generate_complex_heatmaps_chip_scaled_region(ear_initiation_only_ear_H2AZ_matrix_normalized, "H2NA", ear_initiation_only_ear_split, H2A_colors)
ear_initiation_only_ear_H3K27me3_plot <- generate_complex_heatmaps_chip_scaled_region(ear_initiation_only_ear_H3K27me3_matrix_normalized, "H3K27me3", ear_initiation_only_ear_split, H3K27me3_colors)
ear_initiation_only_ear_H3K36me3_plot <- generate_complex_heatmaps_chip_scaled_region(ear_initiation_only_ear_H3K36me3_matrix_normalized, "H3K36me3", ear_initiation_only_ear_split, H3K36me3_colors)
ear_initiation_only_ear_H3K4me1_plot <- generate_complex_heatmaps_chip_scaled_region(ear_initiation_only_ear_H3K4me1_matrix_normalized, "H3K4me1",ear_initiation_only_ear_split,  H3K4me1_colors)
ear_initiation_only_ear_H3K4me3_plot <- generate_complex_heatmaps_chip_scaled_region(ear_initiation_only_ear_H3K4me3_matrix_normalized, "H3K4me3",ear_initiation_only_ear_split,  H3K4me3_colors)
ear_initiation_only_ear_H3K56ac_plot <- generate_complex_heatmaps_chip_scaled_region(ear_initiation_only_ear_H3K56ac_matrix_normalized, "H3K56ac",ear_initiation_only_ear_split,  H3K56ac_colors)


ear_initiation_only_ear_plot <- (ear_initiation_only_ear_H3K36me3_plot + 
                                         ear_initiation_only_ear_H3K4me1_plot + 
                                         ear_initiation_only_ear_H3K4me3_plot + 
                                         ear_initiation_only_ear_H3K56ac_plot + 
                                         ear_initiation_only_ear_H2AZ_plot + 
                                         ear_initiation_only_ear_H3K27me3_plot)


units_apart <- c(.7)
draw_units <- rep(units_apart, 6)

ear_initiation_only_ear_plot_draw = grid.grabExpr(draw(ear_initiation_only_ear_plot, ht_gap = unit(draw_units, "cm")))
ggsave(path = "~/Desktop/", filename="ear_initiation_only_ear_plot_draw.pdf", plot=ear_initiation_only_ear_plot_draw, width = 15, height = 18, units = "in")
system("open ~/Desktop/ear_initiation_only_ear_plot_draw.pdf")



# leaf --------------------------------------------------------------------
leaf_elongation_only_leaf_H2AZ_all <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_elongation_only_leaf_H2AZ_all.gz")
leaf_elongation_only_leaf_H3K27me3_all <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_elongation_only_leaf_H3K27me3_all.gz")
leaf_elongation_only_leaf_H3K36me3_all <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_elongation_only_leaf_H3K36me3_all.gz")
leaf_elongation_only_leaf_H3K4me1_all <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_elongation_only_leaf_H3K4me1_all.gz")
leaf_elongation_only_leaf_H3K4me3_all <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_elongation_only_leaf_H3K4me3_all.gz")
leaf_elongation_only_leaf_H3K56ac_all <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_elongation_only_leaf_H3K56ac_all.gz")
leaf_elongation_only_leaf_input_all <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_elongation_only_leaf_input_all.gz")


leaf_elongation_only_leaf_H2AZ_matrix <- generate_tss_heatmap_scaled_region(leaf_elongation_only_leaf_H2AZ_all, "H2AZ")
leaf_elongation_only_leaf_H3K27me3_matrix <- generate_tss_heatmap_scaled_region(leaf_elongation_only_leaf_H3K27me3_all, "H3K27me")
leaf_elongation_only_leaf_H3K36me3_matrix <- generate_tss_heatmap_scaled_region(leaf_elongation_only_leaf_H3K36me3_all, "H3K36me")
leaf_elongation_only_leaf_H3K4me1_matrix <- generate_tss_heatmap_scaled_region(leaf_elongation_only_leaf_H3K4me1_all, "H3K4me1")
leaf_elongation_only_leaf_H3K4me3_matrix <- generate_tss_heatmap_scaled_region(leaf_elongation_only_leaf_H3K4me3_all, "H3K4me3")
leaf_elongation_only_leaf_H3K56ac_matrix <- generate_tss_heatmap_scaled_region(leaf_elongation_only_leaf_H3K56ac_all, "H3K56ac")
leaf_elongation_only_leaf_input_matrix <- generate_tss_heatmap_scaled_region(leaf_elongation_only_leaf_input_all, "input")


leaf_elongation_only_leaf_H2AZ_matrix_normalized <- normalize_to_input(leaf_elongation_only_leaf_H2AZ_matrix, leaf_elongation_only_leaf_input_matrix)
leaf_elongation_only_leaf_H3K27me3_matrix_normalized <- normalize_to_input(leaf_elongation_only_leaf_H3K27me3_matrix, leaf_elongation_only_leaf_input_matrix)
leaf_elongation_only_leaf_H3K36me3_matrix_normalized <- normalize_to_input(leaf_elongation_only_leaf_H3K36me3_matrix, leaf_elongation_only_leaf_input_matrix)
leaf_elongation_only_leaf_H3K4me1_matrix_normalized <- normalize_to_input(leaf_elongation_only_leaf_H3K4me1_matrix, leaf_elongation_only_leaf_input_matrix)
leaf_elongation_only_leaf_H3K4me3_matrix_normalized <- normalize_to_input(leaf_elongation_only_leaf_H3K4me3_matrix, leaf_elongation_only_leaf_input_matrix)
leaf_elongation_only_leaf_H3K56ac_matrix_normalized <- normalize_to_input(leaf_elongation_only_leaf_H3K56ac_matrix, leaf_elongation_only_leaf_input_matrix)


number_on <- nrow(leaf_elongation_only_leaf_H2AZ_matrix_normalized)
values <- c( "Leaf Init Only")
leaf_elongation_only_leaf_split <- rep(values, times = c(number_on))


leaf_elongation_only_leaf_H2AZ_plot <- generate_complex_heatmaps_chip_scaled_region(leaf_elongation_only_leaf_H2AZ_matrix_normalized, "H2NA", leaf_elongation_only_leaf_split, H2A_colors)
leaf_elongation_only_leaf_H3K27me3_plot <- generate_complex_heatmaps_chip_scaled_region(leaf_elongation_only_leaf_H3K27me3_matrix_normalized, "H3K27me3", leaf_elongation_only_leaf_split, H3K27me3_colors)
leaf_elongation_only_leaf_H3K36me3_plot <- generate_complex_heatmaps_chip_scaled_region(leaf_elongation_only_leaf_H3K36me3_matrix_normalized, "H3K36me3", leaf_elongation_only_leaf_split, H3K36me3_colors)
leaf_elongation_only_leaf_H3K4me1_plot <- generate_complex_heatmaps_chip_scaled_region(leaf_elongation_only_leaf_H3K4me1_matrix_normalized, "H3K4me1",leaf_elongation_only_leaf_split,  H3K4me1_colors)
leaf_elongation_only_leaf_H3K4me3_plot <- generate_complex_heatmaps_chip_scaled_region(leaf_elongation_only_leaf_H3K4me3_matrix_normalized, "H3K4me3",leaf_elongation_only_leaf_split,  H3K4me3_colors)
leaf_elongation_only_leaf_H3K56ac_plot <- generate_complex_heatmaps_chip_scaled_region(leaf_elongation_only_leaf_H3K56ac_matrix_normalized, "H3K56ac",leaf_elongation_only_leaf_split,  H3K56ac_colors)


leaf_elongation_only_leaf_plot <- (leaf_elongation_only_leaf_H3K36me3_plot + 
                                         leaf_elongation_only_leaf_H3K4me1_plot + 
                                         leaf_elongation_only_leaf_H3K4me3_plot + 
                                         leaf_elongation_only_leaf_H3K56ac_plot + 
                                         leaf_elongation_only_leaf_H2AZ_plot + 
                                         leaf_elongation_only_leaf_H3K27me3_plot)


units_apart <- c(.7)
draw_units <- rep(units_apart, 6)

leaf_elongation_only_leaf_plot_draw = grid.grabExpr(draw(leaf_elongation_only_leaf_plot, ht_gap = unit(draw_units, "cm")))
ggsave(path = "~/Desktop/", filename="leaf_elongation_only_leaf_plot_draw.pdf", plot=leaf_elongation_only_leaf_plot_draw, width = 15, height = 18, units = "in")
system("open ~/Desktop/leaf_final_on_off.pdf")




# ear --------------------------------------------------------------------
ear_elongation_only_ear_H2AZ_all <- read_gunzipped_file("07.scaled_to_region/CHIP/ear_elongation_only_ear_H2AZ_all.gz")
ear_elongation_only_ear_H3K27me3_all <- read_gunzipped_file("07.scaled_to_region/CHIP/ear_elongation_only_ear_H3K27me3_all.gz")
ear_elongation_only_ear_H3K36me3_all <- read_gunzipped_file("07.scaled_to_region/CHIP/ear_elongation_only_ear_H3K36me3_all.gz")
ear_elongation_only_ear_H3K4me1_all <- read_gunzipped_file("07.scaled_to_region/CHIP/ear_elongation_only_ear_H3K4me1_all.gz")
ear_elongation_only_ear_H3K4me3_all <- read_gunzipped_file("07.scaled_to_region/CHIP/ear_elongation_only_ear_H3K4me3_all.gz")
ear_elongation_only_ear_H3K56ac_all <- read_gunzipped_file("07.scaled_to_region/CHIP/ear_elongation_only_ear_H3K56ac_all.gz")
ear_elongation_only_ear_input_all <- read_gunzipped_file("07.scaled_to_region/CHIP/ear_elongation_only_ear_input_all.gz")


ear_elongation_only_ear_H2AZ_matrix <- generate_tss_heatmap_scaled_region(ear_elongation_only_ear_H2AZ_all, "H2AZ")
ear_elongation_only_ear_H3K27me3_matrix <- generate_tss_heatmap_scaled_region(ear_elongation_only_ear_H3K27me3_all, "H3K27me")
ear_elongation_only_ear_H3K36me3_matrix <- generate_tss_heatmap_scaled_region(ear_elongation_only_ear_H3K36me3_all, "H3K36me")
ear_elongation_only_ear_H3K4me1_matrix <- generate_tss_heatmap_scaled_region(ear_elongation_only_ear_H3K4me1_all, "H3K4me1")
ear_elongation_only_ear_H3K4me3_matrix <- generate_tss_heatmap_scaled_region(ear_elongation_only_ear_H3K4me3_all, "H3K4me3")
ear_elongation_only_ear_H3K56ac_matrix <- generate_tss_heatmap_scaled_region(ear_elongation_only_ear_H3K56ac_all, "H3K56ac")
ear_elongation_only_ear_input_matrix <- generate_tss_heatmap_scaled_region(ear_elongation_only_ear_input_all, "input")


ear_elongation_only_ear_H2AZ_matrix_normalized <- normalize_to_input(ear_elongation_only_ear_H2AZ_matrix, ear_elongation_only_ear_input_matrix)
ear_elongation_only_ear_H3K27me3_matrix_normalized <- normalize_to_input(ear_elongation_only_ear_H3K27me3_matrix, ear_elongation_only_ear_input_matrix)
ear_elongation_only_ear_H3K36me3_matrix_normalized <- normalize_to_input(ear_elongation_only_ear_H3K36me3_matrix, ear_elongation_only_ear_input_matrix)
ear_elongation_only_ear_H3K4me1_matrix_normalized <- normalize_to_input(ear_elongation_only_ear_H3K4me1_matrix, ear_elongation_only_ear_input_matrix)
ear_elongation_only_ear_H3K4me3_matrix_normalized <- normalize_to_input(ear_elongation_only_ear_H3K4me3_matrix, ear_elongation_only_ear_input_matrix)
ear_elongation_only_ear_H3K56ac_matrix_normalized <- normalize_to_input(ear_elongation_only_ear_H3K56ac_matrix, ear_elongation_only_ear_input_matrix)


number_on <- nrow(ear_elongation_only_ear_H2AZ_matrix_normalized)
values <- c( "ear Init Only")
ear_elongation_only_ear_split <- rep(values, times = c(number_on))


ear_elongation_only_ear_H2AZ_plot <- generate_complex_heatmaps_chip_scaled_region(ear_elongation_only_ear_H2AZ_matrix_normalized, "H2NA", ear_elongation_only_ear_split, H2A_colors)
ear_elongation_only_ear_H3K27me3_plot <- generate_complex_heatmaps_chip_scaled_region(ear_elongation_only_ear_H3K27me3_matrix_normalized, "H3K27me3", ear_elongation_only_ear_split, H3K27me3_colors)
ear_elongation_only_ear_H3K36me3_plot <- generate_complex_heatmaps_chip_scaled_region(ear_elongation_only_ear_H3K36me3_matrix_normalized, "H3K36me3", ear_elongation_only_ear_split, H3K36me3_colors)
ear_elongation_only_ear_H3K4me1_plot <- generate_complex_heatmaps_chip_scaled_region(ear_elongation_only_ear_H3K4me1_matrix_normalized, "H3K4me1",ear_elongation_only_ear_split,  H3K4me1_colors)
ear_elongation_only_ear_H3K4me3_plot <- generate_complex_heatmaps_chip_scaled_region(ear_elongation_only_ear_H3K4me3_matrix_normalized, "H3K4me3",ear_elongation_only_ear_split,  H3K4me3_colors)
ear_elongation_only_ear_H3K56ac_plot <- generate_complex_heatmaps_chip_scaled_region(ear_elongation_only_ear_H3K56ac_matrix_normalized, "H3K56ac",ear_elongation_only_ear_split,  H3K56ac_colors)


ear_elongation_only_ear_plot <- (ear_elongation_only_ear_H3K36me3_plot + 
                                     ear_elongation_only_ear_H3K4me1_plot + 
                                     ear_elongation_only_ear_H3K4me3_plot + 
                                     ear_elongation_only_ear_H3K56ac_plot + 
                                     ear_elongation_only_ear_H2AZ_plot + 
                                     ear_elongation_only_ear_H3K27me3_plot)


units_apart <- c(.7)
draw_units <- rep(units_apart, 6)

ear_elongation_only_ear_plot_draw = grid.grabExpr(draw(ear_elongation_only_ear_plot, ht_gap = unit(draw_units, "cm")))
ggsave(path = "~/Desktop/", filename="ear_elongation_only_ear_plot_draw.pdf", plot=ear_elongation_only_ear_plot_draw, width = 15, height = 18, units = "in")
system("open ~/Desktop/ear_final_on_off.pdf")




root_elongation_only_root_H2AZ_all <- read_gunzipped_file("07.scaled_to_region/CHIP/root_elongation_only_root_H2AZ_all.gz")
root_elongation_only_root_H3K27me3_all <- read_gunzipped_file("07.scaled_to_region/CHIP/root_elongation_only_root_H3K27me3_all.gz")
root_elongation_only_root_H3K36me3_all <- read_gunzipped_file("07.scaled_to_region/CHIP/root_elongation_only_root_H3K36me3_all.gz")
root_elongation_only_root_H3K4me1_all <- read_gunzipped_file("07.scaled_to_region/CHIP/root_elongation_only_root_H3K4me1_all.gz")
root_elongation_only_root_H3K4me3_all <- read_gunzipped_file("07.scaled_to_region/CHIP/root_elongation_only_root_H3K4me3_all.gz")
root_elongation_only_root_H3K56ac_all <- read_gunzipped_file("07.scaled_to_region/CHIP/root_elongation_only_root_H3K56ac_all.gz")
root_elongation_only_root_input_all <- read_gunzipped_file("07.scaled_to_region/CHIP/root_elongation_only_root_input_all.gz")


root_elongation_only_root_H2AZ_matrix <- generate_tss_heatmap_scaled_region(root_elongation_only_root_H2AZ_all, "H2AZ")
root_elongation_only_root_H3K27me3_matrix <- generate_tss_heatmap_scaled_region(root_elongation_only_root_H3K27me3_all, "H3K27me")
root_elongation_only_root_H3K36me3_matrix <- generate_tss_heatmap_scaled_region(root_elongation_only_root_H3K36me3_all, "H3K36me")
root_elongation_only_root_H3K4me1_matrix <- generate_tss_heatmap_scaled_region(root_elongation_only_root_H3K4me1_all, "H3K4me1")
root_elongation_only_root_H3K4me3_matrix <- generate_tss_heatmap_scaled_region(root_elongation_only_root_H3K4me3_all, "H3K4me3")
root_elongation_only_root_H3K56ac_matrix <- generate_tss_heatmap_scaled_region(root_elongation_only_root_H3K56ac_all, "H3K56ac")
root_elongation_only_root_input_matrix <- generate_tss_heatmap_scaled_region(root_elongation_only_root_input_all, "input")


root_elongation_only_root_H2AZ_matrix_normalized <- normalize_to_input(root_elongation_only_root_H2AZ_matrix, root_elongation_only_root_input_matrix)
root_elongation_only_root_H3K27me3_matrix_normalized <- normalize_to_input(root_elongation_only_root_H3K27me3_matrix, root_elongation_only_root_input_matrix)
root_elongation_only_root_H3K36me3_matrix_normalized <- normalize_to_input(root_elongation_only_root_H3K36me3_matrix, root_elongation_only_root_input_matrix)
root_elongation_only_root_H3K4me1_matrix_normalized <- normalize_to_input(root_elongation_only_root_H3K4me1_matrix, root_elongation_only_root_input_matrix)
root_elongation_only_root_H3K4me3_matrix_normalized <- normalize_to_input(root_elongation_only_root_H3K4me3_matrix, root_elongation_only_root_input_matrix)
root_elongation_only_root_H3K56ac_matrix_normalized <- normalize_to_input(root_elongation_only_root_H3K56ac_matrix, root_elongation_only_root_input_matrix)


number_on <- nrow(root_elongation_only_root_H2AZ_matrix_normalized)
values <- c( "root Init Only")
root_elongation_only_root_split <- rep(values, times = c(number_on))


root_elongation_only_root_H2AZ_plot <- generate_complex_heatmaps_chip_scaled_region(root_elongation_only_root_H2AZ_matrix_normalized, "H2NA", root_elongation_only_root_split, H2A_colors)
root_elongation_only_root_H3K27me3_plot <- generate_complex_heatmaps_chip_scaled_region(root_elongation_only_root_H3K27me3_matrix_normalized, "H3K27me3", root_elongation_only_root_split, H3K27me3_colors)
root_elongation_only_root_H3K36me3_plot <- generate_complex_heatmaps_chip_scaled_region(root_elongation_only_root_H3K36me3_matrix_normalized, "H3K36me3", root_elongation_only_root_split, H3K36me3_colors)
root_elongation_only_root_H3K4me1_plot <- generate_complex_heatmaps_chip_scaled_region(root_elongation_only_root_H3K4me1_matrix_normalized, "H3K4me1",root_elongation_only_root_split,  H3K4me1_colors)
root_elongation_only_root_H3K4me3_plot <- generate_complex_heatmaps_chip_scaled_region(root_elongation_only_root_H3K4me3_matrix_normalized, "H3K4me3",root_elongation_only_root_split,  H3K4me3_colors)
root_elongation_only_root_H3K56ac_plot <- generate_complex_heatmaps_chip_scaled_region(root_elongation_only_root_H3K56ac_matrix_normalized, "H3K56ac",root_elongation_only_root_split,  H3K56ac_colors)


root_elongation_only_root_plot <- (root_elongation_only_root_H3K36me3_plot + 
                                         root_elongation_only_root_H3K4me1_plot + 
                                         root_elongation_only_root_H3K4me3_plot + 
                                         root_elongation_only_root_H3K56ac_plot + 
                                         root_elongation_only_root_H2AZ_plot + 
                                         root_elongation_only_root_H3K27me3_plot)


units_apart <- c(.7)
draw_units <- rep(units_apart, 6)

root_elongation_only_root_plot_draw = grid.grabExpr(draw(root_elongation_only_root_plot, ht_gap = unit(draw_units, "cm")))
ggsave(path = "~/Desktop/", filename="root_elongation_only_root_plot_draw.pdf", plot=root_elongation_only_root_plot_draw, width = 15, height = 18, units = "in")
system("open ~/Desktop/root_elongation_only_root_plot_draw.pdf")



