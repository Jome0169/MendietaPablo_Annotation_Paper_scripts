library(tidyverse)
library(EnrichedHeatmap)
library(cowplot)
library("grid")
library("ggplotify")
library("here")
library(circlize)

setwd("/Users/feilab/Projects/03.ncRNA_project/02.Analysis/lncRNA_copy_files/2020-01-06_metaplots_lncRNA_RNAOnly_vs_Myset")

read_gunzipped_file <- function(zipped_file){
  
  bw_zipped_file <- read_delim(zipped_file, delim='\t', col_names = FALSE, skip = 1)
  
  final_return <- bw_zipped_file %>% 
    unite(combined_name, X1,X2,X3,X4,X5,X6, sep='_') %>% 
    column_to_rownames(var="combined_name")
  
  return_matrix <- as.matrix(final_return)
  
  return(return_matrix)
  
  
}
generate_tss_heatmap_2kb_up_down <- function(TSS_file, mark_name) {
  
  #TSS_file <- TSS_file[, -c(151:200)] # delete columns 5 through 7
  
  TSS <- as.normalizedMatrix(TSS_file, 
                             k_upstream = 400, 
                             k_downstream = 400, 
                             k_target = 0,
                             extend = c(400, 400), 
                             signal_name = mark_name, 
                             target_name = "TSS",
                             keep = c(0,.99), smooth = FALSE)
  
  
  
  return(TSS)
}
generate_tss_heatmap_r2c2_2kb_up_down <- function(TSS_file, mark_name) {
  
  #TSS_file <- TSS_file[, -c(151:200)] # delete columns 5 through 7
  
  TSS <- as.normalizedMatrix(TSS_file, 
                             k_upstream = 400, 
                             k_downstream = 400, 
                             k_target = 0,
                             extend = c(400, 400), 
                             background = 0, 
                             signal_name = mark_name, 
                             target_name = "TSS",
                             keep = c(0,.99), smooth = FALSE)
  
  
  
  return(TSS)
}
#Function to  generate the Heatmap. 
generate_complex_heatmaps_chip_2kb_up_down <- function(TSS_matrix, mark_name, color_hex) {
  
  #If row title is NOT missing - pass on it and don't add it 
  
  TSS_matrix <- TSS_matrix
  
  
  #Scale colors across matricies the same
  #common_min <-  min(TSS_matrix)
  common_max <-  max(TSS_matrix)
  col_fun <- circlize::colorRamp2(c(0, common_max), c("white", color_hex))
  
  #Get the same max value for each 
  #y_max_val <- round(max(colMeans(TES_matrix, na.rm = TRUE), colMeans(TSS_matrix, na.rm = TRUE))) + 2 
  
  axis_name = c("-2000bp","TSS","2000bp")
  final_graph <-  EnrichedHeatmap(TSS_matrix, cluster_rows = FALSE, column_title =  mark_name, show_heatmap_legend = FALSE, use_raster = TRUE, column_title_gp = gpar(fontsize = 30, fontface = "bold"), 
                                  col = col_fun, axis_name_rot = 90, row_title_rot = 0, axis_name_gp = gpar(fontsize = 18), axis_name = axis_name,pos_line_gp = gpar(lty = 3),
                                  top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4)))) 
  
  
  
  
  
  return(final_graph)
  
  
}
generate_complex_heatmaps_ATAC_2kb_up_down <- function(final_matrix_list, mark_name, color_hex) {
  
  #If row title is NOT missing - pass on it and don't add it 
  
  TSS_matrix <- final_matrix_list
  common_max <-  max(TSS_matrix)
  #col_fun <- circlize::colorRamp2(c(0, common_max), c("white", color_hex))
  col_fun = colorRamp2(quantile(TSS_matrix, c(0, 0.95)), c("white", color_hex))
  
  
  #This is the only difference between ATAC and ChIP-seq metaplots. Basically we can't scale the same way we did 
  #previously. This allows us to scale the top part of the metaplot appropriatly.
  axis_name = c("-2000bp","TSS","2000bp")
  final_graph <-  EnrichedHeatmap(TSS_matrix, cluster_rows = FALSE, column_title =  mark_name, show_heatmap_legend = FALSE, use_raster = TRUE,
                                  column_title_gp = gpar(fontsize = 30, fontface = "bold"), axis_name = axis_name, pos_line_gp = gpar(lty = 3),
                                  col = col_fun, axis_name_rot = 90, row_title_rot = 0, axis_name_gp = gpar(fontsize = 18))
  return(final_graph)
}
generate_complex_heatmaps_R2C2_2kb_up_down <- function(TSS_matrix, mark_name, color_hex) {
  
  #If row title is NOT missing - pass on it and don't add it 
  
  TSS_matrix <- TSS_matrix
  axis_name = c("-2000bp","TSS","2000bp")
  #common_max <-  max(TSS_matrix)
  #col_fun <- circlize::colorRamp2(c(0, common_max), c("white", color_hex))
  col_fun = colorRamp2(quantile(TSS_matrix, c(0, 0.98)), c("white", color_hex))
  final_graph <- EnrichedHeatmap(TSS_matrix, cluster_rows = FALSE, column_title =  mark_name, show_heatmap_legend = FALSE, use_raster = TRUE, 
                                 column_title_gp = gpar(fontsize = 30, fontface = "bold"), axis_name = axis_name,pos_line_gp = gpar(lty = 3),
                                 col = col_fun, axis_name_rot = 90, row_title_rot = 0, axis_name_gp = gpar(fontsize = 18))
  
  
  return(final_graph)
}
generate_final_plot_2kb_up_down <- function(file_1, mark_name, color_hex, mark_type) {
  
  if (mark_type == "R2C2"){
    TSS_heatmap <- generate_tss_heatmap_r2c2_2kb_up_down(file_1, mark_name)
    R2C2_graphed_heatmap <- generate_complex_heatmaps_R2C2_2kb_up_down(TSS_heatmap, mark_name, color_hex)
    return(R2C2_graphed_heatmap)
    
  } else if (mark_type == "ATAC") {
    ATAC_heatmap <- generate_tss_heatmap_2kb_up_down(file_1, mark_name)
    ATAC_graphed_heatmap <- generate_complex_heatmaps_ATAC_2kb_up_down(ATAC_heatmap, mark_name, color_hex)
    return(ATAC_graphed_heatmap)
    
  } else if (mark_type == "CHIP") {
    Chip_heatmap <- generate_tss_heatmap_2kb_up_down(file_1, mark_name)
    Chip_plot_heatmap <- generate_complex_heatmaps_chip_2kb_up_down(Chip_heatmap, mark_name, color_hex)
    return(Chip_plot_heatmap)
  } 
}



###############################3
#### WORK ON ORIGINAL ANNOTATION
###############################3
leaf_H3K56ac_R2C2_coord_og_2kb_up_down <- read_gunzipped_file("/Users/feilab/Projects/04.R2C2/01.analysis/R2C2_metaplots/2020-02-28_R2C2_metaplots/05.metaplot_matrix/CHIP/Zm-B73-REFERENCE.all.1kb_away_from_other_tis_leaf_mod_H3K56ac_TSS.gz")
leaf_H3K4me3_R2C2_coord_og_2kb_up_down <- read_gunzipped_file("/Users/feilab/Projects/04.R2C2/01.analysis/R2C2_metaplots/2020-02-28_R2C2_metaplots/05.metaplot_matrix/CHIP/Zm-B73-REFERENCE.all.1kb_away_from_other_tis_leaf_mod_H3K4me3_TSS.gz")
leaf_R2C2_R2C2_coord_og_2kb_up_down <- read_gunzipped_file("/Users/feilab/Projects/04.R2C2/01.analysis/R2C2_metaplots/2020-02-28_R2C2_metaplots/05.metaplot_matrix/R2C2/Zm-B73-REFERENCE.all.1kb_away_from_other_tis_leaf_mod_R2C2_TSS.gz")
leaf_ATAC_R2C2_coord_og_2kb_up_down <- read_gunzipped_file("/Users/feilab/Projects/04.R2C2/01.analysis/R2C2_metaplots/2020-02-28_R2C2_metaplots/05.metaplot_matrix/ATAC/Zm-B73-REFERENCE.all.1kb_away_from_other_tis_leaf_mod_ATAC_nonscaled_TSS.gz")


#replace all NAs with a zero ---------------------------------------------
leaf_H3K56ac_R2C2_coord_og_2kb_up_down[is.nan(leaf_H3K56ac_R2C2_coord_og_2kb_up_down)] = 0
leaf_H3K4me3_R2C2_coord_og_2kb_up_down[is.nan(leaf_H3K4me3_R2C2_coord_og_2kb_up_down)] = 0
leaf_R2C2_R2C2_coord_og_2kb_up_down[is.nan(leaf_R2C2_R2C2_coord_og_2kb_up_down)] = 0
leaf_ATAC_R2C2_coord_og_2kb_up_down[is.nan(leaf_ATAC_R2C2_coord_og_2kb_up_down)] = 0


leaf_H3K56ac_R2C2_plot_2kb_up_down_og <- generate_final_plot_2kb_up_down(leaf_H3K56ac_R2C2_coord_og_2kb_up_down, "H3K56ac", "#B89BC9", "CHIP")
leaf_H3K4me3_R2C2_plot_2kb_up_down_og <- generate_final_plot_2kb_up_down(leaf_H3K4me3_R2C2_coord_og_2kb_up_down, "H3K4me3", "#99CA3C", "CHIP")
leaf_R2C2_R2C2_plot_2kb_up_down_og <- generate_final_plot_2kb_up_down(leaf_R2C2_R2C2_coord_og_2kb_up_down, "R2C2", "#FFCD05", "R2C2")
leaf_ATAC_R2C2_plot_2kb_up_down_og <- generate_final_plot_2kb_up_down(leaf_ATAC_R2C2_coord_og_2kb_up_down, "ATAC", "#66c2a5", "ATAC")


final_graph_graph_original_coord <- (leaf_R2C2_R2C2_plot_2kb_up_down_og + leaf_ATAC_R2C2_plot_2kb_up_down_og + 
                                       leaf_H3K4me3_R2C2_plot_2kb_up_down_og + leaf_H3K56ac_R2C2_plot_2kb_up_down_og)


units_apart <- c(.7)
draw_units <- rep(units_apart, 4)
R2C2_final_original_coord_2kb_up_down = grid.grabExpr(draw(final_graph_graph_original_coord, ht_gap = unit(draw_units, "cm")))
ggsave(path = "~/Desktop/R2C2_metaplot", filename="original_2kb_up_down.pdf", plot=R2C2_final_original_coord_2kb_up_down, width = 8, height = 10, units = "in")
system("open ~/Desktop/R2C2_metaplot/original_2kb_up_down.pdf")


###############################3###############################3
###############################3###############################3
# Generte Metaplot for Multiple SCALED-REGIONS ----------------------------
###############################3###############################3
###############################3###############################3
#Generate Normalized Matricies for scaled regions
generate_tss_heatmap_scaled_region <- function(TSS_file, mark_name) {
  
  #TSS_file <- TSS_file[, -c(151:200)] # delete columns 5 through 7
  
  TSS <- as.normalizedMatrix(TSS_file, 
                             k_upstream = 400, 
                             k_downstream = 400, 
                             k_target = 200,
                             extend = c(400, 400), 
                             signal_name = mark_name, 
                             target_name = c("TSS","TTS"),
                             keep = c(0,.99), smooth = FALSE)
  
  
  
  return(TSS)
}
generate_tss_heatmap_r2c2_scaled_region <- function(TSS_file, mark_name) {
  
  #TSS_file <- TSS_file[, -c(151:200)] # delete columns 5 through 7
  
  TSS <- as.normalizedMatrix(TSS_file, 
                             k_upstream = 400, 
                             k_downstream = 400, 
                             k_target = 200,
                             extend = c(400, 400), 
                             background = 0, 
                             signal_name = mark_name, 
                             target_name = c("TSS","TTS"),
                             keep = c(0,.99), smooth = FALSE)
  
  
  
  return(TSS)
}

#Generate the plots for scaled region matricies 
generate_complex_heatmaps_chip_scaled_region <- function(TSS_matrix, mark_name, color_hex) {
  
  #If row title is NOT missing - pass on it and don't add it 
  
  TSS_matrix <- TSS_matrix
  
  
  #Scale colors across matricies the same
  #common_min <-  min(TSS_matrix)
  common_max <-  max(TSS_matrix)
  col_fun <- circlize::colorRamp2(c(0, common_max), c("white", color_hex))
  
  #Get the same max value for each 
  #y_max_val <- round(max(colMeans(TES_matrix, na.rm = TRUE), colMeans(TSS_matrix, na.rm = TRUE))) + 2 
  
  axis_name = c("-2000bp","TSS","TTS","2000bp")
  final_graph <-  EnrichedHeatmap(TSS_matrix, cluster_rows = FALSE, column_title =  mark_name, show_heatmap_legend = FALSE, use_raster = TRUE, column_title_gp = gpar(fontsize = 30, fontface = "bold"), 
                                  col = col_fun, axis_name_rot = 90, row_title_rot = 0, axis_name_gp = gpar(fontsize = 18), axis_name = axis_name,pos_line_gp = gpar(lty = 3),
                                  top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4)))) 
  
  
  
  
  
  return(final_graph)
  
  
}
generate_complex_heatmaps_ATAC_scaled_region <- function(final_matrix_list, mark_name, color_hex) {
  
  #If row title is NOT missing - pass on it and don't add it 
  
  TSS_matrix <- final_matrix_list
  common_max <-  max(TSS_matrix)
  #col_fun <- circlize::colorRamp2(c(0, common_max), c("white", color_hex))
  col_fun = colorRamp2(quantile(TSS_matrix, c(0, 0.95)), c("white", color_hex))
  
  axis_name = c("-2000bp","TSS","TTS","2000bp")
  #This is the only difference between ATAC and ChIP-seq metaplots. Basically we can't scale the same way we did 
  #previously. This allows us to scale the top part of the metaplot appropriatly.
  
  final_graph <-  EnrichedHeatmap(TSS_matrix, cluster_rows = FALSE, column_title =  mark_name, show_heatmap_legend = FALSE, use_raster = TRUE,
                                  column_title_gp = gpar(fontsize = 30, fontface = "bold"), axis_name = axis_name, pos_line_gp = gpar(lty = 3),
                                  col = col_fun, axis_name_rot = 90, row_title_rot = 0, axis_name_gp = gpar(fontsize = 18))
  return(final_graph)
}
generate_complex_heatmaps_R2C2_scaled_region <- function(TSS_matrix, mark_name, color_hex) {
  
  #If row title is NOT missing - pass on it and don't add it 
  
  TSS_matrix <- TSS_matrix
  #common_max <-  max(TSS_matrix)
  #col_fun <- circlize::colorRamp2(c(0, common_max), c("white", color_hex))
  col_fun = colorRamp2(quantile(TSS_matrix, c(0, 0.98)), c("white", color_hex))
  axis_name = c("-2000bp","TSS","TTS","2000bp")
  
  
  final_graph <- EnrichedHeatmap(TSS_matrix, cluster_rows = FALSE, column_title =  mark_name, show_heatmap_legend = FALSE, use_raster = TRUE, 
                                 column_title_gp = gpar(fontsize = 30, fontface = "bold"), axis_name = axis_name,pos_line_gp = gpar(lty = 3),
                                 col = col_fun, axis_name_rot = 90, row_title_rot = 0, axis_name_gp = gpar(fontsize = 18))
  
  
  return(final_graph)
}


#Call sub functions depending on mark_type
generate_final_plot_scaled <- function(file_1, mark_name, color_hex, mark_type) {
  
  if (mark_type == "R2C2"){
    TSS_heatmap <- generate_tss_heatmap_r2c2_scaled_region(file_1, mark_name)
    R2C2_graphed_heatmap <- generate_complex_heatmaps_R2C2_scaled_region(TSS_heatmap, mark_name, color_hex)
    return(R2C2_graphed_heatmap)
    
  } else if (mark_type == "ATAC") {
    ATAC_heatmap <- generate_tss_heatmap_scaled_region(file_1, mark_name)
    ATAC_graphed_heatmap <- generate_complex_heatmaps_ATAC_scaled_region(ATAC_heatmap, mark_name, color_hex)
    return(ATAC_graphed_heatmap)
    
  } else if (mark_type == "CHIP") {
    Chip_heatmap <- generate_tss_heatmap_scaled_region(file_1, mark_name)
    Chip_plot_heatmap <- generate_complex_heatmaps_chip_scaled_region(Chip_heatmap, mark_name, color_hex)
    return(Chip_plot_heatmap)
  } 
}


###############################
#### WORK ON ORIGINAL ANNOTATION
###############################
leaf_genes_leaf_H2AZ <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_genes_leaf_H2AZ_all.gz")
leaf_genes_leaf_H3K36me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_genes_leaf_H3K36me3_all.gz")
leaf_genes_leaf_H3K4me1 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_genes_leaf_H3K4me1_all.gz")
leaf_genes_leaf_H3K4me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_genes_leaf_H3K4me3_all.gz")
leaf_genes_leaf_H3K56ac <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_genes_leaf_H3K56ac_all.gz")
leaf_genes_leaf_input <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_genes_leaf_input_all.gz")
leaf_lncRNAs_leaf_H2AZ <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_lncRNAs_leaf_H2AZ_all.gz")
leaf_lncRNAs_leaf_H3K36me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_lncRNAs_leaf_H3K36me3_all.gz")
leaf_lncRNAs_leaf_H3K4me1 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_lncRNAs_leaf_H3K4me1_all.gz")
leaf_lncRNAs_leaf_H3K4me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_lncRNAs_leaf_H3K4me3_all.gz")
leaf_lncRNAs_leaf_H3K56ac <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_lncRNAs_leaf_H3K56ac_all.gz")
leaf_lncRNAs_leaf_input <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_lncRNAs_leaf_input_all.gz")
leaf_mystery_leaf_H2AZ <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_mystery_leaf_H2AZ_all.gz")
leaf_mystery_leaf_H3K36me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_mystery_leaf_H3K36me3_all.gz")
leaf_mystery_leaf_H3K4me1 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_mystery_leaf_H3K4me1_all.gz")
leaf_mystery_leaf_H3K4me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_mystery_leaf_H3K4me3_all.gz")
leaf_mystery_leaf_H3K56ac <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_mystery_leaf_H3K56ac_all.gz")
leaf_mystery_leaf_input <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_mystery_leaf_input_all.gz")



#replace all NAs with a zero
leaf_genes_leaf_H2AZ[is.nan(leaf_genes_leaf_H2AZ)] = 0 
leaf_genes_leaf_H3K36me3[is.nan(leaf_genes_leaf_H3K36me3)] = 0 
leaf_genes_leaf_H3K4me1[is.nan(leaf_genes_leaf_H3K4me1)] = 0 
leaf_genes_leaf_H3K4me3[is.nan(leaf_genes_leaf_H3K4me3)] = 0 
leaf_genes_leaf_H3K56ac[is.nan(leaf_genes_leaf_H3K56ac)] = 0 
leaf_genes_leaf_input[is.nan(leaf_genes_leaf_input)] = 0 
leaf_lncRNAs_leaf_H2AZ[is.nan(leaf_lncRNAs_leaf_H2AZ)] = 0 
leaf_lncRNAs_leaf_H3K36me3[is.nan(leaf_lncRNAs_leaf_H3K36me3)] = 0 
leaf_lncRNAs_leaf_H3K4me1[is.nan(leaf_lncRNAs_leaf_H3K4me1)] = 0 
leaf_lncRNAs_leaf_H3K4me3[is.nan(leaf_lncRNAs_leaf_H3K4me3)] = 0 
leaf_lncRNAs_leaf_H3K56ac[is.nan(leaf_lncRNAs_leaf_H3K56ac)] = 0 
leaf_lncRNAs_leaf_input[is.nan(leaf_lncRNAs_leaf_input)] = 0 
leaf_mystery_leaf_H2AZ[is.nan(leaf_mystery_leaf_H2AZ)] = 0 
leaf_mystery_leaf_H3K36me3[is.nan(leaf_mystery_leaf_H3K36me3)] = 0 
leaf_mystery_leaf_H3K4me1[is.nan(leaf_mystery_leaf_H3K4me1)] = 0 
leaf_mystery_leaf_H3K4me3[is.nan(leaf_mystery_leaf_H3K4me3)] = 0 
leaf_mystery_leaf_H3K56ac[is.nan(leaf_mystery_leaf_H3K56ac)] = 0 
leaf_mystery_leaf_input[is.nan(leaf_mystery_leaf_input)] = 0 




normalize_leaf_gene_leaf_H3K36me3_TSS_row_counts <- nrow(normalized_leaf_genes_leaf_H3K4me1_TES_file)
normalize_leaf_lncRNAs_leaf_H3K36me3_TSS_row_counts <- nrow(normalized_leaf_lncRNAs_leaf_H3K4me1_TES_file)
normalize_leaf_mystery_leaf_H3K36me3_TSS_row_counts <- nrow(normalized_leaf_mystery_leaf_H3K4me1_TSS_file)

values <- c( "lncRNA", "Unannotated", "Gene")
final_split <- rep(values, times = c(normalize_leaf_lncRNAs_leaf_H3K36me3_TSS_row_counts, normalize_leaf_mystery_leaf_H3K36me3_TSS_row_counts, normalize_leaf_gene_leaf_H3K36me3_TSS_row_counts))


# combine_like_marks ------------------------------------------------------

combine_identical_marks_final <- function(genes, lncRNA, unannotated) {
  
  #genes_1 <- map_df(genes) 
  #lncRNA_1 <- map_df(lncRNA)
  #unannotated_1 <- map_df(unannotated)
  
  bind_all_rows <- bind_rows(
    as.data.frame(genes),
    as.data.frame(lncRNA),
    as.data.frame(unannotated))
  
  final_return <- as.matrix(bind_all_rows)
  
  return(final_return)
}






leaf_H3K56ac_R2C2_plot_sclaed_region_og <- generate_final_plot_scaled(leaf_H3K56ac_R2C2_coord_og, "H3K56ac", "#B89BC9", "CHIP")
leaf_H3K4me3_R2C2_plot_sclaed_region_og <- generate_final_plot_scaled(leaf_H3K4me3_R2C2_coord_og, "H3K4me3", "#99CA3C", "CHIP")
leaf_R2C2_R2C2_plot_sclaed_region_og <- generate_final_plot_scaled(leaf_R2C2_R2C2_coord_og, "R2C2", "#FFCD05", "R2C2")
leaf_ATAC_R2C2_plot_sclaed_region_og <- generate_final_plot_scaled(leaf_ATAC_R2C2_coord_og, "ATAC", "#66c2a5", "ATAC")



final_graph_scale_graph_original_coord <- (leaf_R2C2_R2C2_plot_sclaed_region_og + leaf_ATAC_R2C2_plot_sclaed_region_og + 
                                             leaf_H3K4me3_R2C2_plot_sclaed_region_og + leaf_H3K56ac_R2C2_plot_sclaed_region_og)


units_apart <- c(.7)
draw_units <- rep(units_apart, 4)
R2C2_final_scaled_graph_original_coord = grid.grabExpr(draw(final_graph_scale_graph_original_coord, ht_gap = unit(draw_units, "cm")))
ggsave(path = "~/Desktop/R2C2_metaplot", filename="original_2kb_up_down_scaled_region.pdf", plot=R2C2_final_scaled_graph_original_coord, width = 8, height = 10, units = "in")
system("open ~/Desktop/R2C2_metaplot/original_2kb_up_down_scaled_region.pdf")

