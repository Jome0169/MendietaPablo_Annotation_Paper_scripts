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
    select(-X1,-X2,-X3,-X5,-X6) %>% 
    dplyr::rename(combined_name = X4) %>% 
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
                                  col = col_fun, axis_name_rot = 90, row_title_rot = 90, axis_name_gp = gpar(fontsize = 18), axis_name = axis_name,pos_line_gp = gpar(lty = 3),
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
                                  col = col_fun, axis_name_rot = 90, row_title_rot = 90, axis_name_gp = gpar(fontsize = 18))
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
                                 col = col_fun, axis_name_rot = 90, row_title_rot = 90, axis_name_gp = gpar(fontsize = 18))
  
  
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
ggsave(path = "~/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/R2C2_metaplot", filename="original_2kb_up_down.pdf", plot=R2C2_final_original_coord_2kb_up_down, width = 8, height = 10, units = "in")
system("open ~/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/R2C2_metaplot/original_2kb_up_down.pdf")







######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################
###############################3###############################3
###############################3###############################3
# Generte Metaplot for Multiple SCALED-REGIONS ----------------------------
###############################3###############################3
###############################3###############################3
#Generate Normalized Matricies for scaled regions
generate_tss_heatmap_scaled_region <- function(TSS_file, mark_name) {
  
  #TSS_file <- TSS_file[, -c(151:200)] # delete columns 5 through 7
  
  TSS <- as.normalizedMatrix(TSS_file, 
                             k_upstream = 50, 
                             k_downstream = 50, 
                             k_target = 250,
                             extend = c(50, 50), 
                             signal_name = mark_name, 
                             target_name = c("TSS","TTS"),
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
  
  axis_name = c("-1000bp","TSS","TTS","1000bp")
  final_graph <-  EnrichedHeatmap(TSS_matrix,  column_title =  mark_name, use_raster = TRUE, show_heatmap_legend = FALSE,
                                  row_split = row_split_array, cluster_rows = FALSE, 
                                  column_title_gp = gpar(fontsize = 25, fontface = "bold"), 
                                  col = col_fun, axis_name_rot = 90, row_title_rot = 90, axis_name_gp = gpar(fontsize = 18), axis_name = axis_name,pos_line_gp = gpar(lty = 3),
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
  
  final_graph <-  EnrichedHeatmap(TSS_matrix,  column_title =  mark_name, use_raster = TRUE, show_heatmap_legend = FALSE,
                                  row_split = row_split_array, cluster_rows = FALSE, 
                                  column_title_gp = gpar(fontsize = 25, fontface = "bold"), axis_name = axis_name, pos_line_gp = gpar(lty = 3),
                                  col = col_fun, axis_name_rot = 90, row_title_rot = 90, axis_name_gp = gpar(fontsize = 18))
  return(final_graph)
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



leaf_gene_leaf <- nrow(leaf_genes_leaf_H2AZ)
leaf_lncRNAs_leaf <- nrow(leaf_lncRNAs_leaf_H2AZ)
leaf_mystery_leaf <- nrow(leaf_mystery_leaf_H2AZ)
values <- c( "Gene", "lncRNA", "Unannotated")
final_split <- rep(values, times = c(leaf_gene_leaf, leaf_lncRNAs_leaf, leaf_mystery_leaf))


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

H2AZ_combined_scaled <- combine_identical_marks_final(leaf_genes_leaf_H2AZ, leaf_lncRNAs_leaf_H2AZ, leaf_mystery_leaf_H2AZ)
H3K36me3_combined_scaled <- combine_identical_marks_final(leaf_genes_leaf_H3K36me3, leaf_lncRNAs_leaf_H3K36me3, leaf_mystery_leaf_H3K36me3)
H3K4me1_combined_scaled <- combine_identical_marks_final(leaf_genes_leaf_H3K4me1, leaf_lncRNAs_leaf_H3K4me1, leaf_mystery_leaf_H3K4me1)
H3K4me3_combined_scaled <- combine_identical_marks_final(leaf_genes_leaf_H3K4me3, leaf_lncRNAs_leaf_H3K4me3, leaf_mystery_leaf_H3K4me3)
H3K56ac_combined_scaled <- combine_identical_marks_final(leaf_genes_leaf_H3K56ac, leaf_lncRNAs_leaf_H3K56ac, leaf_mystery_leaf_H3K56ac)
input_combined_scaled <- combine_identical_marks_final(leaf_genes_leaf_input, leaf_lncRNAs_leaf_input, leaf_mystery_leaf_input)




H2AZ_combined_scaled_normalized <- (H2AZ_combined_scaled + 1) / (input_combined_scaled + 1)
H3K36me3_combined_scaled_normalized <- (H3K36me3_combined_scaled + 1)/ (input_combined_scaled + 1)
H3K4me1_combined_scaled_normalized <- (H3K4me1_combined_scaled + 1)/ (input_combined_scaled + 1)
H3K4me3_combined_scaled_normalized <- (H3K4me3_combined_scaled + 1)/ (input_combined_scaled + 1)
H3K56ac_combined_scaled_normalized <- (H3K56ac_combined_scaled + 1)/ (input_combined_scaled + 1)


H3K36me3_colors <- "#CB2026"
H3K4me1_colors <-  "#642165"
H3K4me3_colors <- "#B89BC9"
H3K56ac_colors <- "#99CA3C"
ATAC_colors <- "#FFCD05"
H2A_colors <-  "#66c2a5"


H2AZ_combined_scaled_normalized_plots <- generate_final_plot_scaled(H2AZ_combined_scaled_normalized, "H2Az", H2A_colors, final_split, "CHIP")
H3K36me3_combined_scaled_normalized_plots <- generate_final_plot_scaled(H3K36me3_combined_scaled_normalized, 'H3K36me3', H3K36me3_colors, final_split, "CHIP")
H3K4me1_combined_scaled_normalized_plots <- generate_final_plot_scaled(H3K4me1_combined_scaled_normalized, "H3K4me1", H3K4me1_colors, final_split, "CHIP")
H3K4me3_combined_scaled_normalized_plots <- generate_final_plot_scaled(H3K4me3_combined_scaled_normalized, "H3K4me3", H3K4me3_colors, final_split, "CHIP")
H3K56ac_combined_scaled_normalized_plots <- generate_final_plot_scaled(H3K56ac_combined_scaled_normalized, "H3K56ac", H3K56ac_colors, final_split, "CHIP")


leaf_ATAC_R2C2_plot_sclaed_region_og <- generate_final_plot_scaled(leaf_ATAC_R2C2_coord_og, "ATAC", ATAC_colors, "ATAC")


test_final <- (H3K36me3_combined_scaled_normalized_plots +  H3K4me1_combined_scaled_normalized_plots +
H3K4me3_combined_scaled_normalized_plots + 
H3K56ac_combined_scaled_normalized_plots)


units_apart <- c(.7)
draw_units <- rep(units_apart, 3)
R2C2_final_scaled_graph_original_coord = grid.grabExpr(draw(test_final, ht_gap = unit(draw_units, "cm")))
ggsave(path = "~/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/", filename="scalred_region_test.pdf", plot=R2C2_final_scaled_graph_original_coord, width = 8, height = 10, units = "in")
system("open ~/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/scalred_region_test.pdf")




#################################
#### Generate Metplots for Genes which are ON and OFF
################################ 
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

read_in_TPM_file <- function(TPM_file) {
  col_names_TPM_file <- c("gene_ID", "TPM")
  read_TPM_gene <- read_delim(TPM_file, '\t', col_names = col_names_TPM_file) %>% 
    mutate(updated_name = str_c("gene", gene_ID , sep =':')) %>% 
    select(-gene_ID) %>% 
    dplyr::rename(gene_ID = updated_name)
  
  return(read_TPM_gene)
}
generate_ordered_gene_list <- function(loaded_TPM_file, gene_set_on_off) {
  
  
  gather_gene_row_names <- row.names(gene_set_on_off)
  
  filtered_gene_TPMs <- loaded_TPM_file %>% 
    filter(gene_ID %in% gather_gene_row_names) %>% 
    mutate(log2_tpm = log2(TPM +1)) %>% 
    dplyr::arrange(desc(log2_tpm))
    
  sorted_expression <- filtered_gene_TPMs$gene_ID
  return(sorted_expression)
  
  
}


leaf_TPM_by_gene <- read_in_TPM_file("/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/00.data/TPM_vals/RNA_B73_leaf_TPM.surviving.txt")

leaf_off_genes_leaf_H2AZ <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_off_genes_leaf_H2AZ_all.gz")
leaf_off_genes_leaf_H3K27me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_off_genes_leaf_H3K27me3_all.gz")
leaf_off_genes_leaf_H3K36me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_off_genes_leaf_H3K36me3_all.gz")
leaf_off_genes_leaf_H3K4me1 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_off_genes_leaf_H3K4me1_all.gz")
leaf_off_genes_leaf_H3K4me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_off_genes_leaf_H3K4me3_all.gz")
leaf_off_genes_leaf_H3K56ac <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_off_genes_leaf_H3K56ac_all.gz")
leaf_off_genes_leaf_input <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_off_genes_leaf_input_all.gz")
leaf_off_genes_leaf_ATAC_scaled <- read_gunzipped_file("07.scaled_to_region/ATAC/leaf_off_genes_leaf_ATAC_scaled_all.gz")

leaf_on_genes_leaf_H2AZ <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_on_genes_leaf_H2AZ_all.gz")
leaf_on_genes_leaf_H3K27me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_on_genes_leaf_H3K27me3_all.gz")
leaf_on_genes_leaf_H3K36me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_on_genes_leaf_H3K36me3_all.gz")
leaf_on_genes_leaf_H3K4me1 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_on_genes_leaf_H3K4me1_all.gz")
leaf_on_genes_leaf_H3K4me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_on_genes_leaf_H3K4me3_all.gz")
leaf_on_genes_leaf_H3K56ac <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_on_genes_leaf_H3K56ac_all.gz")
leaf_on_genes_leaf_input <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_on_genes_leaf_input_all.gz")
leaf_on_genes_leaf_ATAC_scaled <- read_gunzipped_file("07.scaled_to_region/ATAC/leaf_on_genes_leaf_ATAC_scaled_all.gz")



###################3

number_off <- nrow(leaf_off_genes_leaf_H2AZ)
number_on <- nrow(leaf_on_genes_leaf_H2AZ)
values <- c( "Expressed Genes", "Non-Expressed Genes")
final_split <- rep(values, times = c(number_on, number_off))
leaf_sorted_on_genes <- generate_ordered_gene_list(leaf_TPM_by_gene, leaf_on_genes_leaf_H2AZ)
leaf_sorted_off_genes <- generate_ordered_gene_list(leaf_TPM_by_gene, leaf_off_genes_leaf_H2AZ)



leaf_H2AZ <- combine_identical_marks_on_off(leaf_on_genes_leaf_H2AZ, leaf_off_genes_leaf_H2AZ)
leaf_H3K27me3 <- combine_identical_marks_on_off(leaf_on_genes_leaf_H3K27me3, leaf_off_genes_leaf_H3K27me3)
leaf_H3K36me3 <- combine_identical_marks_on_off(leaf_on_genes_leaf_H3K36me3, leaf_off_genes_leaf_H3K36me3)
leaf_H3K4me1 <- combine_identical_marks_on_off(leaf_on_genes_leaf_H3K4me1, leaf_off_genes_leaf_H3K4me1)
leaf_H3K4me3 <- combine_identical_marks_on_off(leaf_on_genes_leaf_H3K4me3, leaf_off_genes_leaf_H3K4me3)
leaf_H3K56ac <- combine_identical_marks_on_off(leaf_on_genes_leaf_H3K56ac, leaf_off_genes_leaf_H3K56ac)
leaf_input <- combine_identical_marks_on_off(leaf_on_genes_leaf_input, leaf_off_genes_leaf_input)
leaf_ATAC_scaled <- combine_identical_marks_on_off(leaf_on_genes_leaf_ATAC_scaled, leaf_off_genes_leaf_ATAC_scaled)


leaf_H2AZ_normalized <- (leaf_H2AZ + 1)/(leaf_input +1)
leaf_H3K27me3_normalized <- (leaf_H3K27me3 + 1)/(leaf_input +1)
leaf_H3K36me3_normalized <- (leaf_H3K36me3 + 1)/(leaf_input +1)
leaf_H3K4me1_normalized <- (leaf_H3K4me1 + 1)/(leaf_input +1)
leaf_H3K4me3_normalized <- (leaf_H3K4me3 + 1)/(leaf_input +1)
leaf_H3K56ac_normalized <- (leaf_H3K56ac + 1)/(leaf_input +1)
leaf_ATAC_scaled_normalized <- (leaf_ATAC_scaled + 1)/(leaf_input +1)



H3K36me3_colors <- "#CB2026"
H3K4me1_colors <-  "#642165"
H3K4me3_colors <- "#B89BC9"
H3K56ac_colors <- "#99CA3C"
ATAC_colors <- "#FFCD05"
H2A_colors <-  "#66c2a5"
H3K27me3_colors <- "#234bb8"



leaf_H2AZ_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_H2AZ_normalized, "H2Az", H2A_colors, final_split, "CHIP")
leaf_H3K36me3_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_H3K36me3_normalized, 'H3K36me3', H3K36me3_colors, final_split, "CHIP")
leaf_H3K4me1_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_H3K4me1_normalized, "H3K4me1", H3K4me1_colors, final_split, "CHIP")
leaf_H3K4me3_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_H3K4me3_normalized, "H3K4me3", H3K4me3_colors, final_split, "CHIP")
leaf_H3K56ac_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_H3K56ac_normalized, "H3K56ac", H3K56ac_colors, final_split, "CHIP")
leaf_H3K27me3_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_H3K27me3_normalized, "H3K27me3", H3K27me3_colors, final_split, "CHIP")
leaf_ATAC_R2C2_plot_sclaed_region_og <- generate_final_plot_scaled(leaf_ATAC_scaled_normalized, "ATAC", ATAC_colors, final_split, "ATAC")


leaf_final <- (leaf_H3K36me3_combined_scaled_normalized_plots +
  leaf_H3K4me1_combined_scaled_normalized_plots +
  leaf_H3K4me3_combined_scaled_normalized_plots +
  leaf_H3K56ac_combined_scaled_normalized_plots +
  leaf_H3K27me3_combined_scaled_normalized_plots +
  leaf_H2AZ_combined_scaled_normalized_plots + 
  leaf_ATAC_R2C2_plot_sclaed_region_og)


units_apart <- c(.7)
draw_units <- rep(units_apart, 6)

leaf_final_on_off = grid.grabExpr(draw(leaf_final, ht_gap = unit(draw_units, "cm"), row_order = order(c(leaf_sorted_on_genes, leaf_sorted_off_genes))))
ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/", filename="leaf_final_on_off.pdf", plot=leaf_final_on_off, width = 15.5, height = 10, units = "in")
system("open /Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/leaf_final_on_off.pdf")




# root --------------------------------------------------------------------


root_TPM_by_gene <- read_in_TPM_file("/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/00.data/TPM_vals/RNA_B73_root_TPM.surviving.txt")

root_off_genes_root_H2AZ <- read_gunzipped_file("07.scaled_to_region/CHIP/root_off_genes_root_H2AZ_all.gz")
root_off_genes_root_H3K27me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/root_off_genes_root_H3K27me3_all.gz")
root_off_genes_root_H3K36me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/root_off_genes_root_H3K36me3_all.gz")
root_off_genes_root_H3K4me1 <- read_gunzipped_file("07.scaled_to_region/CHIP/root_off_genes_root_H3K4me1_all.gz")
root_off_genes_root_H3K4me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/root_off_genes_root_H3K4me3_all.gz")
root_off_genes_root_H3K56ac <- read_gunzipped_file("07.scaled_to_region/CHIP/root_off_genes_root_H3K56ac_all.gz")
root_off_genes_root_input <- read_gunzipped_file("07.scaled_to_region/CHIP/root_off_genes_root_input_all.gz")
root_off_genes_root_ATAC_scaled <- read_gunzipped_file("07.scaled_to_region/ATAC/root_off_genes_root_ATAC_scaled_all.gz")

root_on_genes_root_H2AZ <- read_gunzipped_file("07.scaled_to_region/CHIP/root_on_genes_root_H2AZ_all.gz")
root_on_genes_root_H3K27me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/root_on_genes_root_H3K27me3_all.gz")
root_on_genes_root_H3K36me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/root_on_genes_root_H3K36me3_all.gz")
root_on_genes_root_H3K4me1 <- read_gunzipped_file("07.scaled_to_region/CHIP/root_on_genes_root_H3K4me1_all.gz")
root_on_genes_root_H3K4me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/root_on_genes_root_H3K4me3_all.gz")
root_on_genes_root_H3K56ac <- read_gunzipped_file("07.scaled_to_region/CHIP/root_on_genes_root_H3K56ac_all.gz")
root_on_genes_root_input <- read_gunzipped_file("07.scaled_to_region/CHIP/root_on_genes_root_input_all.gz")
root_on_genes_root_ATAC_scaled <- read_gunzipped_file("07.scaled_to_region/ATAC/root_on_genes_root_ATAC_scaled_all.gz")


number_off <- nrow(root_off_genes_root_H2AZ)
number_on <- nrow(root_on_genes_root_H2AZ)
values <- c( "Expressed Genes", "Non-Expressed Genes")
final_split <- rep(values, times = c(number_off, number_on))

root_sorted_on_genes <- generate_ordered_gene_list(root_TPM_by_gene, root_on_genes_root_H2AZ)
root_sorted_off_genes <- generate_ordered_gene_list(root_TPM_by_gene, root_off_genes_root_H2AZ)




root_H2AZ <- combine_identical_marks_on_off(root_on_genes_root_H2AZ, root_off_genes_root_H2AZ)
root_H3K27me3 <- combine_identical_marks_on_off(root_on_genes_root_H3K27me3, root_off_genes_root_H3K27me3)
root_H3K36me3 <- combine_identical_marks_on_off(root_on_genes_root_H3K36me3, root_off_genes_root_H3K36me3)
root_H3K4me1 <- combine_identical_marks_on_off(root_on_genes_root_H3K4me1, root_off_genes_root_H3K4me1)
root_H3K4me3 <- combine_identical_marks_on_off(root_on_genes_root_H3K4me3, root_off_genes_root_H3K4me3)
root_H3K56ac <- combine_identical_marks_on_off(root_on_genes_root_H3K56ac, root_off_genes_root_H3K56ac)
root_input <- combine_identical_marks_on_off(root_on_genes_root_input, root_off_genes_root_input)
root_ATAC_scaled <- combine_identical_marks_on_off(root_on_genes_root_ATAC_scaled, root_off_genes_root_ATAC_scaled)


root_H2AZ_normalized <- (root_H2AZ + 1)/(root_input +1)
root_H3K27me3_normalized <- (root_H3K27me3 + 1)/(root_input +1)
root_H3K36me3_normalized <- (root_H3K36me3 + 1)/(root_input +1)
root_H3K4me1_normalized <- (root_H3K4me1 + 1)/(root_input +1)
root_H3K4me3_normalized <- (root_H3K4me3 + 1)/(root_input +1)
root_H3K56ac_normalized <- (root_H3K56ac + 1)/(root_input +1)
root_ATAC_scaled_normalized <- (root_ATAC_scaled + 1)/(root_input +1)

H3K36me3_colors <- "#CB2026"
H3K4me1_colors <-  "#642165"
H3K4me3_colors <- "#B89BC9"
H3K56ac_colors <- "#99CA3C"
ATAC_colors <- "#FFCD05"
H2A_colors <-  "#66c2a5"
H3K27me3_colors <- "#234bb8"

root_H2AZ_combined_scaled_normalized_plots <- generate_final_plot_scaled(root_H2AZ_normalized, "H2Az", H2A_colors, final_split, "CHIP")
root_H3K36me3_combined_scaled_normalized_plots <- generate_final_plot_scaled(root_H3K36me3_normalized, 'H3K36me3', H3K36me3_colors, final_split, "CHIP")
root_H3K4me1_combined_scaled_normalized_plots <- generate_final_plot_scaled(root_H3K4me1_normalized, "H3K4me1", H3K4me1_colors, final_split, "CHIP")
root_H3K4me3_combined_scaled_normalized_plots <- generate_final_plot_scaled(root_H3K4me3_normalized, "H3K4me3", H3K4me3_colors, final_split, "CHIP")
root_H3K56ac_combined_scaled_normalized_plots <- generate_final_plot_scaled(root_H3K56ac_normalized, "H3K56ac", H3K56ac_colors, final_split, "CHIP")
root_H3K27me3_combined_scaled_normalized_plots <- generate_final_plot_scaled(root_H3K27me3_normalized, "H3K27me3", H3K27me3_colors, final_split, "CHIP")
root_ATAC_R2C2_plot_sclaed_region_og <- generate_final_plot_scaled(root_ATAC_scaled_normalized, "ATAC", ATAC_colors, final_split, "ATAC")


root_final <- (root_H3K36me3_combined_scaled_normalized_plots +
                 root_H3K4me1_combined_scaled_normalized_plots +
                 root_H3K4me3_combined_scaled_normalized_plots +
                 root_H3K56ac_combined_scaled_normalized_plots +
                 root_H3K27me3_combined_scaled_normalized_plots +
                 root_H2AZ_combined_scaled_normalized_plots + 
                 root_ATAC_R2C2_plot_sclaed_region_og)

units_apart <- c(.7)
draw_units <- rep(units_apart, 6)


root_final_on_off = grid.grabExpr(draw(root_final, ht_gap = unit(draw_units, "cm"), row_order = order(c(root_sorted_on_genes, root_sorted_off_genes))))
ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/", filename="root_final_on_off.pdf", plot=root_final_on_off, width = 15.5, height = 10, units = "in")
system("open /Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/root_final_on_off.pdf")


# ear ---------------------------------------------------------------------
ear_off_genes_ear_H2AZ <- read_gunzipped_file("07.scaled_to_region/CHIP/ear_off_genes_ear_H2AZ_all.gz")
ear_off_genes_ear_H3K27me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/ear_off_genes_ear_H3K27me3_all.gz")
ear_off_genes_ear_H3K36me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/ear_off_genes_ear_H3K36me3_all.gz")
ear_off_genes_ear_H3K4me1 <- read_gunzipped_file("07.scaled_to_region/CHIP/ear_off_genes_ear_H3K4me1_all.gz")
ear_off_genes_ear_H3K4me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/ear_off_genes_ear_H3K4me3_all.gz")
ear_off_genes_ear_H3K56ac <- read_gunzipped_file("07.scaled_to_region/CHIP/ear_off_genes_ear_H3K56ac_all.gz")
ear_off_genes_ear_input <- read_gunzipped_file("07.scaled_to_region/CHIP/ear_off_genes_ear_input_all.gz")
ear_off_genes_ear_ATAC_scaled <- read_gunzipped_file("07.scaled_to_region/ATAC/ear_off_genes_ear_ATAC_scaled_all.gz")

ear_on_genes_ear_H2AZ <- read_gunzipped_file("07.scaled_to_region/CHIP/ear_on_genes_ear_H2AZ_all.gz")
ear_on_genes_ear_H3K27me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/ear_on_genes_ear_H3K27me3_all.gz")
ear_on_genes_ear_H3K36me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/ear_on_genes_ear_H3K36me3_all.gz")
ear_on_genes_ear_H3K4me1 <- read_gunzipped_file("07.scaled_to_region/CHIP/ear_on_genes_ear_H3K4me1_all.gz")
ear_on_genes_ear_H3K4me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/ear_on_genes_ear_H3K4me3_all.gz")
ear_on_genes_ear_H3K56ac <- read_gunzipped_file("07.scaled_to_region/CHIP/ear_on_genes_ear_H3K56ac_all.gz")
ear_on_genes_ear_input <- read_gunzipped_file("07.scaled_to_region/CHIP/ear_on_genes_ear_input_all.gz")
ear_on_genes_ear_ATAC_scaled <- read_gunzipped_file("07.scaled_to_region/ATAC/ear_on_genes_ear_ATAC_scaled_all.gz")




number_off <- nrow(ear_off_genes_ear_H2AZ)
number_on <- nrow(ear_on_genes_ear_H2AZ)
values <- c( "Expressed Genes", "Non-Expressed Genes")
final_split <- rep(values, times = c(number_off, number_on))
ear_TPM_by_gene <- read_in_TPM_file("/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/00.data/TPM_vals/RNA_B73_ear_TPM.surviving.txt")
ear_sorted_on_genes <- generate_ordered_gene_list(ear_TPM_by_gene, ear_on_genes_ear_H2AZ)
ear_sorted_off_genes <- generate_ordered_gene_list(ear_TPM_by_gene, ear_off_genes_ear_H2AZ)




ear_H2AZ <- combine_identical_marks_on_off(ear_on_genes_ear_H2AZ, ear_off_genes_ear_H2AZ)
ear_H3K27me3 <- combine_identical_marks_on_off(ear_on_genes_ear_H3K27me3, ear_off_genes_ear_H3K27me3)
ear_H3K36me3 <- combine_identical_marks_on_off(ear_on_genes_ear_H3K36me3, ear_off_genes_ear_H3K36me3)
ear_H3K4me1 <- combine_identical_marks_on_off(ear_on_genes_ear_H3K4me1, ear_off_genes_ear_H3K4me1)
ear_H3K4me3 <- combine_identical_marks_on_off(ear_on_genes_ear_H3K4me3, ear_off_genes_ear_H3K4me3)
ear_H3K56ac <- combine_identical_marks_on_off(ear_on_genes_ear_H3K56ac, ear_off_genes_ear_H3K56ac)
ear_input <- combine_identical_marks_on_off(ear_on_genes_ear_input, ear_off_genes_ear_input)
ear_ATAC_scaled <- combine_identical_marks_on_off(ear_on_genes_ear_ATAC_scaled, ear_off_genes_ear_ATAC_scaled)


ear_H2AZ_normalized <- (ear_H2AZ + 1)/(ear_input +1)
ear_H3K27me3_normalized <- (ear_H3K27me3 + 1)/(ear_input +1)
ear_H3K36me3_normalized <- (ear_H3K36me3 + 1)/(ear_input +1)
ear_H3K4me1_normalized <- (ear_H3K4me1 + 1)/(ear_input +1)
ear_H3K4me3_normalized <- (ear_H3K4me3 + 1)/(ear_input +1)
ear_H3K56ac_normalized <- (ear_H3K56ac + 1)/(ear_input +1)
ear_ATAC_scaled_normalized <- (ear_ATAC_scaled + 1)/(ear_input +1)

H3K36me3_colors <- "#CB2026"
H3K4me1_colors <-  "#642165"
H3K4me3_colors <- "#B89BC9"
H3K56ac_colors <- "#99CA3C"
ATAC_colors <- "#FFCD05"
H2A_colors <-  "#66c2a5"
H3K27me3_colors <- "#234bb8"

ear_H2AZ_combined_scaled_normalized_plots <- generate_final_plot_scaled(ear_H2AZ_normalized, "H2Az", H2A_colors, final_split, "CHIP")
ear_H3K36me3_combined_scaled_normalized_plots <- generate_final_plot_scaled(ear_H3K36me3_normalized, 'H3K36me3', H3K36me3_colors, final_split, "CHIP")
ear_H3K4me1_combined_scaled_normalized_plots <- generate_final_plot_scaled(ear_H3K4me1_normalized, "H3K4me1", H3K4me1_colors, final_split, "CHIP")
ear_H3K4me3_combined_scaled_normalized_plots <- generate_final_plot_scaled(ear_H3K4me3_normalized, "H3K4me3", H3K4me3_colors, final_split, "CHIP")
ear_H3K56ac_combined_scaled_normalized_plots <- generate_final_plot_scaled(ear_H3K56ac_normalized, "H3K56ac", H3K56ac_colors, final_split, "CHIP")
ear_H3K27me3_combined_scaled_normalized_plots <- generate_final_plot_scaled(ear_H3K27me3_normalized, "H3K27me3", H3K27me3_colors, final_split, "CHIP")
ear_ATAC_R2C2_plot_sclaed_region_og <- generate_final_plot_scaled(ear_ATAC_scaled_normalized, "ATAC", ATAC_colors, final_split, "ATAC")


ear_final <- (ear_H3K36me3_combined_scaled_normalized_plots +
                 ear_H3K4me1_combined_scaled_normalized_plots +
                 ear_H3K4me3_combined_scaled_normalized_plots +
                 ear_H3K56ac_combined_scaled_normalized_plots +
                 ear_H3K27me3_combined_scaled_normalized_plots +
                 ear_H2AZ_combined_scaled_normalized_plots + 
                 ear_ATAC_R2C2_plot_sclaed_region_og)


units_apart <- c(.7)
draw_units <- rep(units_apart, 6)


ear_final_on_off = grid.grabExpr(draw(ear_final, ht_gap = unit(draw_units, "cm"), row_order = order(c(ear_sorted_on_genes, ear_sorted_off_genes))))
ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/", filename="ear_final_on_off.pdf", plot=ear_final_on_off, width = 15.5, height = 10, units = "in")
system("open /Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/ear_final_on_off.pdf")



