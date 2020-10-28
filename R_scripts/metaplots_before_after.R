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
    select(-X1,-X2,-X3,-X5,-X6) %>% 
    dplyr::rename(combined_name = X4) %>% 
    column_to_rownames(var="combined_name")
  
  return_matrix <- as.matrix(final_return)
  
  return(return_matrix)
  
  
}


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
                             keep = c(0,.95), smooth = TRUE)
  
  
  
  return(TSS)
}

#Generate the plots for scaled region matricies 
generate_complex_heatmaps_chip_scaled_region <- function(TSS_matrix, mark_name, row_split_array, color_hex) {
  
  #If row title is NOT missing - pass on it and don't add it 
  
  TSS_matrix <- TSS_matrix
  TSS_matrix[TSS_matrix < 0] = 0
  quick <- TSS_matrix[TSS_matrix > 0]
  
  
  #Scale colors across matricies the same
  #col_fun = colorRamp2(quantile(TSS_matrix, c(0, .9)), c("white", color_hex))
  #col_fun = circlize::colorRamp2(c(common_min, common_max), c("white", color_hex))
  
  col_fun = colorRamp2(quantile(quick, c(0, .8),  na.rm = TRUE), c("white", color_hex))
  #Get the same max value for each 
  #y_max_val <- round(max(colMeans(TES_matrix, na.rm = TRUE), colMeans(TSS_matrix, na.rm = TRUE))) + 2 
  
  axis_name = c("-1kb","TSS","TTS","1kp")
  
  if (is.null(row_split_array) == FALSE) {
    
  
  
  final_graph <-  EnrichedHeatmap(TSS_matrix,  column_title =  mark_name, use_raster = TRUE, show_heatmap_legend = FALSE,
                                  row_split = row_split_array, cluster_rows = FALSE, 
                                  column_title_gp = gpar(fontsize = 25, fontface = "bold"), 
                                  col = col_fun, row_title_rot = 90, axis_name_gp = gpar(fontsize = 18), axis_name = axis_name,
                                  pos_line_gp = gpar(lty = 3),
                                  top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = c(H3K27me3_colors, "#888A8D"), lwd = 3 )))) 
  } else if (is.null(row_split_array) == TRUE) {
    final_graph <-  EnrichedHeatmap(TSS_matrix,  column_title =  mark_name, use_raster = TRUE, show_heatmap_legend = FALSE,
                                    cluster_rows = FALSE, 
                                    column_title_gp = gpar(fontsize = 25, fontface = "bold"), 
                                    col = col_fun, row_title_rot = 90, axis_name_gp = gpar(fontsize = 18), axis_name = axis_name,
                                    pos_line_gp = gpar(lty = 3),
                                    top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = c(H3K27me3_colors, "#888A8D"), lwd = 3 )))) 
    
    
    
    
  }
  return(final_graph)
  
}

generate_complex_heatmaps_ATAC_scaled_region <- function(final_matrix_list, mark_name, row_split_array, color_hex) {
  
  #If row title is NOT missing - pass on it and don't add it 
  
  TSS_matrix <- final_matrix_list
  TSS_matrix[TSS_matrix < 0] = 0
  #col_fun <- circlize::colorRamp2(c(0, common_max), c("white", color_hex))
  #col_fun = colorRamp2(quantile(TSS_matrix, c(0, .9)), c("white", color_hex))
  col_fun = colorRamp2(quantile(TSS_matrix, c(.1, .85)), c("white", color_hex))
  
  axis_name = c("-1000bp","TSS","TTS","1000bp")
  #This is the only difference between ATAC and ChIP-seq metaplots. Basically we can't scale the same way we did 
  #previously. This allows us to scale the top part of the metaplot appropriatly.
  
  final_graph <-  EnrichedHeatmap(TSS_matrix,  column_title =  mark_name, use_raster = TRUE, show_heatmap_legend = FALSE,
                                  row_split = row_split_array, cluster_rows = FALSE, 
                                  column_title_gp = gpar(fontsize = 25, fontface = "bold"), axis_name = axis_name, 
                                  col = col_fun, axis_name_rot = 90, row_title_rot = 90, axis_name_gp = gpar(fontsize = 18),
                                  pos_line_gp = gpar(lty = 3),
                                  top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4)))) 
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




#################################
#### Generate Metplots for Genes which are ON and OFF
################################ 
combine_identical_marks_on_off <- function(genes_on, genes_off) {
  
  #genes_1 <- map_df(genes) 
  #lncRNA_1 <- map_df(lncRNA)
  #unannotated_1 <- map_df(unannotated)
  
  bind_all_rows <- rbind(as.data.frame(genes_on),as.data.frame(genes_off))
  final_return <- as.matrix(bind_all_rows)
  
  return(final_return)
}
read_in_TPM_file <- function(TPM_file) {
  col_names_TPM_file <- c("gene_ID", "TPM")
  read_TPM_gene <- read_delim(TPM_file, '\t', col_names = col_names_TPM_file) %>% 
    mutate(name = str_c("gene", gene_ID , sep =':')) %>% 
    select(-gene_ID) %>% 
    dplyr::rename(gene_ID = name)
  
  return(read_TPM_gene)
}
alternative_TPM_file <- function(TPM_file) {
  col_names_TPM_file <- c("gene_ID", "Chr", "Start", "End", "Length", "Reads", "TPM")
  read_TPM_gene <- read_delim(TPM_file, '\t', col_names = col_names_TPM_file) %>% 
    mutate(name = str_c("gene", gene_ID , sep =':')) %>% 
    select(-gene_ID) %>% 
    dplyr::rename(gene_ID = name) %>% 
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
generate_ordered_gene_list <- function(loaded_TPM_file, gene_set_on_off) {
  
  
  gather_gene_row_names <- row.names(gene_set_on_off)
  
  filtered_gene_TPMs <- loaded_TPM_file %>% 
    filter(gene_ID %in% gather_gene_row_names) %>% 
    mutate(log2_tpm = log2(TPM +1)) %>% 
    dplyr::arrange(desc(TPM))
  
  sorted_expression <- filtered_gene_TPMs$gene_ID
  return(sorted_expression)
  
  
}


leaf_TPM_by_gene <- alternative_TPM_file("/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/00.data/TPM_vals/RNA_B73_leaf_Aligned.sortedByCoord.out_genes.out")
#read_in_TPM_file("/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/00.data/TPM_vals/RNA_B73_leaf_TPM.surviving.txt")




leaf_original_major_extension_genes_leaf_H3K27me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_major_extension_metaplots_original_leaf_H3K27me3_all.gz")
leaf_original_major_extension_genes_leaf_H3K36me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_major_extension_metaplots_original_leaf_H3K36me3_all.gz")
leaf_original_major_extension_genes_leaf_H3K4me1 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_major_extension_metaplots_original_leaf_H3K4me1_all.gz")
leaf_original_major_extension_genes_leaf_H3K4me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_major_extension_metaplots_original_leaf_H3K4me3_all.gz")
leaf_original_major_extension_genes_leaf_H3K56ac <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_major_extension_metaplots_original_leaf_H3K56ac_all.gz")
leaf_original_major_extension_genes_leaf_input <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_major_extension_metaplots_original_leaf_input_all.gz")




leaf_major_extension_genes_leaf_H3K27me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_passing_after_major_extension_genes_leaf_H3K27me3_all.gz")
leaf_major_extension_genes_leaf_H3K36me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_passing_after_major_extension_genes_leaf_H3K36me3_all.gz")
leaf_major_extension_genes_leaf_H3K4me1 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_passing_after_major_extension_genes_leaf_H3K4me1_all.gz")
leaf_major_extension_genes_leaf_H3K4me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_passing_after_major_extension_genes_leaf_H3K4me3_all.gz")
leaf_major_extension_genes_leaf_H3K56ac <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_passing_after_major_extension_genes_leaf_H3K56ac_all.gz")
leaf_major_extension_genes_leaf_input <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_passing_after_major_extension_genes_leaf_input_all.gz")



#leaf_H2AZ <- combine_identical_marks_on_off(leaf_on_genes_leaf_H2AZ, leaf_off_genes_leaf_H2AZ)
major_extension_leaf_H3K27me3 <- combine_identical_marks_on_off(leaf_original_major_extension_genes_leaf_H3K27me3, leaf_major_extension_genes_leaf_H3K27me3)
major_extension_leaf_H3K36me3 <- combine_identical_marks_on_off(leaf_original_major_extension_genes_leaf_H3K36me3, leaf_major_extension_genes_leaf_H3K36me3)
major_extension_leaf_H3K4me1 <- combine_identical_marks_on_off(leaf_original_major_extension_genes_leaf_H3K4me1, leaf_major_extension_genes_leaf_H3K4me1)
major_extension_leaf_H3K4me3 <- combine_identical_marks_on_off(leaf_original_major_extension_genes_leaf_H3K4me3, leaf_major_extension_genes_leaf_H3K4me3)
major_extension_leaf_H3K56ac <- combine_identical_marks_on_off(leaf_original_major_extension_genes_leaf_H3K56ac, leaf_major_extension_genes_leaf_H3K56ac)
major_extension_leaf_input <- combine_identical_marks_on_off(leaf_original_major_extension_genes_leaf_input, leaf_major_extension_genes_leaf_input)
#leaf_ATAC_scaled <- combine_identical_marks_on_off(leaf_on_genes_leaf_ATAC_scaled, leaf_off_genes_leaf_ATAC_scaled)


annotation_N <- nrow(leaf_major_extension_genes_leaf_H3K36me3)
original_annotation_N <- nrow(leaf_original_major_extension_genes_leaf_H3K36me3)
values <- c("Old Annotations", "Updated Annotation")
final_split <- rep(values, times = c(original_annotation_N, annotation_N))
leaf_sorted_gene_list <- generate_ordered_gene_list(leaf_TPM_by_gene, leaf_H3K36me3)



#leaf_major_extension_H2AZ_normalized <- (leaf_H2AZ)-(leaf_input)
leaf_major_extension_H3K27me3_normalized <- (major_extension_leaf_H3K27me3)-(major_extension_leaf_input)
leaf_major_extension_H3K36me3_normalized <- (major_extension_leaf_H3K36me3)-(major_extension_leaf_input)
leaf_major_extension_H3K4me1_normalized <- (major_extension_leaf_H3K4me1)-(major_extension_leaf_input)
leaf_major_extension_H3K4me3_normalized <- (major_extension_leaf_H3K4me3)-(major_extension_leaf_input)
leaf_major_extension_H3K56ac_normalized <- (major_extension_leaf_H3K56ac)-(major_extension_leaf_input)
#leaf_major_extension_ATAC_scaled_normalized <- (leaf_ATAC_scaled)-(leaf_input)



H3K36me3_colors <- "#CB2026"
H3K4me1_colors <-  "#642165"
H3K4me3_colors <- "#B89BC9"
H3K56ac_colors <- "#99CA3C"
ATAC_colors <- "#FFCD05"
H2A_colors <-  "#66c2a5"
H3K27me3_colors <- "#234bb8"



#leaf_H2AZ_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_H2AZ_normalized, "H2Az", H2A_colors, final_split, "CHIP")
leaf_major_extension_H3K36me3_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_major_extension_H3K36me3_normalized, 'H3K36me3', H3K36me3_colors, final_split, "CHIP")
leaf_major_extension_H3K4me1_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_major_extension_H3K4me1_normalized, "H3K4me1", H3K4me1_colors, final_split, "CHIP")
leaf_major_extension_H3K4me3_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_major_extension_H3K4me3_normalized, "H3K4me3", H3K4me3_colors, final_split, "CHIP")
leaf_major_extension_H3K56ac_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_major_extension_H3K56ac_normalized, "H3K56ac", H3K56ac_colors, final_split, "CHIP")
leaf_major_extension_H3K27me3_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_major_extension_H3K27me3_normalized, "H3K27me3", H3K27me3_colors, final_split, "CHIP")
#leaf_ATAC_R2C2_plot_sclaed_region_og <- generate_final_plot_scaled(leaf_ATAC_scaled_normalized, "ATAC", ATAC_colors, final_split, "ATAC")

leaf_major_extension_H3K36me3_combined_scaled_normalized_plots


leaf_final_major_extension <- (leaf_major_extension_H3K36me3_combined_scaled_normalized_plots + 
                leaf_major_extension_H3K4me1_combined_scaled_normalized_plots + 
                leaf_major_extension_H3K4me3_combined_scaled_normalized_plots + 
                leaf_major_extension_H3K56ac_combined_scaled_normalized_plots)

#leaf_final_major_extension



units_apart <- c(.8)
draw_units <- rep(units_apart, 3)

leaf_final_on_off = grid.grabExpr(draw(leaf_final_major_extension, ht_gap = unit(draw_units, "cm"), padding = unit(c(1, 1, 1, 1), "cm")))
ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure3_5/imgs/metaplot_b4_after_ims", filename="leaf_major_extension.pdf", plot=leaf_final_on_off, width = 6, height = 4, units = "in")
#ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/", filename="leaf_final_on_off_small.pdf", plot=leaf_final_on_off, width = 15.5, height = 10, units = "in")
system("open /Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure3_5/imgs/metaplot_b4_after_ims/leaf_major_extension.pdf")


###############33#########################################3###########3




leaf_original_merged_genes_leaf_H3K27me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_merged_metaplots_original_leaf_H3K27me3_all.gz")
leaf_original_merged_genes_leaf_H3K36me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_merged_metaplots_original_leaf_H3K36me3_all.gz")
leaf_original_merged_genes_leaf_H3K4me1 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_merged_metaplots_original_leaf_H3K4me1_all.gz")
leaf_original_merged_genes_leaf_H3K4me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_merged_metaplots_original_leaf_H3K4me3_all.gz")
leaf_original_merged_genes_leaf_H3K56ac <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_merged_metaplots_original_leaf_H3K56ac_all.gz")
leaf_original_merged_genes_leaf_input <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_merged_metaplots_original_leaf_input_all.gz")



leaf_merged_genes_leaf_H3K27me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_passing_after_merged_genes_leaf_H3K27me3_all.gz")
leaf_merged_genes_leaf_H3K36me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_passing_after_merged_genes_leaf_H3K36me3_all.gz")
leaf_merged_genes_leaf_H3K4me1 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_passing_after_merged_genes_leaf_H3K4me1_all.gz")
leaf_merged_genes_leaf_H3K4me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_passing_after_merged_genes_leaf_H3K4me3_all.gz")
leaf_merged_genes_leaf_H3K56ac <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_passing_after_merged_genes_leaf_H3K56ac_all.gz")
leaf_merged_genes_leaf_input <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_passing_after_merged_genes_leaf_input_all.gz")



#leaf_H2AZ <- combine_identical_marks_on_off(leaf_on_genes_leaf_H2AZ, leaf_off_genes_leaf_H2AZ)
merged_leaf_H3K27me3 <- combine_identical_marks_on_off(leaf_original_merged_genes_leaf_H3K27me3, leaf_merged_genes_leaf_H3K27me3)
merged_leaf_H3K36me3 <- combine_identical_marks_on_off(leaf_original_merged_genes_leaf_H3K36me3, leaf_merged_genes_leaf_H3K36me3)
merged_leaf_H3K4me1 <- combine_identical_marks_on_off(leaf_original_merged_genes_leaf_H3K4me1, leaf_merged_genes_leaf_H3K4me1)
merged_leaf_H3K4me3 <- combine_identical_marks_on_off(leaf_original_merged_genes_leaf_H3K4me3, leaf_merged_genes_leaf_H3K4me3)
merged_leaf_H3K56ac <- combine_identical_marks_on_off(leaf_original_merged_genes_leaf_H3K56ac, leaf_merged_genes_leaf_H3K56ac)
merged_leaf_input <- combine_identical_marks_on_off(leaf_original_merged_genes_leaf_input, leaf_merged_genes_leaf_input)
#leaf_ATAC_scaled <- combine_identical_marks_on_off(leaf_on_genes_leaf_ATAC_scaled, leaf_off_genes_leaf_ATAC_scaled)


annotation_N <- nrow(leaf_merged_genes_leaf_H3K36me3)
original_annotation_N <- nrow(leaf_original_merged_genes_leaf_H3K36me3)
values <- c("Old Annotations", "Updated Annotation")
final_split <- rep(values, times = c(original_annotation_N, annotation_N))
#leaf_sorted_gene_list <- generate_ordered_gene_list(leaf_TPM_by_gene, leaf_H3K36me3)



#leaf_merged_H2AZ_normalized <- (leaf_H2AZ)-(leaf_input)
leaf_merged_H3K27me3_normalized <- (merged_leaf_H3K27me3)-(merged_leaf_input)
leaf_merged_H3K36me3_normalized <- (merged_leaf_H3K36me3)-(merged_leaf_input)
leaf_merged_H3K4me1_normalized <- (merged_leaf_H3K4me1)-(merged_leaf_input)
leaf_merged_H3K4me3_normalized <- (merged_leaf_H3K4me3)-(merged_leaf_input)
leaf_merged_H3K56ac_normalized <- (merged_leaf_H3K56ac)-(merged_leaf_input)
#leaf_merged_ATAC_scaled_normalized <- (leaf_ATAC_scaled)-(leaf_input)



H3K36me3_colors <- "#CB2026"
H3K4me1_colors <-  "#642165"
H3K4me3_colors <- "#B89BC9"
H3K56ac_colors <- "#99CA3C"
ATAC_colors <- "#FFCD05"
H2A_colors <-  "#66c2a5"
H3K27me3_colors <- "#234bb8"



#leaf_H2AZ_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_H2AZ_normalized, "H2Az", H2A_colors, final_split, "CHIP")
leaf_merged_H3K36me3_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_merged_H3K36me3_normalized, 'H3K36me3', H3K36me3_colors, final_split, "CHIP")
leaf_merged_H3K4me1_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_merged_H3K4me1_normalized, "H3K4me1", H3K4me1_colors, final_split, "CHIP")
leaf_merged_H3K4me3_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_merged_H3K4me3_normalized, "H3K4me3", H3K4me3_colors, final_split, "CHIP")
leaf_merged_H3K56ac_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_merged_H3K56ac_normalized, "H3K56ac", H3K56ac_colors, final_split, "CHIP")
leaf_merged_H3K27me3_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_merged_H3K27me3_normalized, "H3K27me3", H3K27me3_colors, final_split, "CHIP")
#leaf_ATAC_R2C2_plot_sclaed_region_og <- generate_final_plot_scaled(leaf_ATAC_scaled_normalized, "ATAC", ATAC_colors, final_split, "ATAC")



leaf_final_merged <- (leaf_merged_H3K36me3_combined_scaled_normalized_plots + 
                                 leaf_merged_H3K4me1_combined_scaled_normalized_plots + 
                                 leaf_merged_H3K4me3_combined_scaled_normalized_plots + 
                                 leaf_merged_H3K56ac_combined_scaled_normalized_plots)

#leaf_final_merged



units_apart <- c(.8)
draw_units <- rep(units_apart, 3)

leaf_final_on_off = grid.grabExpr(draw(leaf_final_merged, ht_gap = unit(draw_units, "cm"), padding = unit(c(1, 1, 1, 1), "cm")))
ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure3_5/imgs/metaplot_b4_after_ims", filename="leaf_merged.pdf", plot=leaf_final_on_off, width = 6, height = 4, units = "in")
#ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/", filename="leaf_final_on_off_small.pdf", plot=leaf_final_on_off, width = 15.5, height = 10, units = "in")
system("open /Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure3_5/imgs/metaplot_b4_after_ims/leaf_merged.pdf")



###############33#########################################3###########3




leaf_original_minor_extension_genes_leaf_H3K27me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_minor_extension_metaplots_original_leaf_H3K27me3_all.gz")
leaf_original_minor_extension_genes_leaf_H3K36me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_minor_extension_metaplots_original_leaf_H3K36me3_all.gz")
leaf_original_minor_extension_genes_leaf_H3K4me1 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_minor_extension_metaplots_original_leaf_H3K4me1_all.gz")
leaf_original_minor_extension_genes_leaf_H3K4me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_minor_extension_metaplots_original_leaf_H3K4me3_all.gz")
leaf_original_minor_extension_genes_leaf_H3K56ac <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_minor_extension_metaplots_original_leaf_H3K56ac_all.gz")
leaf_original_minor_extension_genes_leaf_input <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_minor_extension_metaplots_original_leaf_input_all.gz")



leaf_minor_extension_genes_leaf_H3K27me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_passing_after_minor_extension_genes_leaf_H3K27me3_all.gz")
leaf_minor_extension_genes_leaf_H3K36me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_passing_after_minor_extension_genes_leaf_H3K36me3_all.gz")
leaf_minor_extension_genes_leaf_H3K4me1 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_passing_after_minor_extension_genes_leaf_H3K4me1_all.gz")
leaf_minor_extension_genes_leaf_H3K4me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_passing_after_minor_extension_genes_leaf_H3K4me3_all.gz")
leaf_minor_extension_genes_leaf_H3K56ac <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_passing_after_minor_extension_genes_leaf_H3K56ac_all.gz")
leaf_minor_extension_genes_leaf_input <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_passing_after_minor_extension_genes_leaf_input_all.gz")



#leaf_H2AZ <- combine_identical_marks_on_off(leaf_on_genes_leaf_H2AZ, leaf_off_genes_leaf_H2AZ)
minor_extension_leaf_H3K27me3 <- combine_identical_marks_on_off(leaf_original_minor_extension_genes_leaf_H3K27me3, leaf_minor_extension_genes_leaf_H3K27me3)
minor_extension_leaf_H3K36me3 <- combine_identical_marks_on_off(leaf_original_minor_extension_genes_leaf_H3K36me3, leaf_minor_extension_genes_leaf_H3K36me3)
minor_extension_leaf_H3K4me1 <- combine_identical_marks_on_off(leaf_original_minor_extension_genes_leaf_H3K4me1, leaf_minor_extension_genes_leaf_H3K4me1)
minor_extension_leaf_H3K4me3 <- combine_identical_marks_on_off(leaf_original_minor_extension_genes_leaf_H3K4me3, leaf_minor_extension_genes_leaf_H3K4me3)
minor_extension_leaf_H3K56ac <- combine_identical_marks_on_off(leaf_original_minor_extension_genes_leaf_H3K56ac, leaf_minor_extension_genes_leaf_H3K56ac)
minor_extension_leaf_input <- combine_identical_marks_on_off(leaf_original_minor_extension_genes_leaf_input, leaf_minor_extension_genes_leaf_input)
#leaf_ATAC_scaled <- combine_identical_marks_on_off(leaf_on_genes_leaf_ATAC_scaled, leaf_off_genes_leaf_ATAC_scaled)


annotation_N <- nrow(leaf_minor_extension_genes_leaf_H3K36me3)
original_annotation_N <- nrow(leaf_original_minor_extension_genes_leaf_H3K36me3)
values <- c("Old Annotations", "Updated Annotation")
final_split <- rep(values, times = c(original_annotation_N, annotation_N))
#leaf_sorted_gene_list <- generate_ordered_gene_list(leaf_TPM_by_gene, leaf_H3K36me3)



#leaf_minor_extension_H2AZ_normalized <- (leaf_H2AZ)-(leaf_input)
leaf_minor_extension_H3K27me3_normalized <- (minor_extension_leaf_H3K27me3)-(minor_extension_leaf_input)
leaf_minor_extension_H3K36me3_normalized <- (minor_extension_leaf_H3K36me3)-(minor_extension_leaf_input)
leaf_minor_extension_H3K4me1_normalized <- (minor_extension_leaf_H3K4me1)-(minor_extension_leaf_input)
leaf_minor_extension_H3K4me3_normalized <- (minor_extension_leaf_H3K4me3)-(minor_extension_leaf_input)
leaf_minor_extension_H3K56ac_normalized <- (minor_extension_leaf_H3K56ac)-(minor_extension_leaf_input)
#leaf_minor_extension_ATAC_scaled_normalized <- (leaf_ATAC_scaled)-(leaf_input)



H3K36me3_colors <- "#CB2026"
H3K4me1_colors <-  "#642165"
H3K4me3_colors <- "#B89BC9"
H3K56ac_colors <- "#99CA3C"
ATAC_colors <- "#FFCD05"
H2A_colors <-  "#66c2a5"
H3K27me3_colors <- "#234bb8"



#leaf_H2AZ_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_H2AZ_normalized, "H2Az", H2A_colors, final_split, "CHIP")
leaf_minor_extension_H3K36me3_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_minor_extension_H3K36me3_normalized, 'H3K36me3', H3K36me3_colors, final_split, "CHIP")
leaf_minor_extension_H3K4me1_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_minor_extension_H3K4me1_normalized, "H3K4me1", H3K4me1_colors, final_split, "CHIP")
leaf_minor_extension_H3K4me3_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_minor_extension_H3K4me3_normalized, "H3K4me3", H3K4me3_colors, final_split, "CHIP")
leaf_minor_extension_H3K56ac_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_minor_extension_H3K56ac_normalized, "H3K56ac", H3K56ac_colors, final_split, "CHIP")
leaf_minor_extension_H3K27me3_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_minor_extension_H3K27me3_normalized, "H3K27me3", H3K27me3_colors, final_split, "CHIP")
#leaf_ATAC_R2C2_plot_sclaed_region_og <- generate_final_plot_scaled(leaf_ATAC_scaled_normalized, "ATAC", ATAC_colors, final_split, "ATAC")
leaf_minor_extension_H3K4me1_combined_scaled_normalized_plots


leaf_final_minor_extension <- (leaf_minor_extension_H3K36me3_combined_scaled_normalized_plots + 
                        leaf_minor_extension_H3K4me1_combined_scaled_normalized_plots + 
                        leaf_minor_extension_H3K4me3_combined_scaled_normalized_plots + 
                        leaf_minor_extension_H3K56ac_combined_scaled_normalized_plots)

#leaf_final_minor_extension



units_apart <- c(.8)
draw_units <- rep(units_apart, 3)

leaf_final_on_off = grid.grabExpr(draw(leaf_final_minor_extension, ht_gap = unit(draw_units, "cm"), padding = unit(c(1, 1, 1, 1), "cm")))
ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure3_5/imgs/metaplot_b4_after_ims", filename="leaf_minor_extension.pdf", plot=leaf_final_on_off, width = 6, height = 4, units = "in")
#ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/", filename="leaf_final_on_off_small.pdf", plot=leaf_final_on_off, width = 15.5, height = 10, units = "in")
system("open /Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure3_5/imgs/metaplot_b4_after_ims/leaf_minor_extension.pdf")





#Hyper Large Genes
###############33#########################################3###########3




leaf_original_hyper_large_genes_leaf_H3K27me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_hyper_large_metaplots_original_leaf_H3K27me3_all.gz")
leaf_original_hyper_large_genes_leaf_H3K36me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_hyper_large_metaplots_original_leaf_H3K36me3_all.gz")
leaf_original_hyper_large_genes_leaf_H3K4me1 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_hyper_large_metaplots_original_leaf_H3K4me1_all.gz")
leaf_original_hyper_large_genes_leaf_H3K4me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_hyper_large_metaplots_original_leaf_H3K4me3_all.gz")
leaf_original_hyper_large_genes_leaf_H3K56ac <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_hyper_large_metaplots_original_leaf_H3K56ac_all.gz")
leaf_original_hyper_large_genes_leaf_input <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_hyper_large_metaplots_original_leaf_input_all.gz")



leaf_hyper_large_genes_leaf_H3K27me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_passing_after_hyper_large_genes_leaf_H3K27me3_all.gz")
leaf_hyper_large_genes_leaf_H3K36me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_passing_after_hyper_large_genes_leaf_H3K36me3_all.gz")
leaf_hyper_large_genes_leaf_H3K4me1 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_passing_after_hyper_large_genes_leaf_H3K4me1_all.gz")
leaf_hyper_large_genes_leaf_H3K4me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_passing_after_hyper_large_genes_leaf_H3K4me3_all.gz")
leaf_hyper_large_genes_leaf_H3K56ac <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_passing_after_hyper_large_genes_leaf_H3K56ac_all.gz")
leaf_hyper_large_genes_leaf_input <- read_gunzipped_file("07.scaled_to_region/CHIP/leaf_passing_after_hyper_large_genes_leaf_input_all.gz")



#leaf_H2AZ <- combine_identical_marks_on_off(leaf_on_genes_leaf_H2AZ, leaf_off_genes_leaf_H2AZ)
hyper_large_leaf_H3K27me3 <- combine_identical_marks_on_off(leaf_original_hyper_large_genes_leaf_H3K27me3, leaf_hyper_large_genes_leaf_H3K27me3)
hyper_large_leaf_H3K36me3 <- combine_identical_marks_on_off(leaf_original_hyper_large_genes_leaf_H3K36me3, leaf_hyper_large_genes_leaf_H3K36me3)
hyper_large_leaf_H3K4me1 <- combine_identical_marks_on_off(leaf_original_hyper_large_genes_leaf_H3K4me1, leaf_hyper_large_genes_leaf_H3K4me1)
hyper_large_leaf_H3K4me3 <- combine_identical_marks_on_off(leaf_original_hyper_large_genes_leaf_H3K4me3, leaf_hyper_large_genes_leaf_H3K4me3)
hyper_large_leaf_H3K56ac <- combine_identical_marks_on_off(leaf_original_hyper_large_genes_leaf_H3K56ac, leaf_hyper_large_genes_leaf_H3K56ac)
hyper_large_leaf_input <- combine_identical_marks_on_off(leaf_original_hyper_large_genes_leaf_input, leaf_hyper_large_genes_leaf_input)
#leaf_ATAC_scaled <- combine_identical_marks_on_off(leaf_on_genes_leaf_ATAC_scaled, leaf_off_genes_leaf_ATAC_scaled)


annotation_N <- nrow(leaf_hyper_large_genes_leaf_H3K36me3)
original_annotation_N <- nrow(leaf_original_hyper_large_genes_leaf_H3K36me3)
values <- c("Old Annotations", "Updated Annotation")
final_split <- rep(values, times = c(original_annotation_N, annotation_N))
#leaf_sorted_gene_list <- generate_ordered_gene_list(leaf_TPM_by_gene, leaf_H3K36me3)



#leaf_hyper_large_H2AZ_normalized <- (leaf_H2AZ)-(leaf_input)
leaf_hyper_large_H3K27me3_normalized <- (hyper_large_leaf_H3K27me3)-(hyper_large_leaf_input)
leaf_hyper_large_H3K36me3_normalized <- (hyper_large_leaf_H3K36me3)-(hyper_large_leaf_input)
leaf_hyper_large_H3K4me1_normalized <- (hyper_large_leaf_H3K4me1)-(hyper_large_leaf_input)
leaf_hyper_large_H3K4me3_normalized <- (hyper_large_leaf_H3K4me3)-(hyper_large_leaf_input)
leaf_hyper_large_H3K56ac_normalized <- (hyper_large_leaf_H3K56ac)-(hyper_large_leaf_input)
#leaf_hyper_large_ATAC_scaled_normalized <- (leaf_ATAC_scaled)-(leaf_input)



H3K36me3_colors <- "#CB2026"
H3K4me1_colors <-  "#642165"
H3K4me3_colors <- "#B89BC9"
H3K56ac_colors <- "#99CA3C"
ATAC_colors <- "#FFCD05"
H2A_colors <-  "#66c2a5"
H3K27me3_colors <- "#234bb8"



#leaf_H2AZ_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_H2AZ_normalized, "H2Az", H2A_colors, final_split, "CHIP")
leaf_hyper_large_H3K36me3_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_hyper_large_H3K36me3_normalized, 'H3K36me3', H3K36me3_colors, final_split, "CHIP")
leaf_hyper_large_H3K4me1_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_hyper_large_H3K4me1_normalized, "H3K4me1", H3K4me1_colors, final_split, "CHIP")
leaf_hyper_large_H3K4me3_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_hyper_large_H3K4me3_normalized, "H3K4me3", H3K4me3_colors, final_split, "CHIP")
leaf_hyper_large_H3K56ac_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_hyper_large_H3K56ac_normalized, "H3K56ac", H3K56ac_colors, final_split, "CHIP")
leaf_hyper_large_H3K27me3_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_hyper_large_H3K27me3_normalized, "H3K27me3", H3K27me3_colors, final_split, "CHIP")
#leaf_ATAC_R2C2_plot_sclaed_region_og <- generate_final_plot_scaled(leaf_ATAC_scaled_normalized, "ATAC", ATAC_colors, final_split, "ATAC")
leaf_hyper_large_H3K4me1_combined_scaled_normalized_plots


leaf_final_hyper_large <- (leaf_hyper_large_H3K36me3_combined_scaled_normalized_plots + 
                                 leaf_hyper_large_H3K4me1_combined_scaled_normalized_plots + 
                                 leaf_hyper_large_H3K4me3_combined_scaled_normalized_plots + 
                                 leaf_hyper_large_H3K56ac_combined_scaled_normalized_plots)

#leaf_final_hyper_large



units_apart <- c(.8)
draw_units <- rep(units_apart, 3)

leaf_final_on_off = grid.grabExpr(draw(leaf_final_hyper_large, ht_gap = unit(draw_units, "cm"), padding = unit(c(1, 1, 1, 1), "cm")))
ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure3_5/imgs/metaplot_b4_after_ims", filename="leaf_hyper_large.pdf", plot=leaf_final_on_off, width = 6, height = 4, units = "in")
#ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/", filename="leaf_final_on_off_small.pdf", plot=leaf_final_on_off, width = 15.5, height = 10, units = "in")
system("open /Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure3_5/imgs/metaplot_b4_after_ims/leaf_hyper_large.pdf")










#####################################333###########################3333###############33333
#Novel Section 

leaf_novel_genes_leaf_H3K36me3 <- read_gunzipped_file("07.scaled_to_region/CHIP//leaf_novel_passing.final_leaf_H3K36me3_all.gz")
leaf_novel_genes_leaf_H3K4me1 <- read_gunzipped_file("07.scaled_to_region/CHIP//leaf_novel_passing.final_leaf_H3K4me1_all.gz")
leaf_novel_genes_leaf_H3K4me3 <- read_gunzipped_file("07.scaled_to_region/CHIP//leaf_novel_passing.final_leaf_H3K4me3_all.gz")
leaf_novel_genes_leaf_H3K56ac <- read_gunzipped_file("07.scaled_to_region/CHIP//leaf_novel_passing.final_leaf_H3K56ac_all.gz")
leaf_novel_genes_leaf_input <- read_gunzipped_file("07.scaled_to_region/CHIP//leaf_novel_passing.final_leaf_input_all.gz")



leaf_novel_genes_leaf_H3K36me3_normalized <- leaf_novel_genes_leaf_H3K36me3 - leaf_novel_genes_leaf_input
leaf_novel_genes_leaf_H3K4me1_normalized <- leaf_novel_genes_leaf_H3K4me1 - leaf_novel_genes_leaf_input
leaf_novel_genes_leaf_H3K4me3_normalized <- leaf_novel_genes_leaf_H3K4me3 - leaf_novel_genes_leaf_input
leaf_novel_genes_leaf_H3K56ac_normalized <- leaf_novel_genes_leaf_H3K56ac -leaf_novel_genes_leaf_input 


leaf_novel_genes_leaf_H3K36me3_normalized_plots <- generate_final_plot_scaled(leaf_novel_genes_leaf_H3K36me3_normalized, 'H3K36me3', H3K36me3_colors, NULL, "CHIP")
leaf_novel_genes_leaf_H3K36me3_normalized_plots

leaf_novel_genes_leaf_H3K4me1_normalized_plots <- generate_final_plot_scaled(leaf_novel_genes_leaf_H3K4me1_normalized, "H3K4me1", H3K4me1_colors, NULL, "CHIP")
leaf_novel_genes_leaf_H3K4me3_normalized_plots <- generate_final_plot_scaled(leaf_novel_genes_leaf_H3K4me3_normalized, "H3K4me3", H3K4me3_colors, NULL, "CHIP")
leaf_novel_genes_leaf_H3K56ac_normalized_plots <- generate_final_plot_scaled(leaf_novel_genes_leaf_H3K56ac_normalized, "H3K56ac", H3K56ac_colors, NULL, "CHIP")



leaf_novel_normalized_plots <- (leaf_novel_genes_leaf_H3K36me3_normalized_plots + 
                                  leaf_novel_genes_leaf_H3K4me1_normalized_plots + 
                                  leaf_novel_genes_leaf_H3K4me3_normalized_plots + 
                                  leaf_novel_genes_leaf_H3K56ac_normalized_plots )

draw(leaf_final_minor_extension, ht_gap = unit(draw_units, "cm"), padding = unit(c(1, 1, 1, 1), "cm"))




units_apart <- c(.8)
draw_units <- rep(units_apart, 3)

leaf_final_on_off = grid.grabExpr(draw(leaf_novel_normalized_plots, ht_gap = unit(draw_units, "cm"), padding = unit(c(1, 1, 1, 1), "cm")))
ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure3_5/imgs/metaplot_b4_after_ims", filename="leaf_novel_regions.pdf", plot=leaf_final_on_off, width = 6, height = 4, units = "in")
#ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/", filename="leaf_final_on_off_small.pdf", plot=leaf_final_on_off, width = 15.5, height = 10, units = "in")
system("open /Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure3_5/imgs/metaplot_b4_after_ims/leaf_novel_regions.pdf")




# Pull Top Annotations For Graphic  ---------------------------------------




leaf_novel_annotaion_draw <- draw(leaf_novel_normalized_plots, ht_gap = unit(draw_units, "cm"), padding = unit(c(1, 1, 1, 1), "cm"))
leaf_final_merged_draw = draw(leaf_final_merged, ht_gap = unit(draw_units, "cm"), padding = unit(c(1, 1, 1, 1), "cm"))
leaf_final_major_extension_draw = draw(leaf_final_major_extension, ht_gap = unit(draw_units, "cm"), padding = unit(c(1, 1, 1, 1), "cm"))
leaf_final_minor_extension_draw = draw(leaf_final_minor_extension, ht_gap = unit(draw_units, "cm"), padding = unit(c(1, 1, 1, 1), "cm"))
leaf_final_hyper_large_draw = draw(leaf_final_hyper_large, ht_gap = unit(draw_units, "cm"), padding = unit(c(1, 1, 1, 1), "cm"))


generate_top_annotation <- function(drawn_annotation_type, save_file_name) {
  
  
  add_anno_enriched = function(ht_list, name, ri, ci) {
    pushViewport(viewport(layout.pos.row = ri, layout.pos.col = ci))
    extract_anno_enriched(ht_list, name, newpage = FALSE)
    upViewport()
  }
  
  
  
  pdf(save_file_name, width = 10, height = 3)
  
  
  pushViewport(viewport(layout = grid.layout(nr = 1, nc = 4)))
  add_anno_enriched(drawn_annotation_type, 1,     1, 1)
  add_anno_enriched(drawn_annotation_type, 2,    1, 2)
  add_anno_enriched(drawn_annotation_type, 3,    1, 3)
  add_anno_enriched(drawn_annotation_type, 4, 1, 4)

  dev.off()

}



generate_top_annotation(leaf_final_merged_draw, "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure3_5/imgs/metaplot_b4_after_ims/top_merged_metaplot_before_after.pdf")
generate_top_annotation(leaf_novel_annotaion_draw,"/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure3_5/imgs/metaplot_b4_after_ims/top_novel_annotation_metaplot_before_after.pdf")
generate_top_annotation(leaf_final_major_extension_draw,"/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure3_5/imgs/metaplot_b4_after_ims/top_major_extension_metaplot_before_after.pdf")
generate_top_annotation(leaf_final_minor_extension_draw,"/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure3_5/imgs/metaplot_b4_after_ims/top_minor_extension_metaplot_before_after.pdf")
generate_top_annotation(leaf_final_hyper_large_draw,"/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure3_5/imgs/metaplot_b4_after_ims/top_hyper_large_metaplot_before_after.pdf")





