library(tidyverse)
library(EnrichedHeatmap)
library(cowplot)
library("grid")
library("ggplotify")
library("here")
library(circlize)


H3K36me3_colors <- "#CB2026"
H3K4me1_colors <-  "#642165"
H3K4me3_colors <- "#B89BC9"
H3K56ac_colors <- "#99CA3C"
ATAC_colors <- "#FFCD05"
H2A_colors <-  "#66c2a5"
H3K27me3_colors <- "#234bb8"



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
  quick <- TSS_matrix[leaf_H3K36me3_normalized > 0]

  
  #Scale colors across matricies the same
  #col_fun = colorRamp2(quantile(TSS_matrix, c(0, .9)), c("white", color_hex))
  #col_fun = circlize::colorRamp2(c(common_min, common_max), c("white", color_hex))
  
  col_fun = colorRamp2(quantile(quick, c(.1, .8)), c("white", color_hex))
  
  
  #Get the same max value for each 
  #y_max_val <- round(max(colMeans(TES_matrix, na.rm = TRUE), colMeans(TSS_matrix, na.rm = TRUE))) + 2 
  
  axis_name = c("-1kb","TSS","TTS","1kb")
  final_graph <-  EnrichedHeatmap(TSS_matrix,  column_title =  mark_name, use_raster = TRUE, show_heatmap_legend = FALSE,
                                  row_split = row_split_array, cluster_rows = FALSE, 
                                  column_title_gp = gpar(fontsize = 25, fontface = "bold"), 
                                  col = col_fun, axis_name_rot = 90, row_title_rot = 90, axis_name_gp = gpar(fontsize = 18), axis_name = axis_name,
                                  pos_line_gp = gpar(lty = 3),
                                  top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = c(H3K27me3_colors, "#888A8D"))))) 
  
  return(final_graph)
  
  
}
generate_complex_heatmaps_ATAC_scaled_region <- function(final_matrix_list, mark_name, row_split_array, color_hex) {
  
  #If row title is NOT missing - pass on it and don't add it 
  
  TSS_matrix <- final_matrix_list
  TSS_matrix[TSS_matrix < 0] = 0
  #col_fun <- circlize::colorRamp2(c(0, common_max), c("white", color_hex))
  #col_fun = colorRamp2(quantile(TSS_matrix, c(0, .9)), c("white", color_hex))
  col_fun = colorRamp2(quantile(TSS_matrix, c(.1, .85)), c("white", color_hex))
  
  axis_name = c("-1kb","TSS","TTS","1kb")
  #This is the only difference between ATAC and ChIP-seq metaplots. Basically we can't scale the same way we did 
  #previously. This allows us to scale the top part of the metaplot appropriatly.
  
  final_graph <-  EnrichedHeatmap(TSS_matrix,  column_title =  mark_name, use_raster = TRUE, show_heatmap_legend = FALSE,
                                  row_split = row_split_array, cluster_rows = FALSE, 
                                  column_title_gp = gpar(fontsize = 25, fontface = "bold"), axis_name = axis_name, 
                                  col = col_fun, axis_name_rot = 90, row_title_rot = 90, axis_name_gp = gpar(fontsize = 18),
                                  pos_line_gp = gpar(lty = 3),
                                  top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = c("#888A8D",H3K27me3_colors))))) 
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
  read_TPM_gene <- read_delim(TPM_file, '\t', col_names = col_names_TPM_file, col_types = "cd") %>% 
    mutate(updated_name = str_c("gene", gene_ID , sep =':')) %>% 
    select(-gene_ID) %>% 
    dplyr::rename(gene_ID = updated_name)
  
  return(read_TPM_gene)
}
alternative_TPM_file <- function(TPM_file) {
  col_names_TPM_file <- c("gene_ID", "Chr", "Start", "End", "Length", "Reads", "TPM")
  read_TPM_gene <- read_delim(TPM_file, '\t', col_names = col_names_TPM_file, col_types = "cdddddd") %>% 
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



#leaf_on_genes_leaf_H2AZ <- read_gunzipped_file("07.scaled_to_region/CHIP/all_genes_leaf_leaf_H2AZ_all.gz")
#leaf_on_genes_leaf_H3K27me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/all_genes_leaf_leaf_H3K27me3_all.gz")
#leaf_on_genes_leaf_H3K36me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/all_genes_leaf_leaf_H3K36me3_all.gz")
#leaf_on_genes_leaf_H3K4me1 <- read_gunzipped_file("07.scaled_to_region/CHIP/all_genes_leaf_leaf_H3K4me1_all.gz")
#leaf_on_genes_leaf_H3K4me3 <- read_gunzipped_file("07.scaled_to_region/CHIP/all_genes_leaf_leaf_H3K4me3_all.gz")
#leaf_on_genes_leaf_H3K56ac <- read_gunzipped_file("07.scaled_to_region/CHIP/all_genes_leaf_leaf_H3K56ac_all.gz")
#leaf_on_genes_leaf_input <- read_gunzipped_file("07.scaled_to_region/CHIP/all_genes_leaf_leaf_input_all.gz")
#leaf_on_genes_leaf_ATAC_scaled <- read_gunzipped_file("07.scaled_to_region/ATAC/all_genes_leaf_leaf_ATAC_scaled_all.gz")



###################3



#Test
#############
#leaf_TPM_by_gene_v_2 <- alternative_TPM_file("/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/00.data/TPM_vals/RNA_B73_leaf_Aligned.sortedByCoord.out_genes.out")
##leaf_TPM_by_gene_v_2 = leaf_TPM_by_gene_v_2[-1,]
#leaf_TPM_by_gene_v_2$TPM <- as.numeric(leaf_TPM_by_gene_v_2$TPM)
#
#
###
#test_test <- head(leaf_on_genes_leaf_H3K36me3, n = 5000)
#
##leaf_sorted_on_genes <- generate_ordered_gene_list(leaf_TPM_by_gene, test_test)
#  
#gather_gene_row_names <- row.names(test_test)
#
#summed_vals <- (rowSums(test_test))
#
#leaf_TPM_by_gene_sorted <- leaf_TPM_by_gene_v_2 %>% 
#  filter(gene_ID %in% gather_gene_row_names) %>% 
#  mutate(log2_tpm = log2(TPM +1)) %>% 
#  dplyr::arrange(desc(log2_tpm))
#
#
#sorted_expression <- leaf_TPM_by_gene_sorted$gene_ID
#
#
#
#common_max <-  max(colMeans(leaf_on_genes_leaf_H3K36me3, na.rm = TRUE))
#common_min <- min(leaf_on_genes_leaf_H3K36me3)
#
#
##col_fun = colorRamp2(quantile(TSS_matrix, c(0, .9)), c("white", color_hex))
#col_fun = circlize::colorRamp2(c(common_min, common_max), c("white", H3K36me3_colors))
#
#col_fun = colorRamp2(quantile(leaf_on_genes_leaf_H3K36me3, c(0, 1)), c("white", H3K36me3_colors))
#
#
#
#quantile(leaf_on_genes_leaf_H3K36me3, c(0, 1))
#
#
#test_normalized <- generate_tss_heatmap_scaled_region(test_test, "TEST")
#
#
#
#quantile(leaf_on_genes_leaf_H3K36me3, c(0, .8))
#
#test_graph <- EnrichedHeatmap(test_normalized,  column_title =  "K4me1", use_raster = TRUE, show_heatmap_legend = FALSE,
#                                cluster_rows = FALSE, col = col_fun, 
#                                column_title_gp = gpar(fontsize = 25, fontface = "bold"), 
#                                axis_name_rot = 90, row_title_rot = 90, axis_name_gp = gpar(fontsize = 18), 
#                                axis_name = c("-1000bp","TSS","TTS","1000bp"),
#                                pos_line_gp = gpar(lty = 3),
#                                top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4)))) 
#
#catch <- draw(test_graph)
##row_order(catch)
#
#catch_2 <- draw(test_graph, row_order = sorted_expression)
##row_order(catch_2)
#
#############




leaf_H2AZ <- combine_identical_marks_on_off(leaf_on_genes_leaf_H2AZ, leaf_off_genes_leaf_H2AZ)
leaf_H3K27me3 <- combine_identical_marks_on_off(leaf_on_genes_leaf_H3K27me3, leaf_off_genes_leaf_H3K27me3)
leaf_H3K36me3 <- combine_identical_marks_on_off(leaf_on_genes_leaf_H3K36me3, leaf_off_genes_leaf_H3K36me3)
leaf_H3K4me1 <- combine_identical_marks_on_off(leaf_on_genes_leaf_H3K4me1, leaf_off_genes_leaf_H3K4me1)
leaf_H3K4me3 <- combine_identical_marks_on_off(leaf_on_genes_leaf_H3K4me3, leaf_off_genes_leaf_H3K4me3)
leaf_H3K56ac <- combine_identical_marks_on_off(leaf_on_genes_leaf_H3K56ac, leaf_off_genes_leaf_H3K56ac)
leaf_input <- combine_identical_marks_on_off(leaf_on_genes_leaf_input, leaf_off_genes_leaf_input)
leaf_ATAC_scaled <- combine_identical_marks_on_off(leaf_on_genes_leaf_ATAC_scaled, leaf_off_genes_leaf_ATAC_scaled)



number_off <- nrow(leaf_off_genes_leaf_H3K36me3)
number_on <- nrow(leaf_on_genes_leaf_H3K36me3)
values <- c( "Expressed Genes", "Non-Expressed Genes")
final_split <- rep(values, times = c(number_on, number_off))
leaf_sorted_gene_list <- generate_ordered_gene_list(leaf_TPM_by_gene, leaf_H3K36me3)



leaf_H2AZ_normalized <- (leaf_H2AZ)-(leaf_input)
leaf_H3K27me3_normalized <- (leaf_H3K27me3)-(leaf_input)
leaf_H3K36me3_normalized <- (leaf_H3K36me3)-(leaf_input)
leaf_H3K4me1_normalized <- (leaf_H3K4me1)-(leaf_input)
leaf_H3K4me3_normalized <- (leaf_H3K4me3)-(leaf_input)
leaf_H3K56ac_normalized <- (leaf_H3K56ac)-(leaf_input)
leaf_ATAC_scaled_normalized <- (leaf_ATAC_scaled)-(leaf_input)


  



leaf_H2AZ_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_H2AZ_normalized, "H2Az", H2A_colors, final_split, "CHIP")
leaf_H3K36me3_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_H3K36me3_normalized, 'H3K36me3', H3K36me3_colors, final_split, "CHIP")
leaf_H3K4me1_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_H3K4me1_normalized, "H3K4me1", H3K4me1_colors, final_split, "CHIP")
leaf_H3K4me3_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_H3K4me3_normalized, "H3K4me3", H3K4me3_colors, final_split, "CHIP")
leaf_H3K56ac_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_H3K56ac_normalized, "H3K56ac", H3K56ac_colors, final_split, "CHIP")
leaf_H3K27me3_combined_scaled_normalized_plots <- generate_final_plot_scaled(leaf_H3K27me3_normalized, "H3K27me3", H3K27me3_colors, final_split, "CHIP")
leaf_ATAC_R2C2_plot_sclaed_region_og <- generate_final_plot_scaled(leaf_ATAC_scaled_normalized, "ATAC", ATAC_colors, final_split, "ATAC")



#draw(leaf_H3K36me3_combined_scaled_normalized_plots, row_order = leaf_sorted_gene_list)
#draw(leaf_H3K36me3_combined_scaled_normalized_plots, ht_gap = unit(draw_units, "cm"), row_order = combined, padding = unit(c(1, 1, 1, 1), "cm"))

#leaf_ATAC_R2C2_plot_sclaed_region_og
#leaf_partition = paste0("cluster", kmeans(leaf_combined_all_matricies, centers = val)$cluster)
#leaf_partition_hm <- Heatmap(final_split, col = structure(2:4, names = c("On", "Off")), show_row_names = FALSE, width = unit(3, "mm"))
#
#
#
#leaf_partition = paste0("cluster", kmeans(leaf_combined_all_matricies, centers = val)$cluster)
#leaf_partition_hm <- Heatmap(leaf_partition, col = structure(1:val, names = paste0("cluster", 1:val)), name = "partition",
#                             show_row_names = FALSE, width = unit(3, "mm"))




leaf_final <- (leaf_H3K36me3_combined_scaled_normalized_plots +
  leaf_H3K4me1_combined_scaled_normalized_plots +
  leaf_H3K4me3_combined_scaled_normalized_plots +
  leaf_H3K56ac_combined_scaled_normalized_plots +
  leaf_H3K27me3_combined_scaled_normalized_plots +
  leaf_H2AZ_combined_scaled_normalized_plots)# + 
  #leaf_ATAC_R2C2_plot_sclaed_region_og)


units_apart <- c(.8)
draw_units <- rep(units_apart, 6)

leaf_final_on_off = grid.grabExpr(draw(leaf_final, ht_gap = unit(draw_units, "cm"), row_order = leaf_sorted_gene_list), padding = unit(c(1, 1, 1, 1), "cm"))
leaf_final_on_off

ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/", filename="leaf_final_on_off.pdf", plot=leaf_final_on_off, width = 7.5, height = 4, units = "in")
#ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/", filename="leaf_final_on_off_small.pdf", plot=leaf_final_on_off, width = 15.5, height = 10, units = "in")
system("open /Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/leaf_final_on_off.pdf")




# root --------------------------------------------------------------------


root_TPM_by_gene <- alternative_TPM_file("/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/00.data/TPM_vals/RNA_B73_root_TPM.surviving.txt")
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




root_H2AZ <- combine_identical_marks_on_off(root_on_genes_root_H2AZ, root_off_genes_root_H2AZ)
root_H3K27me3 <- combine_identical_marks_on_off(root_on_genes_root_H3K27me3, root_off_genes_root_H3K27me3)
root_H3K36me3 <- combine_identical_marks_on_off(root_on_genes_root_H3K36me3, root_off_genes_root_H3K36me3)
root_H3K4me1 <- combine_identical_marks_on_off(root_on_genes_root_H3K4me1, root_off_genes_root_H3K4me1)
root_H3K4me3 <- combine_identical_marks_on_off(root_on_genes_root_H3K4me3, root_off_genes_root_H3K4me3)
root_H3K56ac <- combine_identical_marks_on_off(root_on_genes_root_H3K56ac, root_off_genes_root_H3K56ac)
root_input <- combine_identical_marks_on_off(root_on_genes_root_input, root_off_genes_root_input)
root_ATAC_scaled <- combine_identical_marks_on_off(root_on_genes_root_ATAC_scaled, root_off_genes_root_ATAC_scaled)




number_off <- nrow(root_off_genes_root_H3K36me3)
number_on <- nrow(root_on_genes_root_H3K36me3)
values <- c( "Expressed Genes", "Non-Expressed Genes")
final_split <- rep(values, times = c(number_on, number_off))
root_sorted_gene_list <- generate_ordered_gene_list(root_TPM_by_gene, root_H3K36me3)



root_H2AZ_normalized <- (root_H2AZ + 1)-(root_input +1)
root_H3K27me3_normalized <- (root_H3K27me3 + 1)-(root_input +1)
root_H3K36me3_normalized <- (root_H3K36me3 + 1)-(root_input +1)
root_H3K4me1_normalized <- (root_H3K4me1 + 1)-(root_input +1)
root_H3K4me3_normalized <- (root_H3K4me3 + 1)-(root_input +1)
root_H3K56ac_normalized <- (root_H3K56ac + 1)-(root_input +1)
root_ATAC_scaled_normalized <- (root_ATAC_scaled + 1)-(root_input +1)

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
                 root_H2AZ_combined_scaled_normalized_plots)# + 
                 #root_ATAC_R2C2_plot_sclaed_region_og)




root_final_on_off = grid.grabExpr(draw(root_final, ht_gap = unit(draw_units, "cm"), row_order = root_sorted_gene_list))
ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/", filename="root_final_on_off.pdf", plot=root_final_on_off, width = 7.5, height = 4, units = "in")
ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/", filename="root_final_on_off_small.pdf", plot=root_final_on_off, width = 15.5, height = 10, units = "in")
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




ear_TPM_by_gene <- alternative_TPM_file("/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/00.data/TPM_vals/RNA_B73_ear_TPM.surviving.txt")



ear_H2AZ <- combine_identical_marks_on_off(ear_on_genes_ear_H2AZ, ear_off_genes_ear_H2AZ)
ear_H3K27me3 <- combine_identical_marks_on_off(ear_on_genes_ear_H3K27me3, ear_off_genes_ear_H3K27me3)
ear_H3K36me3 <- combine_identical_marks_on_off(ear_on_genes_ear_H3K36me3, ear_off_genes_ear_H3K36me3)
ear_H3K4me1 <- combine_identical_marks_on_off(ear_on_genes_ear_H3K4me1, ear_off_genes_ear_H3K4me1)
ear_H3K4me3 <- combine_identical_marks_on_off(ear_on_genes_ear_H3K4me3, ear_off_genes_ear_H3K4me3)
ear_H3K56ac <- combine_identical_marks_on_off(ear_on_genes_ear_H3K56ac, ear_off_genes_ear_H3K56ac)
ear_input <- combine_identical_marks_on_off(ear_on_genes_ear_input, ear_off_genes_ear_input)
ear_ATAC_scaled <- combine_identical_marks_on_off(ear_on_genes_ear_ATAC_scaled, ear_off_genes_ear_ATAC_scaled)



number_off <- nrow(ear_off_genes_ear_H3K36me3)
number_on <- nrow(ear_on_genes_ear_H3K36me3)
values <- c( "Expressed Genes", "Non-Expressed Genes")
final_split <- rep(values, times = c(number_on, number_off))
ear_sorted_gene_list <- generate_ordered_gene_list(ear_TPM_by_gene, ear_H3K36me3)





ear_H2AZ_normalized <- (ear_H2AZ + 1)-(ear_input +1)
ear_H3K27me3_normalized <- (ear_H3K27me3 + 1)-(ear_input +1)
ear_H3K36me3_normalized <- (ear_H3K36me3 + 1)-(ear_input +1)
ear_H3K4me1_normalized <- (ear_H3K4me1 + 1)-(ear_input +1)
ear_H3K4me3_normalized <- (ear_H3K4me3 + 1)-(ear_input +1)
ear_H3K56ac_normalized <- (ear_H3K56ac + 1)-(ear_input +1)
ear_ATAC_scaled_normalized <- (ear_ATAC_scaled + 1)-(ear_input +1)

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
                 ear_H2AZ_combined_scaled_normalized_plots)# + 
                 #ear_ATAC_R2C2_plot_sclaed_region_og)



ear_final_on_off = grid.grabExpr(draw(ear_final, ht_gap = unit(draw_units, "cm"), row_order = ear_sorted_gene_list))
ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/", filename="ear_final_on_off.pdf", plot=ear_final_on_off, width = 7.5, height = 4, units = "in")
#ggsave(path = "/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/", filename="ear_final_on_off_small.pdf", plot=ear_final_on_off, width = 15.5, height = 10, units = "in")
system("open /Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2/imgs/metaplots/ear_final_on_off.pdf")



