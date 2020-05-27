library(VennDiagram)
library("extrafont")
library(patchwork)
library(dplyr)
library(cowplot)
library(eulerr)
library(ggplot2)
library(devtools)

setwd("/Users/feilab/Projects/03.ncRNA_project/03.Figures/Figure2")

read_in_ext_init_vals <- function(file_name, mark_type, tissue_type){
  
  
  hist_mod_file <- read_delim(file_name,'\t', col_names = TRUE) %>% 
    mutate(mark_name = mark_type) %>% 
    mutate(tis = tissue_type)
  
  return(hist_mod_file)
  
}
generate_ven_diagram <- function(ext_init_val){
  
  H3K36me3_colors <- "#CB2026"
  H3K4me1_colors <-  "#642165"
  H3K4me3_colors <- "#B89BC9"
  H3K56ac_colors <- "#99CA3C"

  if (ext_init_val$mark_name == "Extension") {
    
    #Manuall need to check this, otherwise colors are flipped
    K36_val <- ext_init_val$H3K36me3
    K4_val <- ext_init_val$H3K4me1
    
    if (K36_val > K4_val) {
      print("K36 GREATER, FLIPPING")
      p_venn <- draw.pairwise.venn(
        area1 = ext_init_val$H3K36me3,
        area2 = ext_init_val$H3K4me1,
        cross.area = ext_init_val$both,
        ext.line.lty = "solid",
        category = c("H3K36me3","H3K4me1"),
        fill = c(H3K36me3_colors, H3K4me1_colors),
        alpha = c(0.7, 0.7),
        scale = TRUE,
        ind = TRUE,
        cex = 2, #Size of the values in the diagram
        cat.cex = 3, #Size of the Titles
        
        #WHEN THERE IS MORE K36 MAKRS YOU 
        #HAVE TO FLIP FOR CONSISTEN NAMEING
        #THE BELOW TWO LINES FLIP THE TEXT AND THE ORDER
        
        
        cat.pos = c(30,-30),
        inverted = TRUE,
        cat.dist = c(0.08, 0.08),
        fontfamily = rep("sans", 3),
        cat.fontfamily = rep("sans", 2),
        #cat.just = list(c(0, 0), c(.5, .5)),
        cat.prompts = TRUE,
        margin = 0.05,
        ext.text = TRUE, #If small as area, place labe outside
        ext.percent = .4, #If small as area, place labe outside
        ext.dist = .03,
        ext.length = 0.75,
        ext.line.lwd = 2,
        ext.pos = c(45,-45))
      
    }
    
    else if (K36_val < K4_val) {

      p_venn <- draw.pairwise.venn(
        area1 = ext_init_val$H3K36me3,
        area2 = ext_init_val$H3K4me1,
        cross.area = ext_init_val$both,
        ext.line.lty = "solid",
        category = c("H3K36me3","H3K4me1"),
        fill = c(H3K36me3_colors, H3K4me1_colors),
        alpha = c(0.7, 0.7),
        scale = TRUE,
        ind = TRUE,
        cex = 2, #Size of the values in the diagram
        cat.cex = 3, #Size of the Titles
        cat.pos = c(-30,30),
        cat.dist = c(0.08, 0.08), 
        fontfamily = rep("sans", 3),
        cat.fontfamily = rep("sans", 2),
        #cat.just = list(c(0, 0), c(.5, .5)),
        cat.prompts = TRUE,
        margin = 0.05,
        ext.text = TRUE, #If small as area, place labe outside
        ext.percent = .40, #If small as area, place labe outside
        ext.dist = .03,
        ext.length = 0.75,
        ext.line.lwd = 2,
        ext.pos = c(45,-45)) 
    }
    
    return(p_venn)
    
    
  } else if (ext_init_val$mark_name == "Initiation"){
    
    p_venn <- draw.pairwise.venn(
      area1 = ext_init_val$H3K4me3,
      area2 = ext_init_val$H3K56ac,
      cross.area = ext_init_val$both,
      ext.line.lty = "solid",
      category = c("H3K4me3","H3K56ac"),
      fill = c(H3K4me3_colors, H3K56ac_colors),
      alpha = c(0.7, 0.7),
      scale = TRUE,
      ind = TRUE,
      cex = 2, #Size of the values in the diagram
      cat.cex = 3, #Size of the Titles
      cat.pos = c(-30,30),
      cat.dist = c(0.08, 0.08), 
      fontfamily = rep("sans", 3),
      cat.fontfamily = rep("sans", 2),
      
      #cat.just = list(c(0, 0), c(.5, .5)),
      cat.prompts = TRUE,
      margin = 0.05,
      ext.text = TRUE, #If small as area, place labe outside
      ext.percent = .40, #If small as area, place labe outside
      ext.dist = .03,
      ext.length = 0.75,
      ext.line.lwd = 2,
      ext.pos = c(45,-45))
  
      
    
    
    return(p_venn)
  }

  
}


leaf_extension <- "01.ven_diagram_data/leaf_extension.txt"
leaf_initiation <- "01.ven_diagram_data/leaf_initiation.txt"
ear_extension <-"01.ven_diagram_data/ear_extension.txt"
ear_initiation <- "01.ven_diagram_data/ear_initiation.txt"
root_extension <- "01.ven_diagram_data/root_extension.txt"
root_initiation <-"01.ven_diagram_data/root_initiation.txt"


leaf_ext <- read_in_ext_init_vals(leaf_extension, "Extension", "leaf")
leaf_init <- read_in_ext_init_vals(leaf_initiation, "Initiation", "leaf")
ear_ext <- read_in_ext_init_vals(ear_extension, "Extension", "ear")
ear_init <- read_in_ext_init_vals(ear_initiation, "Initiation","ear")
root_ext <- read_in_ext_init_vals(root_extension, "Extension","root")
root_init <- read_in_ext_init_vals(root_initiation, "Initiation","root")

leaf_ext_ven <- generate_ven_diagram(leaf_ext)
leaf_init_ven <- generate_ven_diagram(leaf_init)
ear_ext_ven <- generate_ven_diagram(ear_ext)
ear_init_ven <- generate_ven_diagram(ear_init)
root_ext_ven <- generate_ven_diagram(root_ext)
root_init_ven <- generate_ven_diagram(root_init)




leaf_ext_ven_plot <- ggdraw(leaf_ext_ven) +theme(plot.background = element_rect(fill = NA),plot.margin = margin(12, 12, 12, 12))
leaf_init_ven_plot <- ggdraw(leaf_init_ven) +theme(plot.background = element_rect(fill = NA),plot.margin = margin(12, 12, 12, 12))
ear_ext_ven_plot <- ggdraw(ear_ext_ven) +theme(plot.background = element_rect(fill = NA),plot.margin = margin(12, 12, 12, 12))
ear_init_ven_plot <- ggdraw(ear_init_ven) +theme(plot.background = element_rect(fill = NA),plot.margin = margin(12, 12, 12, 12))
root_ext_ven_plot <- ggdraw(root_ext_ven) +theme(plot.background = element_rect(fill = NA),plot.margin = margin(12, 12, 12, 12))
root_init_ven_plot <- ggdraw(root_init_ven) +theme(plot.background = element_rect(fill = NA),plot.margin = margin(12, 12, 12, 12))

ggsave(plot= leaf_ext_ven_plot, filename = "imgs/leaf_ext_ven_plot.ven.svg", width = 7, height =6 )
ggsave(plot= leaf_init_ven_plot, filename = "imgs/leaf_init_ven_plot.ven.svg", width = 7, height =6 )
ggsave(plot= ear_ext_ven_plot, filename = "imgs/ear_ext_ven_plot.ven.svg", width = 7, height =6 )
ggsave(plot= ear_init_ven_plot, filename = "imgs/ear_init_ven_plot.ven.svg", width = 7, height =6 )
ggsave(plot= root_ext_ven_plot, filename = "imgs/root_ext_ven_plot.ven.svg", width = 7, height =6 )
ggsave(plot= root_init_ven_plot, filename = "imgs/root_init_ven_plot.ven.svg", width = 7, height =6 )




leaf_combined <- leaf_init_ven_plot/leaf_ext_ven_plot
ear_combined <- ear_init_ven_plot/ear_ext_ven_plot
root_combined <- root_init_ven_plot/root_ext_ven_plot

ggsave(plot= leaf_combined, filename = "imgs/leaf_combined.ven.svg", width = 6 , height =12 )
ggsave(plot= ear_combined, filename = "imgs/ear_combined.ven.svg", width = 6 , height =12 )
ggsave(plot= root_combined, filename = "imgs/root_combined.ven.svg", width = 6 , height =12 )




leaf_combined_side_by_side <- leaf_init_ven_plot|leaf_ext_ven_plot
ear_combined_side_by_side <- ear_init_ven_plot|ear_ext_ven_plot
root_combined_side_by_side <- root_init_ven_plot|root_ext_ven_plot

ggsave(plot= leaf_combined_side_by_side, filename = "imgs/leaf_combined_side_by_side.ven.svg", width = 12 , height =6 )
ggsave(plot= ear_combined_side_by_side, filename = "imgs/ear_combined_side_by_side.ven.svg", width = 12 , height =6 )
ggsave(plot= root_combined_side_by_side, filename = "imgs/root_combined_side_by_side.ven.svg", width = 12 , height =6 )



# old test ----------------------------------------------------------------


#p_venn_leaf <- draw.pairwise.venn(
#  area1 = leaf_ext$H3K36me3,
#  area2 = leaf_ext$H3K4me1,
#  cross.area = leaf_ext$both,
#  category = c("H3K36me3","H3K4me1"),
#  fill = c(H3K36me3_colors, H3K4me1_colors),
#  alpha = c(0.7, 0.7),
#  scale = TRUE,
#  ind = TRUE,
#  cex = 2, #Size of the values in the diagram
#  cat.cex = 3, #Size of the Titles
#  cat.pos = c(30,-30),
#  cat.dist = c(0.08, 0.08), 
#  #cat.just = list(c(0, 0), c(.5, .5)),
#  cat.prompts = TRUE,
#  margin = 0.05,
#  ext.text = TRUE, #If small as area, place labe outside
#  ext.percent = .15, #If small as area, place labe outside
#  ext.dist = -0.08,
#  ext.length = 0.75,
#  ext.line.lwd = 2,
#  ext.line.lty = "dashed",
#  ext.pos = 45
#
#)
#
#grid.newpage()
#p_venn_leaf
#
#ear_ext
#leaf_ext
#
#grid.newpage()
#p_venn_ear <- draw.pairwise.venn(
#  area1 = ear_ext$H3K36me3,
#  area2 = ear_ext$H3K4me1,
#  cross.area = ear_ext$both,
#  category = c("H3K36me3","H3K4me1"),
#  fill = c(H3K36me3_colors, H3K4me1_colors),
#  alpha = c(0.7, 0.7),
#  scale = TRUE,
#  ind = TRUE,
#  cex = 2, #Size of the values in the diagram
#  cat.cex = 3, #Size of the Titles
#  cat.pos = c(-30,30),
#  cat.dist = c(0.08, 0.08), 
#  #cat.just = list(c(0, 0), c(.5, .5)),
#  cat.prompts = TRUE,
#  margin = 0.05,
#  #rotation = 1,
#  inverted = TRUE,
#  ext.text = TRUE, #If small as area, place labe outside
#  ext.percent = .15, #If small as area, place labe outside
#  ext.dist = -0.08,
#  ext.length = 0.75,
#  ext.line.lwd = 2,
#  ext.line.lty = "dashed",
#  ext.pos = 45
#  
#  







