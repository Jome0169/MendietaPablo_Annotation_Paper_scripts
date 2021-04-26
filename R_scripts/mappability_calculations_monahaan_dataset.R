library(tidyverse)

setwd("/Users/feilab/Projects/03.ncRNA_project/02.Analysis/Monnahan_et_al/Mendieta_Intersection")


read_file_function <- function(file_name, group_name){
  
  read_file <- read_delim(file_name, '\t', col_names = c("chrom", "start", "stop", "gene_name", "score", "strand", "other", "other_1", "gene_name_1", "other_2", "map_score"))
  
  read_file_add_group <- read_file  
  
  
  return(read_file_add_group)
  
}
read_paired_named <- function(file_name) {
  
  
  read_file <- read_delim(file_name, '\t', col_names = c("gene_1", "gene_2")) %>% 
    mutate(id = row_number())
  return(read_file)
  
  
}


gene_mappable_bed <- read_file_function("Monnahan_et_al_validated_B73_genes.merged.Not_Mendieta_intersect.mappability_values.bed")




paired_gene_names <- read_paired_named("Monnahan_et_al_validated_B73_genes.merged.Not_Mendieta_intersect.paired_gene_names_only.txt")
paired_gene_names


paired_gene_names_final <- paired_gene_names %>% 
  pivot_longer(c("gene_1", "gene_2"), names_to = "gene_class", values_to = "gene_ID", names_repair="minimal") %>% 
  mutate(gene_name = str_c("gene", gene_ID, sep = ":"))
  
  

gene_mappable_only_name_score <- gene_mappable_bed %>% 
  select(gene_name, map_score)


final_joined <- left_join(paired_gene_names_final, gene_mappable_only_name_score, by = "gene_name")

ggplot(final_joined, aes(gene_class, map_score)) + 
  geom_boxplot() +
  geom_point(aes(color = id)) +
  geom_line(aes(group = id))


group_mean <- final_joined %>%
  group_by(id) %>% 
  summarise(mean_map = mean(map_score))
  

hist(group_mean$mean_map)

group_mean %>% 
  filter(mean_map < .78) %>% 
  count()


#Since gene pairs were split up into two - this 
gene_mappable_bed %>% 
  filter(map_score <= .80) %>% 
  summarise(n())

ggplot(gene_mappable_bed, aes(map_score)) + geom_histogram()
mappable_score <- ggplot(gene_mappable_bed, aes(map_score)) + stat_ecdf(geom = "step")
mappable_score


ggsave(path = "~/Desktop/", filename="imgs/all_genes_mappability.png", plot=mappable_score)



read_paired_named()