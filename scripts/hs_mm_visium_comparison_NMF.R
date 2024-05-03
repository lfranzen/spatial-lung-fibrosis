#' [hs_mm_visium_comparison_NMF.R]
#'
#' Comparison of NMF and regions of AbBa cells in both IPF and BLM Visium data 
#'
#' L. Franzén [lovisa.franzen@scilifelab.se]

#### Set up ####
##### Define params. ####
SPECIES <- "mouse"
SPECIES <- "human"
DIR_ROOT <- getwd()
DIR_DATA <- file.path(DIR_ROOT, "data")
DIR_RES <- file.path(DIR_ROOT, "results", "translational")
DIR_FIG <- file.path(DIR_RES, "figures")
DIR_FIG_OUT <- file.path(DIR_FIG, "radial_distance")
DIR_OBJ_OUT <- file.path(DIR_RES, "objects", "radial_distance")
dir.create(DIR_FIG_OUT); dir.create(DIR_OBJ_OUT)

##### Load libs ####
library(tidyverse)
library(dplyr)
library(tidyr)
library(STutility)
library(readxl)
library(patchwork)
library(gplots)
library(pheatmap)

library(reshape2)
library(heatmaply)
library(rlang)


##### Other ####
source(file.path(DIR_ROOT, "scripts", "custom_functions.R"))
source(file.path(DIR_ROOT, "scripts", "custom_colors.R"))
theme_custom <- theme(axis.title.x = element_blank())
fig_res <- 300

#### Gene conversion data ####
gene_conv_df <- read.csv(file.path(DIR_DATA, "misc", "orthogene_conv_combined_hs_mm.csv"), row.names = 1)
gene_conv <- setNames(gene_conv_df$symbol_mm, nm = gene_conv_df$symbol_hs)
gene_conv_mm <- setNames(gene_conv_df$symbol_hs, nm = gene_conv_df$symbol_mm)


#### NMF gene loading tables ####
n_factors <- 30

#' Create csv tables of top genes per factor
hs_nmf_genes_list <- lapply(1:n_factors, function(f){
  df <- readxl::read_xlsx(path = file.path(DIR_ROOT, "results", "human", "objects/A_NMF30", "hs_visium_A_nmf_1-30_top100_gene_loadings.xlsx"), sheet = paste0("Sheet",f))
})
hs_nmf_genes <- do.call(rbind, hs_nmf_genes_list)
write.csv(hs_nmf_genes, file.path(DIR_ROOT, "results", "human", "objects/A_NMF30", "hs_visium_A_nmf_1-30_top100_gene_loadings.csv"), row.names = F)

mm_nmf_genes_list <- lapply(1:n_factors, function(f){
  df <- readxl::read_xlsx(path = file.path(DIR_ROOT, "results", "mouse", "objects/NMF30", "mm_visium_nmf30_top100_gene_loadings.xlsx"), sheet = paste0("Sheet",f))
})
mm_nmf_genes <- do.call(rbind, mm_nmf_genes_list)
write.csv(mm_nmf_genes, file.path(DIR_ROOT, "results", "mouse", "objects/NMF30", "mm_visium_nmf30_top100_gene_loadings.csv"), row.names = F)

mm_nmf_d7_genes_list <- lapply(1:n_factors, function(f){
  df <- readxl::read_xlsx(path = file.path(DIR_ROOT, "results", "mouse", "objects/NMF30_d7", "mm_visium_nmf30_d7_top100_gene_loadings.xlsx"), sheet = paste0("Sheet",f))
})
mm_nmf_d7_genes <- do.call(rbind, mm_nmf_d7_genes_list)
write.csv(mm_nmf_d7_genes, file.path(DIR_ROOT, "results", "mouse", "objects/NMF30_d7", "mm_visium_nmf30_d7_top100_gene_loadings.csv"), row.names = F)

mm_nmf_d21_genes_list <- lapply(1:n_factors, function(f){
  df <- readxl::read_xlsx(path = file.path(DIR_ROOT, "results", "mouse", "objects/NMF30_d21", "mm_visium_nmf30_d21_top100_gene_loadings.xlsx"), sheet = paste0("Sheet",f))
})
mm_nmf_d21_genes <- do.call(rbind, mm_nmf_d21_genes_list)
write.csv(mm_nmf_d21_genes, file.path(DIR_ROOT, "results", "mouse", "objects/NMF30_d21", "mm_visium_nmf30_d21_top100_gene_loadings.csv"), row.names = F)

#' Load csv tables of top genes per factor
hs_nmf_genes <- read.csv(file.path(DIR_ROOT, "results", "human", "objects/A_NMF30", "hs_visium_A_nmf_1-30_top100_gene_loadings.csv"))
mm_nmf_genes <- read.csv(file.path(DIR_ROOT, "results", "mouse", "objects/NMF30", "mm_visium_nmf30_top100_gene_loadings.csv"))
mm_nmf_d7_genes <- read.csv(file.path(DIR_ROOT, "results", "mouse", "objects/NMF30_d7", "mm_visium_nmf30_d7_top100_gene_loadings.csv"))
mm_nmf_d21_genes <- read.csv(file.path(DIR_ROOT, "results", "mouse", "objects/NMF30_d21", "mm_visium_nmf30_d21_top100_gene_loadings.csv"))


#### NMF30 comparisons ####
# Define function to calculate Jaccard similarity index
jaccard_similarity <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  jaccard_index <- intersection / union
  return(jaccard_index)
}

ComputeJaccardSimilary <- function(file1, file2){
  # Extract relevant columns from File 1
  file1 <- file1 %>% select(gene, gene_loading, factor)
  
  # Extract relevant columns from File 2
  file2 <- file2 %>% select(gene, gene_loading, factor)
  
  # Initialize empty dataframe to store Jaccard similarity index values
  similarity_df <- data.frame()
  
  # Loop through factors in File 1
  for (factor1 in unique(file1$factor)) {
    # Filter File 1 for the current factor
    file1_factor <- file1 %>% filter(factor == factor1)
    
    # Extract gene names for the current factor in File 1
    gene_names_file1 <- file1_factor$gene
    
    # Extract gene loading values for the current factor in File 1
    gene_loadings_file1 <- file1_factor$gene_loading
    
    # Loop through factors in File 2
    for (factor2 in unique(file2$factor)) {
      # Filter File 2 for the current factor
      file2_factor <- file2 %>% filter(factor == factor2)
      
      # Extract gene names for the current factor in File 2
      gene_names_file2 <- file2_factor$gene
      
      # Extract gene loading values for the current factor in File 2
      gene_loadings_file2 <- file2_factor$gene_loading
      
      # Compute Jaccard similarity index for the current factor pair
      jaccard_index <- jaccard_similarity(gene_names_file1, gene_names_file2)
      
      # Append the Jaccard similarity index value to the similarity dataframe
      similarity_df <- rbind(similarity_df, data.frame(Factor1 = factor1, 
                                                       Factor2 = factor2, 
                                                       Jaccard_Similarity = jaccard_index))
    }
  }
  
  return(similarity_df)
}


#' Compare gene similarities between hs-NMF30 and mm-NMF30_d21
#' 
#' Done initially by Martina

#' Select only shared genes
hs_nmf_genes_subset <- subset(hs_nmf_genes, gene %in% gene_conv_df$symbol_hs)
mm_nmf_d21_genes_subset <- subset(mm_nmf_d21_genes, gene %in% gene_conv_df$symbol_mm)
# mm_nmf_genes_subset <- subset(mm_nmf_genes, gene %in% gene_conv_df$symbol_mm)

nrow(mm_nmf_d21_genes_subset);nrow(mm_nmf_d21_genes)
nrow(hs_nmf_genes_subset);nrow(hs_nmf_genes)

mm_nmf_d21_genes_subset$mm_gene <- mm_nmf_d21_genes_subset$gene
mm_nmf_d21_genes_subset$gene <- gene_conv_mm[mm_nmf_d21_genes_subset$gene]
# mm_nmf_genes_subset$human_gene <- gene_conv_mm[mm_nmf_genes_subset$gene]

# Read File 1
# file1 <- read.csv("hs_visium_nmf30_top100_gene_loadings.csv")
file1 <- hs_nmf_genes_subset

# Read File 2
# file2 <- read.csv("mm_visium_nmf30_d21_top100_gene_loadings_translated.csv")
file2 <- mm_nmf_d21_genes_subset
# file2 <- mm_nmf_genes_subset

# Run analysis
similarity_df <- ComputeJaccardSimilary(file1, file2)

# Convert similarity dataframe to wide format for creating heatmap
similarity_wide <- pivot_wider(similarity_df, names_from = Factor2, values_from = Jaccard_Similarity) %>% 
  column_to_rownames("Factor1")

similarity_wide <- t(similarity_wide)
rownames(similarity_wide) <- paste0("mm_F", rownames(similarity_wide))
colnames(similarity_wide) <- paste0("hs_F", colnames(similarity_wide))

# Convert the similarity matrix to a numeric matrix
similarity_matrix <- as.matrix(similarity_wide)
similarity_matrix_d21 <- similarity_matrix

# Plot heatmaps
similarity_matrix_plot <- similarity_matrix_d21
rownames(similarity_matrix_plot) <- gsub("mm_", "", rownames(similarity_matrix_plot))
colnames(similarity_matrix_plot) <- gsub("hs_", "", colnames(similarity_matrix_plot))

# All factors
pdf(file = file.path(DIR_RES, "figures", "NMF", "Hs_NMF30_vs_Mm_d21NMF30_jaccardHeatmap_full.pdf"), height = 8, width = 8)
pheatmap::pheatmap(similarity_matrix_plot, 
                   color = c("grey98", col_scale_acton), 
                   cellwidth = 10, cellheight = 10, 
                   border_color = NA,
                   # gaps_col = 0, gaps_row = 0, 
                   treeheight_row = 14, treeheight_col = 14)
dev.off()


# Filtered factors (Fig 4h)
# similarity_matrix_plot_filtered <- similarity_matrix_plot[rowSums(similarity_matrix_plot)>1, colSums(similarity_matrix_plot)>1]
similarity_matrix_plot_filtered <- similarity_matrix_plot[apply(similarity_matrix_plot,1, max)>0.1, apply(similarity_matrix_plot, 2, max)>0.1]

pdf(file = file.path(DIR_RES, "figures", "NMF", "Hs_NMF30_vs_Mm_d21NMF30_jaccardHeatmap_filtered.pdf"), height = 8, width = 8)
pheatmap::pheatmap(similarity_matrix_plot_filtered, 
                   color = c("grey98", col_scale_acton), 
                   cellwidth = 10, cellheight = 10,
                   border_color = NA,
                   # gaps_col = 0, gaps_row = 0, 
                   treeheight_row = 14, treeheight_col = 14)
dev.off()



#' Compare gene similarities between hs-NMF30 and mm-NMF30_d7
#' 
#' Done initially by Martina

#' Select only shared genes
hs_nmf_genes_subset <- subset(hs_nmf_genes, gene %in% gene_conv_df$symbol_hs)
mm_nmf_d7_genes_subset <- subset(mm_nmf_d7_genes, gene %in% gene_conv_df$symbol_mm)

nrow(mm_nmf_d7_genes_subset);nrow(mm_nmf_d7_genes)
nrow(hs_nmf_genes_subset);nrow(hs_nmf_genes)

mm_nmf_d7_genes_subset$mm_gene <- mm_nmf_d7_genes_subset$gene
mm_nmf_d7_genes_subset$gene <- gene_conv_mm[mm_nmf_d7_genes_subset$gene]

# Read File 1
file1 <- hs_nmf_genes_subset

# Read File 2
file2 <- mm_nmf_d7_genes_subset

# Run analysis
similarity_df <- ComputeJaccardSimilary(file1, file2)

# Convert similarity dataframe to wide format for creating heatmap
similarity_wide <- pivot_wider(similarity_df, names_from = Factor2, values_from = Jaccard_Similarity) %>% 
  column_to_rownames("Factor1")

similarity_wide <- t(similarity_wide)
rownames(similarity_wide) <- paste0("mm_F", rownames(similarity_wide))
colnames(similarity_wide) <- paste0("hs_F", colnames(similarity_wide))

# Convert the similarity matrix to a numeric matrix
similarity_matrix <- as.matrix(similarity_wide)
similarity_matrix_d7 <- similarity_matrix

# Plot heatmaps
similarity_matrix_plot <- similarity_matrix_d7
rownames(similarity_matrix_plot) <- gsub("mm_", "", rownames(similarity_matrix_plot))
colnames(similarity_matrix_plot) <- gsub("hs_", "", colnames(similarity_matrix_plot))

# All factors
pdf(file = file.path(DIR_RES, "figures", "NMF", "Hs_NMF30_vs_Mm_d7NMF30_jaccardHeatmap_full.pdf"), height = 8, width = 8)
pheatmap::pheatmap(similarity_matrix_plot, 
                   color = c("grey98", col_scale_acton), 
                   cellwidth = 10, cellheight = 10,
                   border_color = NA,
                   # gaps_col = 0, gaps_row = 0, 
                   treeheight_row = 14, treeheight_col = 14)
dev.off()


# Filtered factors
# similarity_matrix_plot_filtered <- similarity_matrix_plot[rowSums(similarity_matrix_plot)>1, colSums(similarity_matrix_plot)>1]
similarity_matrix_plot_filtered <- similarity_matrix_plot[apply(similarity_matrix_plot,1, max)>0.1, apply(similarity_matrix_plot, 2, max)>0.1]

pdf(file = file.path(DIR_RES, "figures", "NMF", "Hs_NMF30_vs_Mm_d7NMF30_jaccardHeatmap_filtered.pdf"), height = 8, width = 8)
pheatmap::pheatmap(similarity_matrix_plot_filtered, 
                   color = c("grey98", col_scale_acton), 
                   cellwidth = 10, cellheight = 10,
                   border_color = NA,
                   # gaps_col = 0, gaps_row = 0, 
                   treeheight_row = 14, treeheight_col = 14)
dev.off()



#' Concatenate gene similarities between hs-NMF30 and mm-NMF30_d21/_d7
rownames(similarity_matrix_d7) <- paste0("d7_", rownames(similarity_matrix_d7))
rownames(similarity_matrix_d21) <- paste0("d21_", rownames(similarity_matrix_d21))

similarity_matrix_concat <- rbind(similarity_matrix_d7, similarity_matrix_d21)
mm_groups <- data.frame(row.names = rownames(similarity_matrix_concat),
                        day = strsplit(rownames(similarity_matrix_concat), "_") %>% unlist() %>% grep(pattern="d", value = T))
group_colors <- list(
  day = c(d7 = "#B1E4C2FF", d21 = "#357BA2FF")
)

pdf(file = file.path(DIR_RES, "figures", "NMF", "Hs_NMF30_vs_Mm_d7NMF30_d21NMF30_jaccardHeatmap.pdf"), height = 8, width = 12)
pheatmap::pheatmap(t(similarity_matrix_concat), 
                   color = c("grey98", col_scale_acton), 
                   cellwidth = 12, 
                   cellheight = 12,
                   border_color = NA,
                   annotation_col = mm_groups, 
                   annotation_colors = group_colors,
                   treeheight_row = 10, 
                   treeheight_col = 30)
dev.off()

# Filtered factors
similarity_matrix_concat_filtered <- similarity_matrix_concat[apply(similarity_matrix_concat,1, max)>0.1, apply(similarity_matrix_concat, 2, max)>0.1]

pdf(file = file.path(DIR_RES, "figures", "NMF", "Hs_NMF30_vs_Mm_d7NMF30_d21NMF30_jaccardHeatmap_filtered.pdf"), height = 8, width = 8)
pheatmap::pheatmap(t(similarity_matrix_concat_filtered), 
                   color = c("grey98", col_scale_acton), 
                   cellwidth = 12, 
                   cellheight = 12,
                   border_color = NA,
                   annotation_col = mm_groups, 
                   annotation_colors = group_colors,
                   treeheight_row = 10, 
                   treeheight_col = 30)
dev.off()


#### hsNMF-hsNMF comparison ####
n_top_gene <- 25
hs_nmf_genes_top <- hs_nmf_genes %>% filter(rank <= n_top_gene)

# Run analysis
similarity_df <- ComputeJaccardSimilary(file1 = hs_nmf_genes_top, file2 = hs_nmf_genes_top)

# Convert similarity dataframe to wide format for creating heatmap
similarity_wide <- pivot_wider(similarity_df, names_from = Factor2, values_from = Jaccard_Similarity) %>% 
  column_to_rownames("Factor1")

similarity_wide <- t(similarity_wide)
rownames(similarity_wide) <- paste0("F", rownames(similarity_wide))
colnames(similarity_wide) <- paste0("F", colnames(similarity_wide))

# Convert the similarity matrix to a numeric matrix
similarity_matrix <- as.matrix(similarity_wide)

# Plot heatmaps
similarity_matrix_plot <- similarity_matrix

# All factors
pdf(file = file.path(DIR_RES, "figures", "NMF", "Hs_NMF30_top25_genes_jaccardHeatmap_full.pdf"), height = 8, width = 8)
pheatmap::pheatmap(similarity_matrix_plot, 
                   color = c("grey98", col_scale_acton), 
                   cellwidth = 10, cellheight = 10,
                   border_color = NA,
                   # gaps_col = 0, gaps_row = 0, 
                   treeheight_row = 14, treeheight_col = 14)
dev.off()
