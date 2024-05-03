#' [hs_mm_visium_comparison_gene_conversion.R]
#'
#' Comparison of regions of TLS in both IPF and BLM Visium data 
#'
#'
#' Mar-Apr 2023, L. Franzén [lovisa.franzen@scilifelab.se]

#### Set up ####
##### Define params. ####
DIR_ROOT <- getwd()
DIR_DATA <- file.path(DIR_ROOT, "data")

##### Load libs ####
library(tidyverse)
library(dplyr)

#### Gene conversion data ####
#' Generated with package `orthogene`, using the function (orthogene_1.0.2)
#' `convert_orthologs(input_species = "human", output_species = "mouse")`,
#' with the list `hs_visium_gene_list.csv` as input (generated from hs visium
#' count assay rownames). Generated 2023-03-29.
#' Genes remaining: 13,330
gene_conv_df_hs_mm <- read.csv(file.path(DIR_DATA, "misc", "orthogene_conv_hs_to_mm.csv"))  # genes: 13330

#' Generated with package `orthogene`, using the function (orthogene_1.0.2)
#' `convert_orthologs(input_species = "mouse", output_species = "human")`,
#' with the list `mm_visium_gene_list.csv` as input (generated from mm visium
#' count assay rownames). Generated 2023-04-06.
#' Total genes dropped after convert_orthologs :
#'  1,568 / 15,362 (10%)
#' Total genes remaining after convert_orthologs :
#'  13,794 / 15,362 (90%)
gene_conv_df_mm_hs <- read.csv(file.path(DIR_DATA, "misc", "orthogene_conv_mm_to_hs.csv"))

#' Combine lists
gene_conv_df_comb <- bind_rows(gene_conv_df_hs_mm, gene_conv_df_mm_hs); dim(gene_conv_df_comb)
gene_conv_df_comb$symbol_hs_mm <- paste0(gene_conv_df_comb$symbol_hs, "_", gene_conv_df_comb$symbol_mm)
gene_conv_df_comb <- gene_conv_df_comb[!duplicated(gene_conv_df_comb$symbol_hs_mm),]; dim(gene_conv_df_comb)

# Save output table
write.csv(gene_conv_df_comb, file.path(DIR_DATA, "misc", "orthogene_conv_combined_hs_mm.csv"))
