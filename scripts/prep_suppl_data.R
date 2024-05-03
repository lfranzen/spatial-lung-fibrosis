#' Prep Supplementary Tables
#' 
#' Script for prepping of additional supplementary files

SPECIES <- "mouse"
SPECIES <- "human"

DIR_ROOT <- getwd()
DIR_DATA <- file.path(DIR_ROOT, "data")
DIR_RES <- file.path(DIR_ROOT, "results", "misc")
DIR_OBJ_OUT <- file.path(DIR_RES, "objects")
dir.create(DIR_OBJ_OUT)

##### Load libs ####
library(tidyverse)
library(writexl)

# NMF Top 100 gene loadings
n_factors <- 30

nmf_hs_list <- lapply(1:n_factors, function(f){
  df <- readxl::read_xlsx(path = file.path(DIR_ROOT, "results", "human", "objects/A_NMF30", "hs_visium_A_nmf_1-30_top100_gene_loadings.xlsx"), sheet = paste0("Sheet",f))
})
nmf_hs <- do.call(rbind, nmf_hs_list)
nmf_hs <- nmf_hs[, c("factor", "gene", "gene_loading", "gene_loading_scaled", "rank")]

nmf_mm_all_list <- lapply(1:n_factors, function(f){
  df <- readxl::read_xlsx(path = file.path(DIR_ROOT, "results", "mouse", "objects/NMF30", "mm_visium_nmf30_top100_gene_loadings.xlsx"), sheet = paste0("Sheet",f))
})
nmf_mm_all <- do.call(rbind, nmf_mm_all_list)
nmf_mm_all <- nmf_mm_all[, c("factor", "gene", "gene_loading", "rank")]

nmf_mm_d7_list <- lapply(1:n_factors, function(f){
  df <- readxl::read_xlsx(path = file.path(DIR_ROOT, "results", "mouse", "objects/NMF30_d7", "mm_visium_nmf30_d7_top100_gene_loadings.xlsx"), sheet = paste0("Sheet",f))
})
nmf_mm_d7 <- do.call(rbind, nmf_mm_d7_list)
nmf_mm_d7 <- nmf_mm_d7[, c("factor", "gene", "gene_loading", "rank")]

nmf_mm_d21_list <- lapply(1:n_factors, function(f){
  df <- readxl::read_xlsx(path = file.path(DIR_ROOT, "results", "mouse", "objects/NMF30_d21", "mm_visium_nmf30_d21_top100_gene_loadings.xlsx"), sheet = paste0("Sheet",f))
})
nmf_mm_d21 <- do.call(rbind, nmf_mm_d21_list)
nmf_mm_d21 <- nmf_mm_d21[, c("factor", "gene", "gene_loading", "rank")]

nmf_data_out <- list(
  Hs_NMF30 = nmf_hs,
  Mm_NMF30_all = nmf_mm_all,
  Mm_NMF30_d7 = nmf_mm_d7,
  Mm_NMF30_d21 = nmf_mm_d21
)

writexl::write_xlsx(nmf_data_out, path = file.path(DIR_OBJ_OUT, "NMF_top_100_genes.xlsx"))

# Save version with all genes concatenated
nmf_hs$gene_loading_scaled <- NULL
nmf_hs$subset <- "Human (Hs) HC and IPF"
nmf_mm_all$subset <- "Mouse (Mm) veh control and BLM treated (all)"
nmf_mm_d7$subset <- "Mouse (Mm) veh control and BLM treated at day 7 (d7)"
nmf_mm_d21$subset <- "Mouse (Mm) veh control and BLM treated at day 21 (d7)"

nmf_master <- bind_rows(nmf_hs, nmf_mm_all, nmf_mm_d7, nmf_mm_d21)

writexl::write_xlsx(nmf_master, path = file.path(DIR_OBJ_OUT, "NMF_top_100_genes_concatenated.xlsx"))
write.csv(nmf_master, file = file.path(DIR_OBJ_OUT, "NMF_top_100_genes_concatenated.csv"))

