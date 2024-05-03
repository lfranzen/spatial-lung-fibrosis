#' Creation of new Supplementary Table with DEA results from 
#' IPF-fibrotic or IPF-alveolar vs HC-alveolar
#' together with distance correlation and p values.
#' 

library(dplyr)

DIR_DATA <- "~/Documents/PhD_Lung/spatial-lung-semla/data"
DIR_RES <- "~/Documents/PhD_Lung/spatial-lung-semla/results"

dea_res_fib <- read.csv(file.path(DIR_DATA, "hs_visium_pseudobulk_dea_res_IPF_region_fibrosis_vs_alveolar_CTRL.csv"), row.names = 1)
dea_res_alv <- read.csv(file.path(DIR_DATA, "hs_visium_pseudobulk_dea_res_IPF_region_alv_vs_alv_CTRL.csv"), row.names = 1)
cor_res <- read.csv(file.path(DIR_RES, "hs_DEA_fibrosis_distance_shared_up_DEG_expression_cor500um.csv"))

dim(cor_res)


# Create base gene table
gene_df <- dea_res_fib |> select(gene_hs, species)
colnames(gene_df) <- c("gene", "species")

# filter dea_res tables
cols_keep <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "design")

# Add column prefix for analysis
dea_res_fib <- dea_res_fib |> select(all_of(cols_keep))
colnames(dea_res_fib) <- paste0("IPF_fib_vs_HC_alv.", colnames(dea_res_fib))

dea_res_alv <- dea_res_alv |> select(all_of(cols_keep))
colnames(dea_res_alv) <- paste0("IPF_alv_vs_HC_alv.", colnames(dea_res_alv))

# Join tables
gene_df_dea <- cbind(gene_df, dea_res_fib, dea_res_alv)


# Prepare cor table
cor_res <- cor_res |> select(gene, cor, pval, pval_FDR)
colnames(cor_res) <- c("gene", "dist.PearsonCor", "dist.pval", "dist.pval_FDR")

# Join tables
gene_df_dea2 <- merge(x = gene_df_dea, y = cor_res, by = "gene", all = TRUE)
# gene_df_dea2$gene <- paste0("'", gene_df_dea2$gene)

write.csv(gene_df_dea2, file.path(DIR_RES, "SupplTable_DEA_res_Distance_cor.csv"), row.names = F)


