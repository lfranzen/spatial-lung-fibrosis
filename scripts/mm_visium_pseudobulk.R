#' [mm_visium_pseudobulk.R]
#'
#' Preform pseudobulk analysis comparisons
#'
#'
#' Jan 2023, L. Franzén [lovisa.franzen@scilifelab.se]

#### Set up ####
##### Define params. ####
set.seed(1)
SPECIES <- "mouse"
DIR_ROOT <- "/home/st-analysis_home/lovisa.franzen/analysis/lung/spatial-lung-fibrosis"  #getwd()
DIR_DATA <- file.path(DIR_ROOT, "data", SPECIES, "visium")
DIR_RES <- file.path(DIR_ROOT, "results", SPECIES)
DIR_FIG_DEA <- file.path(DIR_RES, "figures", "DEA")
DIR_OBJ_DEA <- file.path(DIR_RES, "objects", "DEA")
dir.create(DIR_FIG_DEA)
dir.create(DIR_OBJ_DEA)
fig_res <- 300

##### Load libs ####
library(STutility)
library(RColorBrewer)
library(patchwork)
library(dplyr)
library(magrittr)
library(tidyverse)
library(writexl)
library(DESeq2)
library(UpSetR)
library(gprofiler2)


##### Other ####
# source(file.path(DIR_ROOT, "scripts", "colors.R"))
source(file.path(DIR_ROOT, "scripts", "custom_functions.R"))
source(file.path(DIR_ROOT, "scripts", "custom_colors.R"))
theme_custom <- theme(axis.title.x = element_blank())

##### Read objects ####
fname <- paste0("mm_visium_preproc_se_obj.rds")
se <- readRDS(file = file.path(DIR_RES, "objects", fname))

metadata <- se@meta.data %>% 
  select(sample_name, animal, replicate, day, condition) %>% 
  distinct(sample_name, .keep_all = TRUE) %>% 
  tibble::remove_rownames()
rownames(metadata) <- metadata$sample_name


##### Custom functions ####
GeneratePseudoBulk <- function(
    object,
    sample_column_name = "sample_name"
){
  # Get spot names for each sample
  s_names <- object@meta.data[[sample_column_name]] %>% unique()
  message(paste0("Generating pseudo-bulk data for samples ", paste(s_names, collapse = ", ")))
  spots_list <- lapply(s_names, function(s){
    object@meta.data[object@meta.data[[sample_column_name]] == s, ] %>% rownames()
  }) %>% 
    setNames(s_names)
  
  # Get count data
  count_data <- object@assays$RNA@counts
  
  # Bulk data per sample
  bulk_data_list <- lapply(s_names, function(s){
    sample_count_data <- Matrix::rowSums(count_data[, spots_list[[s]]])
    sample_bulk_data <- data.frame(sample_count_data)
    colnames(sample_bulk_data) <- s
    sample_bulk_data
  })
  
  # Join together and return
  bulk_data <- bind_cols(bulk_data_list)
  return(bulk_data)
}


#### Generate pseudobulk per animal and per region ####
##### Per histopath region #####
# Save animal metadata
mdat <- metadata %>% 
  distinct(animal, .keep_all = T) %>% 
  select(-c(sample_name, replicate)) %>% 
  remove_rownames() %>% 
  group_by(condition, day) %>% 
  mutate(rep = row_number()) %>% 
  ungroup() %>% 
  as.data.frame()

rownames(mdat) <- mdat$animal
write.csv(mdat, file.path(DIR_OBJ_DEA, "mm_visium_pseudobulk_metadata_per_animal.csv"), row.names = T)
mdat <- read.csv(file.path(DIR_OBJ_DEA, "mm_visium_pseudobulk_metadata_per_animal.csv"), row.names = 1)

# Generate pseudo-bulk for all and for each annotated region
anno_regions <- se$annotation %>% unique(); anno_regions
anno_regions <- anno_regions[!is.na(anno_regions)]

regions_test <- c("all", anno_regions)
bulk_data_list <- setNames(lapply(regions_test, function(r){
  message(paste0("Subsetting data for region ", r))
  if (r == "all") {
    bulk_dat <- GeneratePseudoBulk(se, sample_column_name = "animal")
  } else {
    se_subset <- SubsetSTData(se, annotation == r)
    bulk_dat <- GeneratePseudoBulk(se_subset, sample_column_name = "animal")
  }
  return(bulk_dat)
}), nm = regions_test)

# Export pseudobulk data per animal
for(r in names(bulk_data_list)){
  d <- bulk_data_list[[r]]
  r_fname <- gsub(" ", "", r)
  r_fname <- gsub("\\/", "", r_fname)
  message(paste0("Export data for ", r_fname))
  write.csv(d, file.path(DIR_OBJ_DEA, paste0("mm_visium_pseudobulk_data_per_animal_", r_fname, ".csv")), row.names = T)
}

##### d21-NMF30-F14hi #####
se_f14high <- readRDS(file = file.path(DIR_RES, "objects/NMF30_d21", "mm_visium_nmf30_d21_f14high_subset.rds"))

# F14hi spots
bulk_dat_f14hi <- GeneratePseudoBulk(se_f14high, sample_column_name = "animal")
write.csv(bulk_dat_f14hi, file.path(DIR_OBJ_DEA, paste0("mm_visium_pseudobulk_data_per_animal_d21-NMF30-F14hi.csv")), row.names = T)

# F14hi-C0 spots
se_f14high_C0 <- SubsetSTData(se_f14high, f14_subclusters == 0)
bulk_dat_f14hi_C0 <- GeneratePseudoBulk(se_f14high_C0, sample_column_name = "animal")
write.csv(bulk_dat_f14hi_C0, file.path(DIR_OBJ_DEA, paste0("mm_visium_pseudobulk_data_per_animal_d21-NMF30-F14hi-C0.csv")), row.names = T)


#### All: CTRL vs BLM, per day ####
bulk_all <- read.csv(file.path(DIR_OBJ_DEA, "mm_visium_pseudobulk_data_per_animal_all.csv"), row.names = 1)

# Column data
mat_all <- metadata %>% 
  distinct(animal, .keep_all = T) %>% 
  select(-c(sample_name, replicate)) %>% 
  remove_rownames() %>% 
  group_by(condition, day) %>% 
  mutate(rep = row_number()) %>% 
  ungroup() %>% 
  as.data.frame()
rownames(mat_all) <- mat_all$animal

##### DEA: Day 7 ##### 
de_des_d7 <- subset(mat_all, day == "d7") %>% select(-c(animal, day))
de_des_d7$condition <- relevel(factor(de_des_d7$condition), ref = "control")
d_d7 <- bulk_all[, rownames(de_des_d7)]

# DESeq2
dds_all_d7 <- DESeqDataSetFromMatrix(countData = d_d7,
                                         colData = de_des_d7,
                                         design = ~ condition)
dds_all_d7 <- DESeq(dds_all_d7)

# Results
# resultsNames(dds_all_d7) # lists the coefficients
res_all_d7 <- results(dds_all_d7, name="condition_bleomycin_vs_control")
plotMA(res_all_d7, ylim=c(-3,3), main = "MA-plot")

res_all_d7_df <- as.data.frame(res_all_d7)
res_all_d7_sign <- res_all_d7_df %>% 
  dplyr::arrange(padj, abs(log2FoldChange)) %>%
  dplyr::filter(padj<0.01)


##### DEA: Day 21 ##### 
de_des_d21 <- subset(mat_all, day == "d21") %>% select(-c(animal, day))
de_des_d21$condition <- relevel(factor(de_des_d21$condition), ref = "control")
d_d21 <- bulk_all[, rownames(de_des_d21)]

# DESeq2
dds_all_d21 <- DESeqDataSetFromMatrix(countData = d_d21,
                                          colData = de_des_d21,
                                          design = ~ condition)
dds_all_d21 <- DESeq(dds_all_d21)

# Results
# resultsNames(dds_all_d21) # lists the coefficients
res_all_d21 <- results(dds_all_d21, name="condition_bleomycin_vs_control")
plotMA(res_all_d21, ylim=c(-3,3), main = "MA-plot")

res_all_d21_df <- as.data.frame(res_all_d21)
res_all_d21_sign <- res_all_d21_df %>% 
  dplyr::arrange(padj, abs(log2FoldChange)) %>%
  dplyr::filter(padj<0.01)


######  Export results ###### 
write.csv(res_all_d7_df, file.path(DIR_OBJ_DEA, "mm_visium_pseudobulk_anno_all_day7_deseq2_res.csv"))
write.csv(res_all_d21_df, file.path(DIR_OBJ_DEA, "mm_visium_pseudobulk_anno_all_day21_deseq2_res.csv"))

res_all_d7_df <- read.csv(file.path(DIR_OBJ_DEA, "mm_visium_pseudobulk_anno_all_day7_deseq2_res.csv"), row.names = 1) %>% 
  rownames_to_column(var = "gene") %>% 
  arrange(desc(abs(log2FoldChange)))
res_all_d21_df <- read.csv(file.path(DIR_OBJ_DEA, "mm_visium_pseudobulk_anno_all_day21_deseq2_res.csv"), row.names = 1) %>% 
  rownames_to_column(var = "gene") %>% 
  arrange(desc(abs(log2FoldChange)))

dea_res_export <- setNames(list(subset(res_all_d7_df, padj<0.05), 
                                subset(res_all_d21_df, padj<0.05)), 
                           nm = c("d7_all_degs", "d21_all_degs"))

fname <- "mm_visium_pseudobulk_anno_all_day_ctrl_vs_blm_degs"
write_xlsx(
  x = dea_res_export,
  path = file.path(DIR_OBJ_DEA, paste0(fname, ".xlsx")),
  col_names = TRUE,
  format_headers = TRUE,
  use_zip64 = FALSE
)



#### Alveoli: CTRL vs BLM, per day ####

##### Prep data ##### 
# se_alveoli <- SubsetSTData(se, annotation == "Normal Alveolar and Other")
# bulk_alveoli <- GeneratePseudoBulk(se_alveoli, sample_column_name = "animal") # Generate bulk data (ignore replicate)
bulk_alveoli <- read.csv(file.path(DIR_OBJ_DEA, "mm_visium_pseudobulk_data_per_animal_NormalAlveolarandOther.csv"), row.names = 1)

# Column data
mat_alv <- metadata %>% 
  distinct(animal, .keep_all = T) %>% 
  select(-c(sample_name, replicate)) %>% 
  remove_rownames() %>% 
  group_by(condition, day) %>% 
  mutate(rep = row_number()) %>% 
  ungroup() %>% 
  as.data.frame()

rownames(mat_alv) <- mat_alv$animal


##### DEA: Day 7 ##### 
de_des_d7 <- subset(mat_alv, day == "d7") %>% select(-c(animal, day))
de_des_d7$condition <- relevel(factor(de_des_d7$condition), ref = "control")
d_d7 <- bulk_alveoli[, rownames(de_des_d7)]

# DESeq2
dds_alveoli_d7 <- DESeqDataSetFromMatrix(countData = d_d7,
                                         colData = de_des_d7,
                                         design = ~ condition)
dds_alveoli_d7 <- DESeq(dds_alveoli_d7)

# Results
# resultsNames(dds_alveoli_d7) # lists the coefficients
res_alveoli_d7 <- results(dds_alveoli_d7, name="condition_bleomycin_vs_control")
plotMA(res_alveoli_d7, ylim=c(-3,3), main = "MA-plot")

res_alveoli_d7_df <- as.data.frame(res_alveoli_d7)
res_alveoli_d7_sign <- res_alveoli_d7_df %>% 
  dplyr::arrange(padj, abs(log2FoldChange)) %>%
  dplyr::filter(padj<0.01)


##### DEA: Day 21 ##### 
de_des_d21 <- subset(mat_alv, day == "d21") %>% select(-c(animal, day))
de_des_d21$condition <- relevel(factor(de_des_d21$condition), ref = "control")
d_d21 <- bulk_alveoli[, rownames(de_des_d21)]

# DESeq2
dds_alveoli_d21 <- DESeqDataSetFromMatrix(countData = d_d21,
                                         colData = de_des_d21,
                                         design = ~ condition)
dds_alveoli_d21 <- DESeq(dds_alveoli_d21)

# Results
# resultsNames(dds_alveoli_d21) # lists the coefficients
res_alveoli_d21 <- results(dds_alveoli_d21, name="condition_bleomycin_vs_control")
plotMA(res_alveoli_d21, ylim=c(-3,3), main = "MA-plot")

res_alveoli_d21_df <- as.data.frame(res_alveoli_d21)
res_alveoli_d21_sign <- res_alveoli_d21_df %>% 
  dplyr::arrange(padj, abs(log2FoldChange)) %>%
  dplyr::filter(padj<0.01)


######  Export results ###### 
dea_res_export <- setNames(list(subset(res_alveoli_d7_df, padj<0.05), 
                                subset(res_alveoli_d21_df, padj<0.05)), 
                           nm = c("d7_alveoli_degs", "d21_alveoli_degs"))

fname <- "mm_visium_pseudobulk_anno_alveoli_day_ctrl_vs_blm_degs"
write_xlsx(
  x = dea_res_export,
  path = file.path(DIR_OBJ_DEA, paste0(fname, ".xlsx")),
  col_names = TRUE,
  format_headers = TRUE,
  use_zip64 = FALSE
)

write.csv(res_alveoli_d7_df, file.path(DIR_OBJ_DEA, "mm_visium_pseudobulk_anno_alveoli_day7_deseq2_res.csv"))
write.csv(res_alveoli_d21_df, file.path(DIR_OBJ_DEA, "mm_visium_pseudobulk_anno_alveoli_day21_deseq2_res.csv"))



##### Plot results ##### 

###### Day 7 ######

####### Box plot top DEGs ####### 
# Box plot of top res_alveoli_d7_sign
d7_g_plot <- bind_rows(
  res_alveoli_d7_df %>% 
    dplyr::filter(padj<0.01, log2FoldChange > 2) %>% 
    dplyr::arrange(desc(log2FoldChange), padj) %>%
    head(10),
  res_alveoli_d7_df %>% 
    dplyr::filter(padj<0.01, log2FoldChange < -2) %>% 
    dplyr::arrange((log2FoldChange), padj) %>%
    head(10),
)

g_plot <- rownames(d7_g_plot)
d_plot <- t(d_d7[g_plot, ])
d_plot <- cbind(d_plot, de_des_d7)
plot_list <- lapply(g_plot, function(g){
  p <- ggplot(d_plot, aes_string(x = "condition", y = paste0("log2(", g,"+1)"), fill = "condition")) +
    geom_boxplot() +
    geom_point() +
    labs(y="log2(count+1)", title=g) +
    scale_fill_manual(values = cols_cond) +
    theme_classic() +
    theme(legend.position = "none", 
          axis.title.x = element_blank(), 
          axis.text = element_text(size=10),
          plot.title = element_text(hjust=0.5, size=10, face = "bold"))
})

p_box_top_d7 <- wrap_plots(plot_list, ncol = 5) + 
  patchwork::plot_annotation(title = "Top DEGs in day 7: CTRL vs BLM", theme = theme(plot.title = element_text(hjust=0.5)))


######## Volcano ####### 
res_alveoli_d7_df$significant_0.05 <- res_alveoli_d7_df$padj < 0.05
res_alveoli_d7_df$sign_up_down <- ifelse(test = (res_alveoli_d7_df$significant_0.05 & res_alveoli_d7_df$log2FoldChange > 0),
                                          yes = "sign_up",
                                          no = ifelse((res_alveoli_d7_df$significant_0.05 & res_alveoli_d7_df$log2FoldChange < 0),
                                                      yes = "sign_down", no = NA))
cols_sign_up_down <- setNames(c("#D7191C", "#2B83BA", "grey40"),
                              nm = c("sign_up", "sign_down", "not_sign"))

res_alveoli_d7_df$gene_labels <- ifelse((res_alveoli_d7_df$padj < 0.01 & abs(res_alveoli_d7_df$log2FoldChange) > 5) |
                                          (-log10(res_alveoli_d7_df$padj) > 40),
                                   rownames(res_alveoli_d7_df), ""); unique(res_alveoli_d7_df$gene_labels) %>% length()

# glabs_d7 <- res_alveoli_d7_df[res_alveoli_d7_df$gene_labels != "" & !is.na(res_alveoli_d7_df$gene_labels), "gene_labels"]

minmax_xilm <- round( max(abs(res_alveoli_d7_df$log2FoldChange)) + 1 )

p_volcano_d7 <- ggplot(res_alveoli_d7_df, aes(x = log2FoldChange, y = -log10(padj), color = sign_up_down)) +
  # geom_hline(yintercept = -log10(0.05), linetype="dashed", color = colors_multi[3], alpha=.8) +
  geom_hline(yintercept = 0, linetype="solid", color = "grey90") +
  geom_hline(yintercept = 2, linetype="dashed", color = "grey90") +
  geom_vline(xintercept = c(-2, 2), linetype="dashed", color = "grey90") +
  geom_point(size = .5, alpha = 0.8) +
  ggrepel::geom_text_repel(data = subset(res_alveoli_d7_df, gene_labels != "" & !is.na(gene_labels)), 
                           mapping = aes(x = log2FoldChange, y = -log10(padj), label = gene_labels), 
                           color = "grey20", 
                           size = 3, 
                           max.overlaps = 20) +
  # scale_color_manual(values = c("grey20", "orange")) +
  scale_color_manual(values = cols_sign_up_down) +
  xlim(c(-minmax_xilm, minmax_xilm)) +
  labs(title = "Day 7 DEGs CTRL vs BLM") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1), hjust = 0.5),
        axis.title = element_text(size = rel(1)),
        axis.line = element_line(colour = "black"))

p_volcano_d7


###### Day 21 ######

####### Box plot top DEGs ####### 
# Box plot of top res_alveoli_d21_sign
d21_g_plot <- bind_rows(
  res_alveoli_d21_df %>% 
    dplyr::filter(padj<0.01, log2FoldChange > 2) %>% 
    dplyr::arrange(desc(log2FoldChange), padj) %>%
    head(10),
  res_alveoli_d21_df %>% 
    dplyr::filter(padj<0.01, log2FoldChange < -2) %>% 
    dplyr::arrange((log2FoldChange), padj) %>%
    head(10),
)

g_plot <- rownames(d21_g_plot)
d_plot <- t(d_d21[g_plot, ])
d_plot <- cbind(d_plot, de_des_d21)
plot_list <- lapply(g_plot, function(g){
  p <- ggplot(d_plot, aes_string(x = "condition", y = paste0("log2(", g,"+1)"), fill = "condition")) +
    geom_boxplot() +
    geom_point() +
    labs(y="log2(count+1)", title=g) +
    scale_fill_manual(values = cols_cond) +
    theme_classic() +
    theme(legend.position = "none", 
          axis.title.x = element_blank(), 
          axis.text = element_text(size=10),
          plot.title = element_text(hjust=0.5, size=10, face = "bold"))
})
p_box_top_d21 <- wrap_plots(plot_list, ncol = 5) + 
  patchwork::plot_annotation(title = "Top DEGs in day 21: CTRL vs BLM", theme = theme(plot.title = element_text(hjust=0.5)))


####### Volcano ####### 
res_alveoli_d21_df$significant_0.05 <- res_alveoli_d21_df$padj < 0.05
res_alveoli_d21_df$sign_up_down <- ifelse(test = (res_alveoli_d21_df$significant_0.05 & res_alveoli_d21_df$log2FoldChange > 0),
                                          yes = "sign_up",
                                          no = ifelse((res_alveoli_d21_df$significant_0.05 & res_alveoli_d21_df$log2FoldChange < 0),
                                                      yes = "sign_down", no = NA)
                                          )
cols_sign_up_down <- setNames(c("#D7191C", "#2B83BA", "grey40"),
                              nm = c("sign_up", "sign_down", "not_sign"))

res_alveoli_d21_df$gene_labels <- ifelse((res_alveoli_d21_df$padj < 0.01 & abs(res_alveoli_d21_df$log2FoldChange) > 5) |
                                           (-log10(res_alveoli_d21_df$padj) > 40),
                                         rownames(res_alveoli_d21_df), ""); unique(res_alveoli_d21_df$gene_labels) %>% length()

minmax_xilm <- round( max(abs(res_alveoli_d21_df$log2FoldChange)) + 1 )

p_volcano_d21 <- ggplot(res_alveoli_d21_df, aes(x = log2FoldChange, y = -log10(padj), color = sign_up_down)) +
  geom_hline(yintercept = 0, linetype="solid", color = "grey90") +
  geom_hline(yintercept = 2, linetype="dashed", color = "grey90") +
  geom_vline(xintercept = c(-2, 2), linetype="dashed", color = "grey90") +
  geom_point(size = .5, alpha = 0.8) +
  ggrepel::geom_text_repel(data = subset(res_alveoli_d21_df, gene_labels != "" & !is.na(gene_labels)),
                           mapping = aes(x = log2FoldChange, y = -log10(padj), label = gene_labels),
                           color = "grey20",
                           size = 3,
                           max.overlaps = 20) +
  # scale_color_manual(values = c("grey20", "orange")) +
  scale_color_manual(values = cols_sign_up_down) +
  xlim(c(-minmax_xilm, minmax_xilm)) +
  labs(title = "Day 21 DEGs CTRL vs BLM") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1), hjust = 0.5),
        axis.title = element_text(size = rel(1)),
        axis.line = element_line(colour = "black"), 
        aspect.ratio = 1)

p_volcano_d21


###### Export plots ######
pdf(file = file.path(DIR_FIG_DEA, "mm_visium_pseudobulk_anno_alveoli_day_ctrl_vs_blm.pdf"), width = 10, height = 7)
p_box_top_d7
p_box_top_d21
p_volcano_d7 | p_volcano_d21
dev.off()


##### Day 7 vs Day 21 overlap ##### 
# Cutoffs: padj < 0.01, abs(logfc) > 1
sign_gene_list <- setNames(
  lapply(list(res_alveoli_d7_df, res_alveoli_d21_df), function(d){
    subset(d, padj < 0.01 & abs(log2FoldChange) > 1) %>% rownames()
  }), 
  nm = c("d7", "d21")
)

###### Upset plot ###### 
upset(fromList(sign_gene_list), order.by = "freq")

# list unique and shared genes
unique_d7 <- setdiff(sign_gene_list[["d7"]], sign_gene_list[["d21"]])
unique_d21 <- setdiff(sign_gene_list[["d21"]], sign_gene_list[["d7"]])
shared_d7d21 <- intersect(sign_gene_list[["d7"]], sign_gene_list[["d21"]])


######  gProfiler ###### 
library(gprofiler2)

org <- "mmusculus"
fname <- "mm_visium_pseudobulk_anno_alveoli_day_ctrl_vs_blm_deg_overlap_gProfiler"
gene_query_list <- list(unique_d7, unique_d21, shared_d7d21)
names(gene_query_list) <- c("unique_d7", "unique_d21", "shared_d7d21")
gostres_list <- list()

for (f in names(gene_query_list)){
  message(f)
  gostres_list[[f]] <- gprofiler2::gost(query = gene_query_list[[f]],
                                        organism = org, 
                                        ordered_query = F)
  }

gostres_export_list <- list()
for (f in names(gostres_list)){
  gostres_export_list[[f]] <- gostres_list[[f]]$result
  gostres_export_list[[f]]$query <- f
}
write_xlsx(
  x = gostres_export_list,
  path = file.path(DIR_OBJ_DEA, paste0(fname, ".xlsx")),
  col_names = TRUE,
  format_headers = TRUE,
  use_zip64 = FALSE
)

pdf(file = file.path(DIR_FIG_DEA, paste0(fname, ".pdf")), width = 10, height = 8, useDingbats = F)
for(f in names(gostres_list)){
  message(f)
  top_sign_terms <- gostres_list[[f]]$result %>%
    group_by(source) %>%
    arrange(p_value) %>%
    top_n(n=2, wt = -p_value)
  p <- gostplot(gostres_list[[f]], capped = F, interactive = F) +
    labs(title = f) +
    theme(plot.title = element_text(vjust = -0.5))
  publish_gostplot(p, highlight_terms = top_sign_terms$term_id)
}
dev.off()





