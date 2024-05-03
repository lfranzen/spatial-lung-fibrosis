#' [hs_visium_pseudobulk.R]
#'
#' Preform pseudobulk analysis comparisons
#'
#'
#' Jan 2023, L. Franzén [lovisa.franzen@scilifelab.se]

#### Set up ####
##### Define params. ####
set.seed(1)
SPECIES <- "human"
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
fname <- paste0("hs_visium_preproc_A_se_obj.rds")
se <- readRDS(file = file.path(DIR_RES, "objects", fname))

metadata <- se@meta.data %>% 
  select(sample_name, subject_alias, replicate, tissue_alias, fibrotic_extent_score_by_pathologist_0.3, condition) %>% 
  distinct(sample_name, .keep_all = TRUE) %>% 
  remove_rownames()
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


#### Generate pseudobulk per donor and per region ####
# Save donor metadata
mdat <- metadata %>% 
  distinct(subject_alias, .keep_all = T) %>% 
  select(-c(sample_name, replicate, fibrotic_extent_score_by_pathologist_0.3, tissue_alias)) %>% 
  remove_rownames() %>% 
  group_by(condition) %>% 
  mutate(rep = row_number()) %>% 
  ungroup() %>% 
  arrange(condition) %>% 
  as.data.frame()
rownames(mdat) <- mdat$subject_alias
write.csv(mdat, file.path(DIR_OBJ_DEA, "hs_visium_pseudobulk_metadata_per_donor.csv"), row.names = T)

# Generate pseudo-bulk for all and for each annotated region
anno_regions <- se$annotation %>% unique(); anno_regions
anno_regions <- anno_regions[!is.na(anno_regions)]

regions_test <- c(anno_regions)
bulk_data_list <- setNames(lapply(regions_test, function(r){
  message(paste0("Subsetting data for region ", r))
  if (r == "all") {
    bulk_dat <- GeneratePseudoBulk(se, sample_column_name = "subject_alias")
  } else {
    se_subset <- SubsetSTData(se, annotation == r)
    bulk_dat <- GeneratePseudoBulk(se_subset, sample_column_name = "subject_alias")
  }
  bulk_dat <- bulk_dat[, colnames(bulk_dat) %>% sort()]
  return(bulk_dat)
}), nm = regions_test)

# Export pseudobulk data per donor
for(r in names(bulk_data_list)){
  d <- bulk_data_list[[r]]
  r_fname <- gsub(" ", "", r)
  r_fname <- gsub("\\/", "", r_fname)
  message(paste0("Export data for ", r_fname))
  write.csv(d, file.path(DIR_OBJ_DEA, paste0("hs_visium_pseudobulk_data_per_donor_", r_fname, ".csv")), row.names = T)
}


##### NMF30-F14hi #####
se_f14high <- readRDS(file = file.path(DIR_RES, "objects", "hs_visium_preproc_A_se_obj_f14high_subset.rds"))

# F14hi spots
se_f14high <- SubsetSTData(se_f14high, condition == "IPF")
bulk_dat_f14hi <- GeneratePseudoBulk(se_f14high, sample_column_name = "subject_alias")
write.csv(bulk_dat_f14hi, file.path(DIR_OBJ_DEA, paste0("hs_visium_pseudobulk_data_per_donor_NMF30-F14hi.csv")), row.names = T)

# F14hi-C0 spots
se_f14high_C0 <- SubsetSTData(se_f14high, f14_subclusters == 0)
bulk_dat_f14hi_C0 <- GeneratePseudoBulk(se_f14high_C0, sample_column_name = "subject_alias")
write.csv(bulk_dat_f14hi_C0, file.path(DIR_OBJ_DEA, paste0("hs_visium_pseudobulk_data_per_donor_NMF30-F14hi-C0.csv")), row.names = T)



#### HC vs IPF ####
# Bulk data per donor

# Generate bulk data (ignore replicate)
bulk_all <- GeneratePseudoBulk(se, sample_column_name = "subject_alias")
bulk_all <- bulk_all[, c(paste0("HC_", 1:4), paste0("IPF_", 1:4))]

write.csv(bulk_all, file = file.path(DIR_OBJ_DEA, "hs_visium_pseudobulk_data_grouped_by_donor_all.csv"))

# Column data
mat_all <- metadata %>% 
  distinct(subject_alias, .keep_all = T) %>% 
  select(-c(sample_name, replicate, fibrotic_extent_score_by_pathologist_0.3, tissue_alias)) %>% 
  remove_rownames() %>% 
  group_by(condition) %>% 
  mutate(rep = row_number()) %>% 
  ungroup() %>% 
  as.data.frame()
rownames(mat_all) <- mat_all$subject_alias
mat_all$condition <- relevel(factor(mat_all$condition), ref = "control")
mat_all <- mat_all[colnames(bulk_all), ]


##### DEA: HC vs IPF, all donors##### 
# de_des_all <- mat_all %>% select(-c(subject_alias, day))
# de_des_d7$condition <- relevel(factor(de_des_d7$condition), ref = "control")
# d_d7 <- bulk_alveoli[, rownames(de_des_d7)]

# DESeq2
dds_all <- DESeqDataSetFromMatrix(countData = bulk_all,
                                  colData = mat_all,
                                  design = ~ condition)
dds_all <- DESeq(dds_all)


# Results
resultsNames(dds_all) # lists the coefficients
res_all <- results(dds_all, name="condition_IPF_vs_control")
plotMA(res_all, ylim=c(-3,3), main = "MA-plot")

res_all_df <- as.data.frame(res_all)
res_all_sign <- res_all_df %>% 
  dplyr::arrange(padj, abs(log2FoldChange)) %>%
  dplyr::filter(padj<0.01)

write.csv(res_all_df, file.path(DIR_OBJ_DEA, "hs_visium_pseudobulk_data_grouped_by_donor_deseq2_res.csv"))


####### Box plot top DEGs ####### 
# Box plot of top res_alveoli_d21_sign
g_plot <- bind_rows(
  res_all_df %>% 
    dplyr::filter(padj<0.01, log2FoldChange > 2) %>% 
    dplyr::arrange(desc(log2FoldChange), padj) %>%
    head(10),
  res_all_df %>% 
    dplyr::filter(padj<0.01, log2FoldChange < -2) %>% 
    dplyr::arrange((log2FoldChange), padj) %>%
    head(10),
  ) %>% 
  rownames(g_plot)

d_plot <- t(bulk_all[g_plot, ])
d_plot <- cbind(d_plot, mat_all)

plot_list <- lapply(g_plot, function(g){
  p <- ggplot(d_plot, aes_string(x = "condition", y = paste0("log2(", g,"+1)"), fill = "condition")) +
    geom_boxplot(color="black", size=0.2) +
    geom_point(size=1) +
    labs(y="log2(count+1)", title=g) +
    scale_fill_manual(values = cols_cond) +
    theme_classic() +
    theme(legend.position = "none", 
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size=8),
          axis.text = element_text(size=8),
          plot.title = element_text(hjust=0.5, size=8, face = "bold"))
})
p_box_top_all <- wrap_plots(plot_list, ncol = 5) + 
  patchwork::plot_annotation(title = "Top DEGs: HC vs IPF", theme = theme(plot.title = element_text(hjust=0.5)))

pdf(file = file.path(DIR_FIG_DEA, "hs_visium_pseudobulk_data_grouped_by_donor_deseq2_res_top10_boxplot.pdf"), width = 6, height = 6)
p_box_top_all
dev.off()


######## Fig. 1D: DEA Volcano ####### 
# res_all_df <- read.csv(file.path(DIR_OBJ_DEA, "hs_visium_pseudobulk_data_grouped_by_donor_deseq2_res.csv"), row.names = 1)

res_all_df$significant_0.05 <- res_all_df$padj < 0.05
res_all_df$sign_up_down <- ifelse(test = (res_all_df$significant_0.05 & res_all_df$log2FoldChange > 0),
                                         yes = "sign_up",
                                         no = ifelse((res_all_df$significant_0.05 & res_all_df$log2FoldChange < 0),
                                                     yes = "sign_down", no = NA))
cols_sign_up_down <- setNames(c(cols_cond[["IPF"]], cols_cond[["control"]], "grey40"),
                              nm = c("sign_up", "sign_down", "not_sign"))

res_all_df$gene_labels <- ifelse((res_all_df$padj < 0.0001 & abs(res_all_df$log2FoldChange) > 6) |
                                   (-log10(res_all_df$padj) > 15),
                                 rownames(res_all_df), ""); unique(res_all_df$gene_labels) %>% length()

res_all_df$negLog10_padj_capped <- ifelse(-log10(res_all_df$padj)>20, 20, -log10(res_all_df$padj))
minmax_xilm <- round( max(abs(res_all_df$log2FoldChange)) + 1 )

p_volcano_all <- ggplot(res_all_df, aes(x = log2FoldChange, y = negLog10_padj_capped, fill = sign_up_down, color = sign_up_down)) +
  geom_hline(yintercept = c(0, 20), linetype="dashed", color = "grey") +
  geom_point(size = .5, alpha = 0.8, shape=20) +
  scale_color_manual(values = cols_sign_up_down) +
  xlim(c(-minmax_xilm, minmax_xilm)) +
  labs(title = "DEGs HC vs IPF", x="log2(fold-change)", y="-log10(adj. p-value") +
  theme_bw() +
  theme(legend.position = "none",
        aspect.ratio = 1,
        panel.grid = element_blank(), 
        panel.border = element_rect(color = "black"),
        axis.title = element_text(size=8, color = "black"),
        axis.text = element_text(size=8, color = "black"),
        plot.title = element_text(hjust=0.5, size=8))

p_volcano_all_txt <- p_volcano_all +
  ggrepel::geom_text_repel(data = subset(res_all_df, gene_labels != "" & !is.na(gene_labels)),
                           mapping = aes(x = log2FoldChange, y = negLog10_padj_capped, label = gene_labels),
                           color = "black", 
                           size = 2,
                           box.padding = 0.8,
                           max.overlaps = Inf)

p_volcano_all <- ggrastr::rasterize(p_volcano_all, layers = "Point", dpi = 600)
p_volcano_all_txt <- ggrastr::rasterize(p_volcano_all_txt, layers = "Point", dpi = 600)

pdf(file = file.path(DIR_FIG_DEA, "hs_visium_pseudobulk_data_grouped_by_donor_deseq2_res_volcano.pdf"), width = 2, height = 2)  # !Main Figure 3c
p_volcano_all
p_volcano_all_txt
dev.off()


####### Spatial DEGs ####### 
head(res_all_df)
res_all_df$gene <- rownames(res_all_df)

# Genes with text labels in volcano plot
genes_plot <- subset(res_all_df, gene_labels != "" & !is.na(gene_labels))$gene_labels

se$sample_name_new <- gsub("G", "B", se$sample_name)


# subset data
se_plot_subset <- SubsetSTData(se, annotation != "NA")

se_plot_subset$sample_name %>% unique()
samples_keep <- c("HC_1.TH010.B0.1", 
                  "IPF_1.TD012.B2.1", 
                  "IPF_2.TD021B.B2.2", 
                  "IPF_3.TD032A.B2.2",
                  "IPF_4.TD041B.B2.1")
se_plot_subset <- SubsetSTData(se, sample_name %in% samples_keep)

se_plot_subset$sample_name_new2 <- paste0(se_plot_subset$subject_alias, ".B", 
                                          se_plot_subset$grade_tissue_selection,
                                          ".", se_plot_subset$replicate)

# Plot
genes_plot_selected <- c("COL1A1", "COL3A1", "THY1", "SFRP2", "MMP11")


p1 <- ST.FeaturePlot(se_plot_subset, features = "annotation", label.by = "sample_name_new2", 
                     cols = cols_annotation2,
                     ncol = 1) & 
  theme(aspect.ratio = 1, legend.position = "left", 
        plot.title = element_text(hjust=0.5), legend.title = element_blank()) &
  guides(fill = guide_legend(override.aes = list(size = 3)));p1


p2 <- ST.FeaturePlot(se_plot_subset, features = genes_plot_selected, label.by = "sample_name_new2", 
               ncol = 1, cols = col_scale_or, min.cutoff = 0.5, pt.border = F) & 
  theme(aspect.ratio = 1, legend.position = "top", plot.title = element_text(hjust=0.5));p2

pdf(file = file.path(DIR_FIG_DEA, "hs_visium_pseudobulk_data_grouped_by_donor_deseq2_res_spatial_selected.pdf"), 
    width = 15, height = 12)
(p1|p2)+plot_layout(widths = c(1,length(genes_plot_selected)))
dev.off()

res <- 400
png(file = file.path(DIR_FIG_DEA, "hs_visium_pseudobulk_data_grouped_by_donor_deseq2_res_spatial_selected_lowres.png"), 
    width = 15*res, height = 12*res, res = res)
(p1|p2)+plot_layout(widths = c(1,length(genes_plot_selected)))
dev.off()


