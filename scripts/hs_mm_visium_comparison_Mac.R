#' [hs_mm_visium_comparison_Mac.R]
#'
#' Comparison of macrophage populations in both IPF and BLM Visium data 
#'
#'
#' Nov 2023, L. Franzén [lovisa.franzen@scilifelab.se]

#### Set up ####
##### Define params. ####
SPECIES <- "mouse"
SPECIES <- "human"
DIR_ROOT <- getwd()
DIR_DATA <- file.path(DIR_ROOT, "data")
DIR_RES <- file.path(DIR_ROOT, "results", "translational")
DIR_FIG <- file.path(DIR_RES, "figures")
DIR_FIG_OUT <- file.path(DIR_FIG, "Mac")
DIR_OBJ_OUT <- file.path(DIR_RES, "objects", "Mac")
dir.create(DIR_FIG_OUT); dir.create(DIR_OBJ_OUT)


##### Load libs ####
library(tidyverse)
library(dplyr)
library(STutility)
library(patchwork)
library(readxl)
library(writexl)
library(pheatmap)

##### Other ####a
source(file.path(DIR_ROOT, "scripts", "custom_functions.R"))
source(file.path(DIR_ROOT, "scripts", "custom_colors.R"))
theme_custom <- theme(axis.title.x = element_blank())
fig_res <- 300

#### Gene conversion data ####
gene_conv_df <- read.csv(file.path(DIR_DATA, "misc", "orthogene_conv_combined_hs_mm.csv"), row.names = 1)
rownames(gene_conv_df) <- gene_conv_df$symbol_hs_mm
gene_conv <- setNames(gene_conv_df$symbol_mm, nm = gene_conv_df$symbol_hs)
gene_conv_mm <- setNames(gene_conv_df$symbol_hs, nm = gene_conv_df$symbol_mm)

#### Read gene annotation data ####
mm_cell_anno <- read.csv(file.path(DIR_ROOT, "data", "misc", "strunz_cell_type_groups.csv"), sep = ";")
hs_cell_anno <- read.csv(file.path(DIR_ROOT, "data", "misc", "habermann_cell_type_groups.csv"), sep = ";")

rownames(hs_cell_anno) <- hs_cell_anno$cell_name
rownames(mm_cell_anno) <- ifelse(!is.na(mm_cell_anno$cell_name), mm_cell_anno$cell_name, "NA")


#### Factor comparison ####
#' In the NMF analyses, we identified factors in the hs/mm that both overlapped
#' and corresponded to macrophages:
#' mmNMF30-d21 F8
#' hsNMF30 F5
#' 
#' Most likely, these are tissue resident macrophagges, seen in both healthy
#' fibrotic, however they are over-represented in treated/diseased samples
#'

##### NMF gene overlap #####
n_factors <- 30

#' Read tables generated with script 'hs_mm_visium_comparison_NMF.R'
hs_nmf_genes <- read.csv(file.path(DIR_ROOT, "results", "human", "objects/A_NMF30", "hs_visium_A_nmf_1-30_top100_gene_loadings.csv"))
mm_nmf_d21_genes <- read.csv(file.path(DIR_ROOT, "results", "mouse", "objects/NMF30_d21", "mm_visium_nmf30_d21_top100_gene_loadings.csv"))
mm_nmf_d7_genes <- read.csv(file.path(DIR_ROOT, "results", "mouse", "objects/NMF30_d7", "mm_visium_nmf30_d7_top100_gene_loadings.csv"))


#' Mac factor
hs_mac_fac_df <- subset(hs_nmf_genes, factor == 5)
mm_d21_mac_fac_df <- subset(mm_nmf_d21_genes, factor == 8)
mm_d7_mac_fac_df <- subset(mm_nmf_d7_genes, factor == 15)


hs_mac_fac_df$gene_mm <- gene_conv[hs_mac_fac_df$gene]
mm_d21_mac_fac_df$gene_hs <- gene_conv_mm[mm_d21_mac_fac_df$gene]
mm_d7_mac_fac_df$gene_hs <- gene_conv_mm[mm_d7_mac_fac_df$gene]

mm_d21_mac_fac_df$gene_upper <- toupper(mm_d21_mac_fac_df$gene)
mm_d7_mac_fac_df$gene_upper <- toupper(mm_d7_mac_fac_df$gene)

rownames(hs_mac_fac_df) <- hs_mac_fac_df$gene
rownames(mm_d21_mac_fac_df) <- mm_d21_mac_fac_df$gene_upper
rownames(mm_d7_mac_fac_df) <- mm_d7_mac_fac_df$gene_upper


hs_mac_fac_df$gene_loading_scaled <- scales::rescale(c(0, hs_mac_fac_df$gene_loading), to = c(0,1))[-1]
mm_d21_mac_fac_df$gene_loading_scaled <- scales::rescale(c(0, mm_d21_mac_fac_df$gene_loading), to = c(0,1))[-1]
mm_d7_mac_fac_df$gene_loading_scaled <- scales::rescale(c(0, mm_d7_mac_fac_df$gene_loading), to = c(0,1))[-1]


#' Shared genes
shared_genes <- c(intersect(hs_mac_fac_df$gene, mm_d21_mac_fac_df$gene_hs), 
                  intersect(hs_mac_fac_df$gene, mm_d7_mac_fac_df$gene_hs),
                  intersect(mm_d21_mac_fac_df$gene_hs, mm_d7_mac_fac_df$gene_hs)) %>% unique()

# shared_genes <- intersect(hs_mac_fac_df$gene, mm_d21_mac_fac_df$gene_hs)
shared_genes


shared_genes_upper <- c(intersect(hs_mac_fac_df$gene, mm_d21_mac_fac_df$gene_upper), 
                        intersect(hs_mac_fac_df$gene, mm_d7_mac_fac_df$gene_upper),
                        intersect(mm_d21_mac_fac_df$gene_upper, mm_d7_mac_fac_df$gene_upper)) %>% unique()
# shared_genes_upper <- intersect(hs_mac_fac_df$gene, mm_d21_mac_fac_df$gene_upper)

intersect(shared_genes, shared_genes_upper)


#' Plot heatmap
plot_hm_df <- data.frame(gene = shared_genes_upper,
                         hs_mac = hs_mac_fac_df[shared_genes_upper, "gene_loading_scaled"],
                         mm_d7_mac = mm_d7_mac_fac_df[shared_genes_upper, "gene_loading_scaled"],
                         mm_d21_mac = mm_d21_mac_fac_df[shared_genes_upper, "gene_loading_scaled"])
plot_hm_df_long <- pivot_longer(plot_hm_df, cols = c("hs_mac", "mm_d7_mac", "mm_d21_mac"), names_to = "group", values_to = "weight")
gene_order <- plot_hm_df_long %>% 
  filter(group=="hs_mac") %>% 
  arrange(desc(weight)) %>% 
  pull(gene)

plot_hm_df_long$gene <- factor(plot_hm_df_long$gene, levels = gene_order)
plot_hm_df_long$group <- factor(plot_hm_df_long$group, levels = c("hs_mac", "mm_d7_mac", "mm_d21_mac") %>% rev())
plot_hm_df_long$weight <- as.numeric(plot_hm_df_long$weight)

p <- ggplot(plot_hm_df_long, aes(x=gene, y=group, fill=weight)) +
  geom_tile() +
  scale_fill_gradientn(colours = col_scale_mako, na.value = "grey80", limits = c(0,1)) +
  theme_dotplot +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
        axis.text.y = element_text(face = "plain"),
        legend.position = "top");p

#' Export
pdf(file.path(DIR_FIG_OUT, "hs_mm_visium_M2Mac_factor_gene_heatmap_all_shared.pdf"), 
    width = 15, height = 2);p;dev.off()

write.csv(plot_hm_df, file.path(DIR_OBJ_OUT, "hs_mm_visium_M2Mac_factor_gene_heatmap_all_shared.csv"), row.names = F)



##### Mac-hi spot subset ######

###### Human ####### 
VlnPlot(se_hs, features = "factor_5", group.by = "condition")
f <- "factor_5"

#' Define cutoff (99th percentile) based on all data
se_nmf_emb <- se_hs@reductions$NMF@cell.embeddings
f_cutoff <- as.numeric(quantile(se_nmf_emb[, f], c(.99)))
f_metadata_add <- ifelse(se_nmf_emb[,f]>f_cutoff,"high", "low")

#' Add new metadata with cutoff info to data subsets
se_hs <- AddMetaData(se_hs, metadata = f_metadata_add, col.name = "mac_factor")
se_hs_subset <- AddMetaData(se_hs_subset, metadata = f_metadata_add, col.name = "mac_factor")

#' Check markers (should be similar to gene loadings)
# se_hs_subset <- SetIdent(se_hs_subset, value =  "mac_factor")
# hs_mac_markers <- FindMarkers(se_hs_subset, ident.1 = "high", ident.2 = "low", only.pos = T)

#' Plot spatial subset
se_hs_subset$mac_anno <- ifelse(se_hs_subset$mac_factor == "high", "Mac-F-hi", 
                                ifelse(se_hs_subset$annotation %in% c("Diseased", "Suspect Fibrosis"), "fibrosis", "other")
                                )
se_hs_subset$fib_anno <- ifelse(se_hs_subset$annotation %in% c("Diseased", "Suspect Fibrosis"), "fibrosis", "other")

ST.FeaturePlot(se_hs_subset, features = "mac_anno", ncol = 4, cols = c("grey60", "#81E97B", "grey90"), show.sb = F)


#' Subset data
se_hs_subset_mac <- SubsetSTData(se_hs_subset, mac_factor == "high")

se_hs_mac <- SubsetSTData(se_hs, mac_factor == "high")


message(paste("n spots in mac subset:", (se_hs_subset_mac@meta.data %>% nrow())))


###### Mouse ####### 
VlnPlot(se_mmd21, features = "factor_8", group.by = "condition")
f <- "factor_8"

#' d21: Define cutoff (99th percentile) based on all data
se_nmf_emb <- se_mmd21@reductions$NMF@cell.embeddings
f_cutoff <- as.numeric(quantile(se_nmf_emb[, f], c(.99)))
f_metadata_add <- ifelse(se_nmf_emb[,f]>f_cutoff,"high", "low")

#' Add new metadata with cutoff info to data subsets
se_mmd21 <- AddMetaData(se_mmd21, metadata = f_metadata_add, col.name = "mac_factor")
se_mmd21_subset <- AddMetaData(se_mmd21_subset, metadata = f_metadata_add, col.name = "mac_factor")

#' d7: Define cutoff (99th percentile) based on all data
f_d7 <- "factor_15"
se_nmf_emb <- se_mmd7@reductions$NMF@cell.embeddings
f_cutoff <- as.numeric(quantile(se_nmf_emb[, f], c(.99)))
f_metadata_add <- ifelse(se_nmf_emb[,f]>f_cutoff,"high", "low")

#' Add new metadata with cutoff info to data subsets
se_mmd7 <- AddMetaData(se_mmd7, metadata = f_metadata_add, col.name = "mac_factor")


#' d21: Check markers (should be similar to gene loadings)
# se_mmd21_subset <- SetIdent(se_mmd21_subset, value =  "mac_factor")
# hs_mac_markers <- FindMarkers(se_mmd21_subset, ident.1 = "high", ident.2 = "low", only.pos = T)

#' Plot spatial
se_mmd21_subset$mac_anno <- ifelse(se_mmd21_subset$mac_factor == "high", "Mac-F-hi", 
                                ifelse(se_mmd21_subset$annotation %in% c("Suspect Fibrosis/Fibroplasia"), "fibrosis", "other")
)
se_mmd21_subset$fib_anno <- ifelse(se_mmd21_subset$annotation %in% c("Suspect Fibrosis/Fibroplasia"), "fibrosis", "other")

anno_cols <- setNames(c("grey60", "#81E97B", "grey90"), nm = c("fibrosis", "Mac-F-hi", "other"))
ST.FeaturePlot(se_mmd21_subset, features = "mac_anno", ncol = 3, cols = anno_cols, show.sb = F) & theme(aspect.ratio = 1)

VlnPlot(se_mmd21_subset, group.by = "mac_factor", features = c("c2l_Recruited.macrophages", "c2l_Resolution.macrophages"))


#' Plot spatial all (annotated)
#' D21
se_mmd21$mac_anno <- ifelse(se_mmd21$mac_factor == "high", "Mac-F-hi", 
                           ifelse(se_mmd21$annotation %in% c("Suspect Fibrosis/Fibroplasia"), "fibrosis", "other")
                           )
se_mmd21_filt <- SubsetSTData(se_mmd21, annotation != "NA")

cols_mac <- setNames(c("#3F4C90", "grey70", "grey90"), nm = c("Mac-F-hi", "fibrosis", "other"))

pdf(file.path(DIR_FIG_OUT, "mm_visium_M2Mac_factor_hi_anno_spatial.pdf"), 
    width = 6, height = 6)
ST.FeaturePlot(se_mmd21_filt, features = "mac_anno", ncol = 3, cols = cols_mac, 
               show.sb = F, label.by = "sample_name", pt.size = 0.7) & 
  theme(aspect.ratio = 1)
dev.off()

pdf(file.path(DIR_FIG_OUT, "mm_visium_M2Mac_factor_spatial.pdf"), 
    width = 6, height = 6)
ST.FeaturePlot(se_mmd21_filt, features = f, ncol = 3, cols = c("grey90",col_scale_mako), 
               show.sb = F, label.by = "sample_name", pt.size = 0.7) & 
  theme(aspect.ratio = 1)
dev.off()

#' D7
se_mmd7$mac_anno <- ifelse(se_mmd7$mac_factor == "high", "Mac-F-hi", 
                            ifelse(se_mmd7$annotation %in% c("Suspect Fibrosis/Fibroplasia", "Inflammation"), "fibrosis", "other")
                           )
se_mmd7_filt <- SubsetSTData(se_mmd7, annotation != "NA")

cols_mac <- setNames(c("#3F4C90", "grey70", "grey90"), nm = c("Mac-F-hi", "fibrosis", "other"))

pdf(file.path(DIR_FIG_OUT, "mm_visium_d7_M2Mac_factor_hi_anno_spatial.pdf"), 
    width = 6, height = 6)
ST.FeaturePlot(se_mmd7_filt, features = "mac_anno", ncol = 3, cols = cols_mac, 
               show.sb = F, label.by = "sample_name", pt.size = 0.7) & 
  theme(aspect.ratio = 1)
dev.off()

pdf(file.path(DIR_FIG_OUT, "mm_visium_d7_M2Mac_factor_spatial.pdf"), 
    width = 6, height = 6)
ST.FeaturePlot(se_mmd7_filt, features = f_d7, ncol = 3, cols = c("grey90",col_scale_mako), 
               show.sb = F, label.by = "sample_name", pt.size = 0.7) & 
  theme(aspect.ratio = 1)
dev.off()


#' all data: dist cutoff
se_mmd21_subset$F14_C0_close <- ifelse(se_mmd21_subset$r_dist_F14_C0 < 300, "close", 
                                       ifelse(se_mmd21_subset$r_dist_F14_C0 < 1000, "far", "rest"))

VlnPlot(se_mmd21_subset, 
        features = c("factor_8", "factor_22", "c2l_Recruited.macrophages", "c2l_Resolution.macrophages"),
        pt.size = 0, group.by = "F14_C0_close", ncol = 4) #& scale_y_log10()


#' Subset only Mac-F high spots
se_mmd21_subset_mac <- SubsetSTData(se_mmd21_subset, mac_factor == "high")


#### Spatial Mac-factor activity ####
#' HE overlay for selected
hs_sample_select <- c("IPF_1.B3.1", "IPF_2.B3.2", "IPF_3.B2.2")
mm_d21_sample_select <- c("d21_bleo_3", "d21_bleo_4a", "d21_bleo_5b")

se_hs$sample_name2 <- paste0(se_hs$subject_alias, ".B", se_hs$B_tissue_selection, ".", se_hs$replicate)
se_hs_he <- SubsetSTData(se_hs, sample_name2 %in% hs_sample_select)

se_mm_he <- SubsetSTData(se_mmd21_filt, sample_name %in% mm_d21_sample_select)

#' Load images
se_hs_he <- LoadImages(se_hs_he, xdim = 2e3)
se_mm_he <- LoadImages(se_mm_he, xdim = 2e3)

#' Plot humamn
f_hs <- "factor_5"

pdf(file.path(DIR_FIG_OUT, "hs_visium_M2Mac_factor_spatial_HE.pdf"), 
    width = 8, height = 4)
FeatureOverlay(se_hs_he, features = f_hs, sampleids = 1:3, add.alpha = T, 
               ncols = 3, label.by = "sample_name2", 
               max.cutoff = 6, pt.size = 0.8,
               cols = rev(col_scale_mako), value.scale = "all") & 
  theme(aspect.ratio = 1, legend.position = "bottom")
dev.off()


#' Plot mouse d21
f_d21 <- "factor_8"

pdf(file.path(DIR_FIG_OUT, "mm_visium_M2Mac_factor_spatial_HE.pdf"), 
    width = 8, height = 4)
FeatureOverlay(se_mm_he, features = f_d21, sampleids = 1:3, add.alpha = T, 
               ncols = 3, label.by = "sample_name", 
               max.cutoff = 6, pt.size = 0.8,
               cols = rev(col_scale_mako), value.scale = "all") & 
  theme(aspect.ratio = 1, legend.position = "bottom")
dev.off()

