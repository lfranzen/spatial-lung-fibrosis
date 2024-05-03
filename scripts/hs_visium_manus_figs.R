#' [hs_visium_manus_figs.R]
#'
#' Generate figures for manuscript 
#' Human Visium data
#'
#' Oct 2022 - Nov 2023, L. Franzén [lovisa.franzen@scilifelab.se]


#### Set up ####
##### Define params. ####
set.seed(1)
SPECIES <- "human"
DIR_ROOT <- getwd()
DIR_DATA <- file.path(DIR_ROOT, "data", SPECIES, "visium")
DIR_RES <- file.path(DIR_ROOT, "results", SPECIES)
fig_res <- 500
DIR_FIG <- file.path(DIR_RES, "figures", "manus_figs")
dir.create(path = DIR_FIG)


##### Load libs ####
library(STutility)
library(RColorBrewer)
library(patchwork)
library(dplyr)
library(magrittr)
library(tidyr)
library(writexl)
library(readxl)

##### Other ####
source(file.path(DIR_ROOT, "scripts", "custom_colors.R"))
source(file.path(DIR_ROOT, "scripts", "custom_functions.R"))

theme_custom <- theme(axis.title.x = element_blank(), text = element_text(size=10, color="black"))
theme_transparent_bg <- theme(
  panel.background = element_rect(fill='transparent'),
  plot.background = element_rect(fill='transparent', color=NA),
  legend.background = element_rect(fill='transparent', color = NA),
  legend.box.background = element_rect(fill='transparent', color = NA)
)

theme_black_text <- theme(
  text = element_text(color="black"),
  plot.title = element_text(color="black", hjust=0.5, size=10),
  axis.title = element_text(color="black", hjust=0.5, size=10),
  axis.text = element_text(color="black", size=10),
  axis.ticks = element_line(color="black"), 
  axis.line = element_blank(), 
  panel.border = element_rect(color="black"),
  panel.grid = element_blank()
)

##### Read data ####
metadata <- read.table(file.path(DIR_DATA, "hs_visium_metadata_deposit.tsv"), sep = "\t", header = T)
metadata$fibrotic_extent_score_by_pathologist_0.3 <- as.character(metadata$fibrotic_extent_score_by_pathologist_0.3)
rownames(metadata) <- metadata$sample_id
metadata$n <- 1:nrow(metadata)

#'Workflow A
fname <- paste0("hs_visium_preproc_A_se_obj_nmf.rds")
se.subset <- readRDS(file = file.path(DIR_RES, "objects", fname))

se.subset$f14_nbs_clusters2 <- as.character(se.subset$f14_nbs_clusters)
rename_spots_clusters <- grep("^[0-9]$", se.subset$f14_nbs_clusters, value=T)
se.subset$f14_nbs_clusters2[names(rename_spots_clusters)] <- paste0("F14C0_nbs_C", rename_spots_clusters)
se.subset$f14_nbs_clusters2 <- factor(se.subset$f14_nbs_clusters2, 
                                      levels = c(paste0("F14C0_nbs_C", 0:5), "F14_C0", "other"))


#' A: NMF-F14-high
#' R object created in 'hs_visium_nmf.R' script
fname <- paste0("hs_visium_A_nmf_1-", 30, "_f14_subclusters_markers.csv")
makers_f14_subcluster <- read.csv(file = file.path(DIR_RES, "objects", "A_NMF30", fname))
se.f14high <- readRDS(file = file.path(DIR_RES, "objects", "hs_visium_preproc_A_se_obj_f14high_subset.rds"))

# nbs clusters (csv file created in 'hs_visium_nmf.R' script)
fname <- paste0("hs_visium_A_nmf_1-", 30, "_f14_subclusters_nbs_cluster_markers.csv")
makers_f14_subcluster_nbs <- read.csv(file = file.path(DIR_RES, "objects", "A_NMF30", fname))


#' Workflow B
fname <- paste0("hs_visium_preproc_B_se_obj_list.rds")
se_subset_split <- readRDS(file = file.path(DIR_RES, "objects", fname))
subject_names <- sort(names(se_subset_split))


#'Habermann cell2location data
c2l_all <- read.csv(file.path(DIR_ROOT, "data", SPECIES, "sc_deconvolution_habermann", "compiled_all_samples_cell_abundances.csv"), row.names = 1)
cell_type_col_names <- colnames(c2l_all)
cell_type_names <- gsub("c2l_", "", gsub("[.]$", "", gsub("..", "_", cell_type_col_names, fixed = TRUE)))

se.subset <- AddMetaData(se.subset, c2l_all)


habermann_cell_names <- read.csv(file.path(DIR_ROOT, "data", "misc", "habermann_cell_type_groups.csv"), sep = ";")
habermann_cell_names$c2l_col_names <- paste0("c2l_", habermann_cell_names$cell_name)


# saveRDS(se.subset[[]], file.path(DIR_RES, "objects", "hs_visium_preproc_A_se_full_metadata.rds"))


#### Supplementary figures ####

##### H&E image overview ####
se.subset <- LoadImages(se.subset, xdim = 1e3)

pdf(file = file.path(DIR_FIG, "hs_visium_fig_HE_all.pdf"), width = 16, height = 16)    # !Manuscript Extended Data figure 1a
ImagePlot(se.subset, annotate = T, ncols = 5)
dev.off()

metadata[, c("n", "sample_name","subject_alias", "condition")]

# mouse
se <- LoadImages(se, xdim = 1e3)
pdf(file = file.path(DIR_RES, "figures", "mm_visium_fig_HE_all.pdf"), width = 16, height = 16)
ImagePlot(se, annotate = T, ncols = 5)
dev.off()
metadata[, c("n", "sample_name","animal", "condition")]

##### Visium QC Plots ####
qc_stats <- read.csv(file.path(DIR_RES, "objects", "hs_visium_preproc_A_sample_summary.csv"), row.names = 1)
theme_qc <- theme(legend.position = "none", 
        panel.grid = element_blank(), 
        text = element_text(color="black", size=10), 
        plot.title = element_text(size=10, hjust=0.5))

p1 <- ggplot(qc_stats, aes(x=sample_name, y=n_spots, fill=subject_alias)) +
  geom_col(color="black", size=0.2) +
  scale_fill_manual(values = cols_donor) +
  labs(title="number of spots", x="", y="") +
  coord_flip() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 5e3)) +
  theme_linedraw() +
  theme_qc
  
p2 <- ggplot(qc_stats, aes(x=sample_name, y=avg_genes_spot, fill=subject_alias)) +
  geom_col(color="black", size=0.2) +
  scale_fill_manual(values = cols_donor) +
  labs(title="avg. unique genes / spot", x="", y="") +
  coord_flip() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3500)) +
  theme_linedraw() +
  theme_qc +
  theme(axis.text.y = element_blank())

p3 <- ggplot(qc_stats, aes(x=sample_name, y=avg_umis_spot, fill=subject_alias)) +
  geom_col(color="black", size=0.2) +
  scale_fill_manual(values = cols_donor) +
  labs(title="avg. UMIs / spot", x="", y="") +
  coord_flip() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 8500)) +
  theme_linedraw() +
  theme_qc +
  theme(axis.text.y = element_blank())
# p <- p1+p2+p3
# pdf(file = file.path(DIR_FIG, "hs_visium_fig_qc_stats.pdf"), width = 8, height = 4);p;dev.off()

p4 <- ggplot() +
  geom_boxplot(data = qc_stats, mapping = aes(x=condition, y=avg_genes_spot, fill=condition), 
               color="black", size=0.2, outlier.size = 0.5) +
  scale_fill_manual(values = cols_cond) +
  labs(x="", y="") +
  theme_linedraw() +
  theme_qc
p5 <- ggplot() +
  geom_boxplot(data = qc_stats, mapping = aes(x=condition, y=avg_umis_spot, fill=condition), 
               color="black", size=0.2, outlier.size = 0.5) +
  scale_fill_manual(values = cols_cond) +
  labs(x="", y="") +
  theme_linedraw() +
  theme_qc
# p <- p4+p5
# pdf(file = file.path(DIR_FIG, "hs_visium_fig_qc_stats_box.pdf"), width = 4, height = 2);p;dev.off()

panel1 <- (p1/plot_spacer() + patchwork::plot_layout(heights = c(3,1)))
panel2 <- (p2/p4 + patchwork::plot_layout(heights = c(3,1)))
panel3 <- (p3/p5 + patchwork::plot_layout(heights = c(3,1)))
pdf(file = file.path(DIR_FIG, "hs_visium_fig_qc_stats_panel.pdf"), width = 8, height = 6);panel1|panel2|panel3;dev.off()


##### Visium QC spatial ####
p1 <- ST.FeaturePlot(se.subset, features = "nFeature_RNA", ncol=5, pt.size = 0.4, 
                     cols = col_scale_rocket, show.sb = F, label.by = "sample_name") & 
  labs(fill="nFeature_RNA") &
  theme(plot.title = element_text(hjust=0.5, size=10), text = element_text(size=10))
p1 <- ggrastr::rasterize(p1, layers = "Point", dpi = fig_res)
p2 <- ST.FeaturePlot(se.subset, features = "nCount_RNA", ncol=5, pt.size = 0.4,
                     # max.cutoff = 35e3,
                     cols = col_scale_rocket, show.sb = F, label.by = "sample_name") & 
  labs(fill="nCount_RNA") &
  theme(plot.title = element_text(hjust=0.5, size=10), text = element_text(size=10));p2
p2 <- ggrastr::rasterize(p2, layers = "Point", dpi = fig_res)
pdf(file = file.path(DIR_FIG, "hs_visium_fig_qc_stats_spatial.pdf"), width = 8.5, height = 8)
p1;p2
dev.off()


##### Visium QC spatial gene expression ####
PlotViolinGene <- function (seurat.object, gene) {
  p1 <- VlnPlot(seurat.object, features = gene, group.by = "subject_alias", pt.size = 0, cols = cols_donor) & theme_custom & NoLegend()
  p2 <- VlnPlot(seurat.object, features = gene, group.by = "fibrotic_extent_score_by_pathologist_0.3", pt.size = 0, cols = cols_grade) & theme_custom & NoLegend()
  p <- p1/p2
  return(p)
}
p1 <- PlotViolinGene(se.subset, "GAPDH")
p2 <- PlotViolinGene(se.subset, "COL1A1")
p3 <- ST.FeaturePlot(se.subset, features = "COL1A1", ncol=5, pt.size = 0.4, slot = "data",
                     cols = col_scale_rocket, show.sb = F, label.by = "sample_name") & 
  theme(plot.title = element_text(hjust=0.5, size=10), text = element_text(size=10))
p3 <- ggrastr::rasterize(p3, layers = "Point", dpi = fig_res)

pdf(file = file.path(DIR_FIG, "hs_visium_fig_qc_geneexpr_spatial.pdf"), width = 13, height = 8)
p1|p2|p3|patchwork::plot_layout(widths = c(1,1,4))
dev.off()


#### Figure 1 ####
##### 1C: Stats ##### 
sample_stats <- read.csv(file = file.path(DIR_RES, "objects", "hs_visium_preproc_A_sample_summary.csv"))
sample_stats$grade <- strsplit(x = sample_stats$sample_name, split = "\\.") %>% unlist() %>% grep(pattern = "^G[0-9]", value = T) %>% gsub(pattern = "G", replacement = "")
sample_stats$rep <- strsplit(x = sample_stats$sample_name, split = "\\.") %>% unlist() %>% grep(pattern = "^[0-9]", value = T)
sample_stats$subject_alias <- factor(sample_stats$subject_alias, levels = c(paste0("IPF_", 4:1), paste0("HC_", 4:1)))

p1 <- ggplot(sample_stats, aes(x=subject_alias, y=n_spots, fill=grade)) +
  geom_col() +
  scale_fill_manual(values = cols_grade) +
  coord_flip() +
  labs(title="Spots", x="Donor") +
  theme_bw() +
  theme_transparent_bg +
  theme_black_text +
  theme(legend.position = "bottom", 
        axis.title.x = element_blank())

p1.2 <- ggplot(sample_stats, aes(x=subject_alias, y=n_spots, fill=grade, color = sample_name)) +
  geom_col(color="black") +
  scale_fill_manual(values = cols_grade) +
  coord_flip() +
  labs(title="Spots", x="Donor") +
  theme_bw() +
  theme_transparent_bg +
  theme_black_text +
  theme(legend.position = "bottom",
        axis.title.x = element_blank())

p2 <- ggplot(sample_stats, aes(x=subject_alias, y=avg_genes_spot, fill=condition)) +
  geom_boxplot(color="black") +
  geom_point(mapping = aes(fill=grade), size=2, shape=21, color="black") +
  scale_color_manual(values = cols_grade) +
  scale_fill_manual(values = c(cols_cond, cols_grade)) +
  coord_flip() +
  labs(title="Unique genes (avg.)", x="Donor") +
  theme_bw() +
  theme_transparent_bg +
  theme_black_text +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.title.x = element_blank())

pdf(file = file.path(DIR_FIG, "hs_visium_fig1_sample_stats.pdf"), width = 4.25, height = 2.75)
p1|p2
p1.2|p2
dev.off()


##### 1D: DEA ##### 
#' See script 'hs_visium_pseudobulk.R'

##### 1E: A-NMF30 results in one tissue #####
se.he <- SubsetSTData(se.subset, expression = sample_name %in% "IPF_1.TD012.B2.1")
se.he <- LoadImages(se.he, xdim = 2e3)
# ImagePlot(se.he, annotate = F)
# FeatureOverlay(se.he, features = "factor_6", cols = rev(col_scale_mako), add.alpha = T, pt.size = 2) & theme(legend.position = "bottom", plot.subtitle = element_blank())

se.he_crop <- se.he
se.he_crop <- CropImages(se.he_crop, crop.geometry.list = list("1" = "650x650+1000+350"))
ImagePlot(se.he_crop, annotate = F)

factors_plot <- c(
  1, # Ciliated cells
  # 21,  # Goblet
  10,  # Smooth muscle
  6,  # B/Plasma cells
  3 # BV
)

p_list <- lapply(factors_plot, function(f){
  p <- FeatureOverlay(se.he_crop, features = paste0("factor_",f), 
                 cols = rev(col_scale_mako), max.cutoff = 5,
                 add.alpha = T, pt.size = 2) & 
    theme(legend.position = "top", plot.subtitle = element_blank())
  p <- ggrastr::rasterize(p, layers = "Point", dpi = fig_res)
})
p1 <- patchwork::wrap_plots(p_list, nrow = 1) & theme(plot.margin = margin(0, 5, 0, 0))

p_list2 <- lapply(factors_plot, function(f){
  factor_gene_loadings <- read_xlsx(file.path(DIR_RES, "objects", "A_NMF30", paste0("hs_visium_A_nmf_1-30_top100_gene_loadings.xlsx")),sheet = paste0("Sheet",f))
  d_plot <- subset(factor_gene_loadings, rank <= 15)
  ggplot(d_plot, aes(x=reorder(gene, rank), y=gene_loading_scaled)) +
    geom_segment(aes(x = reorder(gene, rank), xend = reorder(gene, rank), y = 0, yend = gene_loading_scaled)) +
    geom_point(color="black", size=1) +
    labs(x="", y="") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.1)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle=45, hjust=1, color="black", size=8, face = "italic"), 
          axis.ticks.x = element_line(color = "black"),
          axis.line.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.text.y = element_blank(),
          panel.grid = element_blank())
})
p2 <- patchwork::wrap_plots(p_list2, nrow = 1) & theme(plot.margin = margin(0, 5, 0, 0))

pdf(file = file.path(DIR_FIG, "hs_visium_fig_A-NMF30_selected_crop.pdf"), width = 10, height = 4.5)
p1/p2 + patchwork::plot_layout(heights = c(4,1))
dev.off()

# se.he_crop2 <- se.he
# se.he_crop2 <- CropImages(se.he_crop2, crop.geometry.list = list("1" = "650x650+600+400"))
# FeatureOverlay(se.he_crop2, features = "factor_14", cols = rev(col_scale_mako), add.alpha = T, pt.size = 2.5) & theme(legend.position = "bottom", plot.subtitle = element_blank())



##### 1F: Correlation matrix Factor-Cell #### 
factor_cell_cor_data <- read.csv(file.path(DIR_RES, "objects", "A_NMF30", "hs_visium_A_nmf_1-30_factorweight_cell2location_habermann_metadata.csv"), 
                                 row.names = 1)

names_factors <- grep("factor_", colnames(factor_cell_cor_data), value = T)
names_cells <- colnames(factor_cell_cor_data)[!colnames(factor_cell_cor_data) %in% names_factors][-c(1,2,3)]

factor_cell_cor_data_filt <- factor_cell_cor_data %>% mutate_at(all_of(names_cells), ~if_else(.<0.25, 0, .)) # filter low values

factor_cell_cor_mat <- factor_cell_cor_data_filt %>% 
  select(-c(sample_name, subject_alias, condition)) %>% 
  cor()
diag(factor_cell_cor_mat) <- NA

names_factors <- grep("factor_", colnames(factor_cell_cor_mat), value = T)
names_cells <- colnames(factor_cell_cor_mat)[!colnames(factor_cell_cor_mat) %in% names_factors]
factor_cell_cor_mat <- factor_cell_cor_mat[names_cells, c(names_factors, names_cells)]

plot_groups <- data.frame(row.names = c(names_factors, names_cells),
                          group = c(
                            rep("NMF_factor", length(names_factors)),
                            rep("Cell_type", length(names_cells))),
                          cell_group = "NA")
plot_groups[names_cells, "cell_group"] <- c(
  "Epithelium",
  "Epithelium",
  "Lymphoid",
  "Epithelium",
  "Epithelium",
  "Epithelium",
  "Endothelium",
  "Mesenchyme",
  "Mesenchyme",
  "Epithelium",
  "Endothelium",
  "Epithelium",
  "Epithelium",
  "Myeloid",
  "Myeloid",
  "Mesenchyme",
  "Myeloid",
  "Mesenchyme",
  "Lymphoid",
  "Mesenchyme",
  "Lymphoid",
  "Epithelium",
  "Myeloid",
  "Lymphoid",
  "Epithelium",
  "Epithelium",
  "Mesenchyme",
  "Lymphoid",
  "Epithelium",
  "Myeloid",
  "Myeloid"
)
# plot_groups <- plot_groups %>% select(c(group, cell_group))
group_colors = list(
  group = c(NMF_factor = "#366A9FFF" , Cell_type = "#B26694"), # cols_donor[2], cols_donor[7]
  cell_group = setNames(c(
    # scico::scico(n = (plot_groups$cell_group %>% unique() %>% length() -1), palette = "roma"),
    RColorBrewer::brewer.pal(n = (plot_groups$cell_group %>% unique() %>% length() -1), name = "Pastel2"),
    "grey90"), 
    nm = plot_groups$cell_group %>% unique() %>% sort())
  )
pal_length <- length(col_scale_div_custom2)

pdf(file = file.path(DIR_FIG, paste0("hs_visium_fig_A-NMF30_Fig1_correlation_cell_cell_factor_heamap.pdf")),
    width = 14, height = 8, useDingbats = F)
abs_max_val <- max(factor_cell_cor_mat, na.rm = T)
ph_breaks <- c(seq(-abs_max_val, 0, length.out=ceiling(pal_length/2) + 1), 
               seq(abs_max_val/pal_length, abs_max_val, length.out=floor(pal_length/2)))
pheatmap::pheatmap(factor_cell_cor_mat,
                   cellwidth = 10, 
                   cellheight = 10, 
                   treeheight_row = 12,  
                   treeheight_col = 10,
                   annotation_row = plot_groups %>% select(cell_group), 
                   annotation_col = plot_groups %>% select(group), 
                   annotation_colors = group_colors,
                   color = col_scale_div_custom2,
                   na_col = "grey80",
                   cutree_cols = 3,
                   cutree_rows = 3,
                   breaks = ph_breaks
                   )
dev.off()

pdf(file = file.path(DIR_FIG, paste0("hs_visium_fig_A-NMF30_Fig1_correlation_cell_factor_heamap.pdf")),
    width = 8, height = 8, useDingbats = F)
factor_cell_cor_mat2 <- factor_cell_cor_mat[names_cells, names_factors]
abs_max_val <- max(factor_cell_cor_mat2, na.rm = T)
ph_breaks <- c(seq(-abs_max_val, 0, length.out=ceiling(pal_length/2) + 1), 
               seq(abs_max_val/pal_length, abs_max_val, length.out=floor(pal_length/2)))
pheatmap::pheatmap(factor_cell_cor_mat2,
                   cellwidth = 10, 
                   cellheight = 10, 
                   treeheight_row = 12,  
                   treeheight_col = 10,
                   annotation_row = plot_groups %>% select(cell_group), 
                   # annotation_col = plot_groups %>% select(group), 
                   annotation_colors = group_colors,
                   color = col_scale_div_custom2,
                   na_col = "grey80",
                   cutree_cols = 3,
                   cutree_rows = 3,
                   breaks = ph_breaks
                   )
dev.off()

pdf(file = file.path(DIR_FIG, paste0("hs_visium_fig_A-NMF30_Fig1_correlation_cell_cell_heamap.pdf")),
    width = 8, height = 8, useDingbats = F)
factor_cell_cor_mat3 <- factor_cell_cor_mat[names_cells, names_cells]
abs_max_val <- max(factor_cell_cor_mat3, na.rm = T)
ph_breaks <- c(seq(-abs_max_val, 0, length.out=ceiling(pal_length/2) + 1), 
               seq(abs_max_val/pal_length, abs_max_val, length.out=floor(pal_length/2)))
pheatmap::pheatmap(factor_cell_cor_mat3,
                   cellwidth = 10, 
                   cellheight = 10, 
                   treeheight_row = 12,  
                   treeheight_col = 10,
                   annotation_row = plot_groups %>% select(cell_group), 
                   # annotation_col = plot_groups %>% select(group), 
                   annotation_colors = group_colors,
                   color = col_scale_div_custom2,
                   na_col = "grey80",
                   cutree_cols = 3,
                   cutree_rows = 3,
                   breaks = ph_breaks
                   )
dev.off()


##### 1G: Factor/Cell/Gene spatial mapping #### 
# Demonstrate spatial overlap of factor/cells to further establish co-localization patterns
se.he <- SubsetSTData(se.subset, expression = sample_name %in% "IPF_1.TD012.B2.1")
se.he <- SubsetSTData(se.subset, expression = sample_name %in% c("HC_1.TH010.B0.1", "IPF_1.TD012.B2.1"))

p0 <- ST.FeaturePlot(se.he, features = "annotation", show.sb = F, label.by = "sample_name", cols = cols_annotation2, ncol = 2) & theme(legend.position = "bottom", aspect.ratio = 1)

p1 <- ST.FeaturePlot(se.he, features = "factor_12", show.sb = F, label.by = "sample_name", cols = c("grey95",col_scale_mako), min.cutoff = 0, max.cutoff = 2, ncol = 2)
p2 <- ST.FeaturePlot(se.he, features = "c2l_AT1", show.sb = F, label.by = "sample_name", cols = c("grey95",col_scale_rocket), min.cutoff = 0, max.cutoff = 3, ncol = 2)
p_top <- (p1|p2) & theme(legend.position = "bottom", aspect.ratio = 1)

p3 <- ST.FeaturePlot(se.he, features = "factor_9", show.sb = F, label.by = "sample_name", cols = c("grey95",col_scale_mako), min.cutoff = 0.5, max.cutoff = 3, ncol = 2)
p4 <- ST.FeaturePlot(se.he, features = "c2l_Endothelial.Cells", show.sb = F, label.by = "sample_name", cols = c("grey95",col_scale_rocket), min.cutoff = 0.5, max.cutoff = 5, ncol = 2)
p_mid <- (p3|p4) & theme(legend.position = "bottom", aspect.ratio = 1)

p5 <- ST.FeaturePlot(se.he, features = "factor_5", show.sb = F, label.by = "sample_name", cols = c("grey95",col_scale_mako), min.cutoff = 0.5, max.cutoff = 5, ncol = 2)
p6 <- ST.FeaturePlot(se.he, features = "c2l_Macrophages", show.sb = F, label.by = "sample_name", cols = c("grey95",col_scale_rocket), min.cutoff = 0.5, max.cutoff = 8, ncol = 2)
p_mid2 <- (p5|p6) & theme(legend.position = "bottom", aspect.ratio = 1)

p7 <- ST.FeaturePlot(se.he, features = "factor_4", show.sb = F, label.by = "sample_name", cols = c("grey95",col_scale_mako), min.cutoff = 0.5, max.cutoff = 5, ncol = 2)
p8 <- ST.FeaturePlot(se.he, features = "c2l_Fibroblasts", show.sb = F, label.by = "sample_name", cols = c("grey95",col_scale_rocket), min.cutoff = 0.5, max.cutoff = 8, ncol = 2)
p_bottom <- (p7|p8) & theme(legend.position = "bottom", aspect.ratio = 1)

pdf(file = file.path(DIR_FIG, paste0("hs_visium_fig_A-NMF30_Fig1_selected_factor_cell_spatial.pdf")),
    width = 14, height = 18, useDingbats = F)
(p0|plot_spacer())/p_top/p_mid/p_mid2/p_bottom
dev.off()


#### Figure 2 #### 

##### 2A-i: NMF gene heatmap ##### 
top_n_genes <- 5

all_factor_top_genes <- apply(se.subset@reductions$NMF@feature.loadings, 2, Scale01) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  pivot_longer(cols = paste0("factor_", 1:30), names_to = "factor_x", values_to = "gene_loading_scaled") %>% 
  mutate(factor_x = factor(factor_x, levels = paste0("factor_", 1:30)),
         factor_n = factor(gsub("factor_", "", factor_x), levels = 1:30)) %>% 
  arrange(factor_x, desc(gene_loading_scaled))


factor_select <- c(1:6,8:10,14,21)
all_factor_top_genes_filt <- all_factor_top_genes %>% filter(factor_n %in% factor_select)

genes_plot <- all_factor_top_genes_filt %>% 
  group_by(factor_x) %>% 
  dplyr::slice_max(order_by = gene_loading_scaled, n = top_n_genes) %>% 
  ungroup() %>% 
  select(gene) %>% 
  unique()

d_plot <- all_factor_top_genes_filt %>% filter(gene %in% genes_plot$gene)
d_plot$gene_order <- 1:nrow(d_plot)
d_plot$gene <- factor(d_plot$gene, levels = genes_plot$gene)


p_gw <- ggplot(d_plot, aes(y=gene, x=reorder(factor_n, desc(factor_n)), fill=gene_loading_scaled)) + 
  geom_tile(color=NA, width=1, height=1) +
  scale_fill_viridis(option = "G", direction = -1) +
  labs(x="Factor") +
  coord_flip() +
  theme_dotplot +
  theme(
    axis.text.y = element_text(face = "plain", size=8),
    axis.title.y = element_text(face = "plain", size=8),
    axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, face="plain", size=8),
    panel.grid = element_blank(),
    legend.position = "top");p_gw

pdf(file = file.path(DIR_FIG, paste0("hs_visium_fig_A-NMF30_Fig1_factor_gene_weight_heamap.pdf")),
    width = 8, height = 3, useDingbats = F)
p_gw
dev.off()


##### 2A-ii: NMF HC/IPF proportions ##### 
se_nmf_emb <- se.subset@reductions$NMF@cell.embeddings
factor_xy <- c(paste0("factor_", 1:30))

factor_cutoff <- lapply(factor_xy, function(x){
  as.numeric(quantile(se_nmf_emb[,x], c(.99)))
}) %>% setNames(nm = factor_xy)

f_metadata <- se.subset@meta.data %>% select(c(sample_name, B_tissue_selection)) # conditon
f_metadata$biopsy <- paste0("B", f_metadata$B_tissue_selection)

for(f in factor_xy){
  f_metadata[[f]] <- "low"
  spots_f_cutoff <- intersect(colnames(se.subset), rownames(se_nmf_emb[se_nmf_emb[,f]>factor_cutoff[[f]],]))
  f_metadata[spots_f_cutoff,f] <- "high"
}
f_metadata$bc <- rownames(f_metadata)


f_metadata_long <- pivot_longer(f_metadata, cols = all_of(factor_xy), names_to = "factor_n", values_to = "cutoff")

f_summary <- f_metadata_long %>% 
  group_by(biopsy, factor_n, cutoff) %>% 
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>% 
  filter(cutoff == "high")

f_summary$factor_n <- factor(f_summary$factor_n, levels=factor_xy %>% rev())
f_summary$biopsy <- factor(f_summary$biopsy, levels=unique(f_summary$biopsy) %>% rev())

factors_plot <- c(paste0("factor_", factor_select)) # selected factors

cols_biopsy <- cols_grade
names(cols_biopsy) <- paste0("B", names(cols_grade))
f_summary$factorF <- gsub("factor_", "F", f_summary$factor_n)
f_summary$factorF <- factor(f_summary$factorF, levels = paste0("F", 1:30) %>% rev())

p_fact_prop <- ggplot(
  data = subset(f_summary, factor_n %in% factors_plot),
  mapping = aes(x=factorF, y = freq, fill=biopsy)
) +
  geom_bar(stat="identity", width=0.85, position = "fill") +
  scale_fill_manual(values = cols_biopsy) +
  scale_x_discrete(position = "top") +
  labs(x="Factor") +
  coord_flip() +
  theme_dotplot +
  theme(
    axis.text.y = element_text(face = "plain", size=8),
    axis.title.y = element_text(face = "plain", size=8),
    axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, face="plain", size=8),
    panel.grid = element_blank(),
    legend.position = "right");p_fact_prop


# pdf(file = file.path(DIR_FIG, paste0("hs_visium_fig_A-NMF30_ExtDataFig2_factor_biopsy_prop_all.pdf")),
#     width = 3, height = 6, useDingbats = F);p_fact_prop;dev.off()
# write.csv(f_summary, file.path(DIR_FIG, "hs_visium_fig_A-NMF30_ExtDataFig2_factor_biopsy_prop_all_data.csv"), row.names = F)

pdf(file = file.path(DIR_FIG, paste0("hs_visium_fig_A-NMF30_Fig1_factor_biopsy_prop.pdf")),
    width = 3, height = 3, useDingbats = F)
p_fact_prop
p_fact_prop + geom_bar(stat="identity", width=0.85, position = "fill", color="black") + theme_void()
dev.off()

pdf(file = file.path(DIR_FIG, paste0("hs_visium_fig_A-NMF30_Fig2a_factor_gene_weight_heamap.pdf")),
    width = 10, height = 3.5, useDingbats = F)
(p_gw|p_fact_prop) + plot_layout(widths = c(10,1))
dev.off()



##### 2B: H&E panels #### 
factors_plot <- paste0("factor_", c(4, 14, 5, 21, 9))

# i-iii. H&E overlay
# Crop windows
se.ipf1_g1 <- SubsetSTData(se.subset, sample_name == "IPF_1.TD011.B1.1")
se.ipf1_g3 <- SubsetSTData(se.subset, sample_name == "IPF_1.TD013.B3.1")
se.ipf4_g3 <- SubsetSTData(se.subset, sample_name == "IPF_4.TD042.B3.1")


se.ipf1_g1 <- LoadImages(se.ipf1_g1, xdim = 2e3)
se.ipf1_g3 <- LoadImages(se.ipf1_g3, xdim = 2e3)
se.ipf4_g3 <- LoadImages(se.ipf4_g3, xdim = 2e3)

crop_size <- "1000x1000"

#' ipf1_g1
ImagePlot(se.ipf1_g1)
se.he_crop <- se.ipf1_g1
se.he_crop <- CropImages(se.he_crop, crop.geometry.list = list("1"=paste0(crop_size, "+750+800")))
# ImagePlot(se.he_crop, annotate = F)

p_ipf1_g1 <- FeatureOverlay(se.he_crop, 
               features = factors_plot[5], 
               label.by = "sample_name",
               cols = rev(col_scale_mako),
               # cols = c("grey95",col_scale_mako), 
               add.alpha = T,
               # min.cutoff = 0.5, 
               max.cutoff = 4.5,
               pt.alpha = 0.8,
               pt.size = 2,
               show.sb = T) &
  theme(aspect.ratio = 1,
        legend.position = "bottom", 
        # plot.title = element_blank(),
        plot.margin = margin(4, 5, 2.5, 5)); p_ipf1_g1


#' ipf1_g3
ImagePlot(se.ipf1_g3)
se.he_crop <- se.ipf1_g3
se.he_crop <- CropImages(se.he_crop, crop.geometry.list = list("1"=paste0(crop_size, "+700+200")))
# ImagePlot(se.he_crop, annotate = F)

p_ipf1_g3 <- FeatureOverlay(se.he_crop, 
                            features = factors_plot[3:4],
                            ncols=2,
                            label.by = "sample_name",
                            cols = rev(col_scale_mako),
                            # cols = c("grey95",col_scale_mako), 
                            add.alpha = T,
                            # min.cutoff = 0.5, 
                            max.cutoff = 6,
                            pt.alpha = 0.8,
                            pt.size = 2,
                            show.sb = T) &
  theme(aspect.ratio = 1,
        legend.position = "bottom", 
        # plot.title = element_blank(),
        plot.margin = margin(4, 5, 2.5, 5));p_ipf1_g3

#' ipf4_g3
ImagePlot(se.ipf4_g3)
se.he_crop <- se.ipf4_g3
se.he_crop <- CropImages(se.he_crop, crop.geometry.list = list("1"=paste0(crop_size, "+700+200")))
# ImagePlot(se.he_crop, annotate = F)

p_ipf4_g3 <- FeatureOverlay(se.he_crop, 
                            features = factors_plot[1:2],
                            ncols=2,
                            label.by = "sample_name",
                            cols = rev(col_scale_mako),
                            # cols = c("grey95",col_scale_mako), 
                            add.alpha = T,
                            # min.cutoff = 0.5, 
                            max.cutoff = 4,
                            pt.alpha = 0.8,
                            pt.size = 2,
                            show.sb = T) &
  theme(aspect.ratio = 1,
        legend.position = "bottom", 
        # plot.title = element_blank(),
        plot.margin = margin(4, 5, 2.5, 5));p_ipf4_g3


#' Export plots
pdf(file = file.path(DIR_FIG, paste0("hs_visium_fig_A-NMF30_Fig2b_factor_spatial.pdf")),
    width = 8, height = 4, useDingbats = F)
ST.FeaturePlot(se.ipf4_g3, features = factors_plot[1:2], cols = c("grey95",col_scale_mako), min.cutoff = 0.5, show.sb = F, label.by = "sample_name") & theme(legend.position = "bottom", aspect.ratio=1)
ST.FeaturePlot(se.ipf1_g3, features = factors_plot[3:4], cols = c("grey95",col_scale_mako), min.cutoff = 0.5, show.sb = F, label.by = "sample_name") & theme(legend.position = "bottom", aspect.ratio=1)
ST.FeaturePlot(se.ipf1_g1, features = factors_plot[5], cols = c("grey95",col_scale_mako), min.cutoff = 0.5, show.sb = F, label.by = "sample_name") & theme(legend.position = "bottom", aspect.ratio=1)
dev.off()

pdf(file = file.path(DIR_FIG, paste0("hs_visium_fig_A-NMF30_Fig2b_factor_spatial_HE_crop.pdf")),
    width = 8, height = 4, useDingbats = F)
p_ipf4_g3
p_ipf1_g3
p_ipf1_g1
dev.off()


# iv. Plot dotplot with factor activity per annotated area
se.subset.anno <- subset(se.subset, annotation != "NA")
se.subset.anno <- AddMetaData(se.subset.anno, se.subset@reductions$NMF@cell.embeddings, col.name = colnames(se.subset@reductions$NMF@cell.embeddings))

pdf(file = file.path(DIR_FIG, paste0("hs_visium_fig_A-NMF30_Fig2b_factor_anno_dotplot.pdf")),
    width = 3.5, height = 3.5, useDingbats = F)
DotPlot(se.subset.anno, features = rev(factors_plot), group.by = "annotation") + 
  scale_colour_gradient2(low = "#49C1ADFF", mid = "grey90", high = "#6D537F") + 
  coord_flip() +
  theme(axis.text.x = element_text(angle=45, hjust=1), 
        axis.title = element_blank(),
        axis.text = element_text(size=10),
        text = element_text(size=10),
        aspect.ratio = 1)
dev.off()


##### 2C: F14 activity rank #### 
#' A-NMF30-F14 Factor rank
rank_dat <- bind_cols(data.frame(factor_14 = se.subset@reductions$NMF@cell.embeddings[,"factor_14"]),
                      se.subset@meta.data)

rank_dat_summary <- rank_dat %>%
  dplyr::group_by(B_tissue_selection, sample_name, condition) %>%
  dplyr::slice_max(order_by = factor_14, n = 100) %>%
  summarise(factor_sum = sum(factor_14),
            factor_avg = mean(factor_14))

p_rank_box <- ggplot(rank_dat_summary, aes(x=condition, y=factor_sum, fill=as.factor(B_tissue_selection))) +
  geom_boxplot(width = 0.8) +
  labs(x="", y="sum", fill="") +
  scale_fill_manual(values = cols_grade) +
  theme_bw() +
  theme_black_text +
  theme(legend.position = "right",
        text = element_text(size=10, color="black"),
        axis.text = element_text(size=10, color="black"))

p_rank <- FactorRankSpotPlot(seurat.object = se.subset, 
                             factor = 14, 
                             top.n.spots = 100, 
                             include.zoom = F,
                             line.alpha = 0.8,
                             group.by = "B_tissue_selection",
                             color.group = cols_grade,
                             split.by = "sample_name" 
) &
  labs(x="spot rank") &
  theme_bw() &
  theme(
    plot.title = element_blank(),
    panel.grid = element_blank(),
    text = element_text(size=10, color="black"),
    axis.text = element_text(size=10, color="black"),
    axis.ticks = element_line(color = "black", size=0.5)
  ) &
  theme_transparent_bg &
  NoLegend();p_rank

pdf(file = file.path(DIR_FIG, "hs_visium_fig_A-NMF30_Fig2_F14_spot_actvity_rank.pdf"), width = 5, height = 2)
(p_rank | p_rank_box) + patchwork::plot_layout(widths = c(3,1))
dev.off()


##### 2D: F14 gene rank ####
#' Gene loadings (file produced in 'hs_visium_nmf.R' script)
factor_gene_loadings <- read_xlsx(file.path(DIR_RES, "objects", "A_NMF30", paste0("hs_visium_A_nmf_1-30_top100_gene_loadings.xlsx")),sheet = paste0("Sheet",14))

# Line with selected genes to highlight
factor_x <- 14
top_n_genes <- 100
feat_loads <- as.data.frame(se.subset@reductions$NMF@feature.loadings[,paste0("factor_",factor_x)])
colnames(feat_loads) <- "gene_loading"
feat_loads$gene <- rownames(feat_loads)
feat_loads$factor <- factor_x
feat_loads <- feat_loads %>%
  dplyr::slice_max(order_by = gene_loading, n = top_n_genes) %>%
  mutate(rank = dense_rank(-gene_loading))
feat_loads$rank <- as.factor(feat_loads$rank)
feat_loads$rank <- factor(feat_loads$rank, levels=(levels(feat_loads$rank)))
feat_loads$gene_loading <- as.numeric(feat_loads$gene_loading)

genes_label <- c("COL3A1", "SERPINH1", "CTHRC1", "POSTN", "TNC", "FN1", "KRT7", "KRT8")
genes_label <- feat_loads[genes_label,] %>% arrange(desc(rank)) %>% rownames()

p_gload <- ggplot(feat_loads, aes(x=rank, y=gene_loading, group=1)) +
  geom_point(stat='summary', fun=sum, size=NA) +
  stat_summary(fun=sum, geom="line", color = "grey10", linewidth=0.75) +
  geom_point(data = feat_loads[genes_label,], 
             shape=21,
             size=1.5, 
             color = "grey10",
             fill = "#F7931E"
  ) + 
  scale_x_discrete(breaks = c(as.character(feat_loads[genes_label,"rank"]), 100)) +
  labs(y="gene weight", x="genes ranked by F14 contribution") +
  theme_classic() +
  theme(panel.grid = element_blank(), 
        text = element_text(size=10, color="black"), 
        axis.text = element_text(size=10, color="black"), 
        axis.line = element_line(color = "black", size=0.5),
        axis.ticks = element_line(color = "black", size=0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
  ) +
  theme_transparent_bg;p_gload

p_txt <- ggplot(feat_loads) +
  annotate(geom = "text", 
           x=-10, y=-10, 
           label=paste(genes_label, collapse = "\n"),
           hjust = 0, size=2.5, color="black") +
  theme_void()

pdf(file = file.path(DIR_FIG, "hs_visium_fig_A-NMF30_Fig2_F14_gene_weight_rank.pdf"), width = 4, height = 2)
(p_gload|plot_spacer()|p_txt) + plot_layout(widths = c(1,-0.9,1))
dev.off()


##### 2E: F14 cell type correlation ####
#' Plotted in 'hs_visium_cell2location_proc_data.R' script (copied code below)
#' 'factor_cell_cor_data_group' object generated in 'hs_visium_cell2location_proc_data.R'

factor_cell_cor_data_group_f14 <- do.call(cbind,
                                          lapply(names(factor_cell_cor_data_group)[2:6], function(n){
                                            d <- t(factor_cell_cor_data_group[[n]])
                                            colnames(d) <- paste0(n, ".", colnames(d))
                                            return(d)
                                          })
) %>% 
  as.data.frame() %>% 
  select(contains("factor_14"))


pal_length <- length(col_scale_div_custom2)
ph_colors <- col_scale_div_custom2
abs_max_val <- max(factor_cell_cor_data_group_f14)
ph_breaks <- c(seq(-abs_max_val, 0, length.out=ceiling(pal_length/2) + 1), 
               seq(abs_max_val/pal_length, abs_max_val, length.out=floor(pal_length/2)))

pdf(file = file.path(DIR_FIG_C2L, paste0("hs_visium_c2l_res_cell_type_A-NMF30_factor_cor_heatmap_per_subject_F14.pdf")), 
    width = 6, height = 6, useDingbats = F)
pheatmap::pheatmap(factor_cell_cor_data_group_f14, 
                   cellwidth = 16, 
                   cellheight = 10, 
                   treeheight_row = 10,  
                   treeheight_col = 4,
                   color = ph_colors,
                   breaks = ph_breaks
)
dev.off()


##### 2F: H&E panels #### 
sample_select <- list(IPF_1 = "IPF_1.TD012.B2.1", 
                      IPF_3a = "IPF_3.TD031.B1.2", 
                      IPF_3b = "IPF_3.TD032B.B3.2", 
                      IPF_4 = "IPF_4.TD042.B3.1")


crop_size <- "1750x1750"
crop_window_list_full <- list("1" = paste0(crop_size, "+100+100"),  # IPF_1
                              "2" = paste0(crop_size, "+100+50"),  # IPF4
                              "3" = paste0(crop_size, "+110+200"), # IPF_3a
                              "4" = paste0(crop_size, "+150+100"))   # IPF_3b

crop_size <- "500x500"
crop_window_list_FF <- list("1" = paste0(crop_size, "+500+900"),  # IPF_1
                            "2" = paste0(crop_size, "+1100+950"),  # IPF4
                            "3" = paste0(crop_size, "+200+250"), # IPF_3a
                            "4" = paste0(crop_size, "+300+1000"))   # IPF_3b


se.he <- SubsetSTData(se.subset, expression = sample_name %in% unlist(sample_select))
se.he <- LoadImages(se.he, xdim = 2e3)

se.he_crop <- se.he
se.he_crop <- CropImages(se.he_crop, crop.geometry.list = crop_window_list_FF)
ImagePlot(se.he_crop, ncols = 4)


# plot Full size 
se.he <- CropImages(se.he, crop.geometry.list = crop_window_list_full)
ImagePlot(se.he, ncols = 4, annotate = F)

pdf(file = file.path(DIR_FIG, "hs_visium_fig_A-NMF30_Fig2_HE_full_HE.pdf"), width = 16, height = 4);ImagePlot(se.he, ncols = 4, annotate = F);dev.off()



###### A-NMF30-F14 #####
# HE
# f14_sample_select <- list(IPF_1 = "IPF_1.TD012.B2.1", 
#                           IPF_3a = "IPF_3.TD031.B1.2", 
#                           IPF_3b = "IPF_3.TD032B.B3.2", 
#                           IPF_4 = "IPF_4.TD042.B3.1")
# 
# se.he <- SubsetSTData(se.subset, expression = sample_name %in% unlist(f14_sample_select))
# se.he <- LoadImages(se.he, xdim = 2e3)
# 
# # Full size
# crop_size <- "1750x1750"
# crop_window_list <- list("1" = paste0(crop_size, "+100+100"),  # IPF_1
#                          "2" = paste0(crop_size, "+100+50"),  # IPF4
#                          "3" = paste0(crop_size, "+110+200"), # IPF_3a
#                          "4" = paste0(crop_size, "+150+100"))   # IPF_3b
# se.he_crop <- se.he
# se.he_crop <- CropImages(se.he_crop, crop.geometry.list = crop_window_list)



f14_max <- 3
p_f14 <- FeatureOverlay(se.he_crop, 
               sampleids = 1:4,
               ncols = 4,
               features = "factor_14", 
               label.by = "sample_name",
               cols = rev(col_scale_mako), 
               min.cutoff = 1,
               max.cutoff = f14_max,
               add.alpha = T, 
               pt.size = 3.5) &
  theme(legend.position = "bottom", 
        aspect.ratio = 1,
        plot.title = element_blank(),
        plot.margin = margin(4, 5, 2.5, 5));p_f14
p_f14 <- ggrastr::rasterize(p_f14, layers = "Point", dpi = fig_res)

pdf(file = file.path(DIR_FIG, "hs_visium_fig_A-NMF30_Fig2_HE_crop_F14.pdf"), width = 16, height = 4);p_f14;dev.off()

# 
# # Crop
# ImagePlot(se.he, indices = 1:4)
# crop_window_list <- list("1" = "700x700+180+1100",  # IPF_1
#                          "2" = "700x700+1000+850",  # IPF4
#                          "3" = "700x700+200+200", # IPF_3a
#                          "4" = "700x700+300+900")   # IPF_3b
# 
# se.he_crop <- se.he
# se.he_crop <- CropImages(se.he_crop, crop.geometry.list = crop_window_list)
# # ImagePlot(se.he_crop, annotate = F)
# p <- FeatureOverlay(se.he_crop, 
#                     sampleids = 1:4,
#                     ncols = 4,
#                     features = "factor_14", 
#                     label.by = "sample_name",
#                     cols = rev(col_scale_mako), 
#                     min.cutoff = 1,
#                     max.cutoff = f14_max,
#                     add.alpha = T, 
#                     pt.size = 2.5) &
#   theme(legend.position = "bottom", 
#         plot.title = element_blank(),
#         plot.margin = margin(4, 5, 2.5, 5))
# pdf(file = file.path(DIR_FIG, "hs_visium_fig_A-NMF30_F14_selected_crop.pdf"), width = 17, height = 5);p;dev.off()


###### Histopath #####
p_histo <- FeatureOverlay(se.he_crop, 
                        sampleids = 1:4,
                        ncols = 4,
                        features = "annotation", 
                        label.by = "sample_name",
                        cols = cols_annotation, 
                        pt.alpha = 0.8, 
                        pt.size = 3) &
  theme(legend.position = "bottom", 
        aspect.ratio = 1,
        plot.title = element_blank(),
        plot.margin = margin(4, 5, 2.5, 5));p_histo
p_histo <- ggrastr::rasterize(p_histo, layers = "Point", dpi = fig_res)

pdf(file = file.path(DIR_FIG, "hs_visium_fig_A-NMF30_Fig2_HE_crop_histopath.pdf"), width = 16, height = 4);p_histo;dev.off()


###### c2l_KRT5..KRT17. #####
p_abba <- FeatureOverlay(se.he_crop, 
                         sampleids = 1:4,
                         ncols = 4,
                         features = "c2l_KRT5..KRT17.", 
                         label.by = "sample_name",
                         cols = rev(col_scale_mako), 
                         add.alpha = T, 
                         pt.size = 3.5) &
  theme(legend.position = "bottom", 
        aspect.ratio = 1,
        plot.title = element_blank(),
        plot.margin = margin(4, 5, 2.5, 5));p_abba
p_abba <- ggrastr::rasterize(p_abba, layers = "Point", dpi = fig_res)

pdf(file = file.path(DIR_FIG, "hs_visium_fig_A-NMF30_Fig2_HE_crop_c2lKRT5KRT17.pdf"), width = 16, height = 4);p_abba;dev.off()



###### c2l_Myofibroblasts #####
p_mf <- FeatureOverlay(se.he_crop, 
                         sampleids = 1:4,
                         ncols = 4,
                         features = "c2l_Myofibroblasts", 
                         label.by = "sample_name",
                         cols = rev(col_scale_mako), 
                         add.alpha = T, 
                         pt.size = 3.5) &
  theme(legend.position = "bottom", 
        aspect.ratio = 1,
        plot.title = element_blank(),
        plot.margin = margin(4, 5, 2.5, 5));p_mf
p_mf <- ggrastr::rasterize(p_mf, layers = "Point", dpi = fig_res)
pdf(file = file.path(DIR_FIG, "hs_visium_fig_A-NMF30_Fig2_HE_crop_c2lMyofibroblasts.pdf"), width = 16, height = 4);p_mf;dev.off()


###### c2l_HAS1.High.fbroblasts #####
grep("HAS1", colnames(se.he_crop[[]]), value = T)
p_hf <- FeatureOverlay(se.he_crop, 
                       sampleids = 1:4,
                       ncols = 4,
                       features = "c2l_HAS1.High.Fibroblasts", 
                       label.by = "sample_name",
                       cols = rev(col_scale_mako), 
                       add.alpha = T, 
                       pt.size = 3.5) &
  theme(legend.position = "bottom", 
        aspect.ratio = 1,
        plot.title = element_blank(),
        plot.margin = margin(4, 5, 2.5, 5));p_hf
p_hf <- ggrastr::rasterize(p_hf, layers = "Point", dpi = fig_res)

pdf(file = file.path(DIR_FIG, "hs_visium_fig_A-NMF30_Fig2_HE_crop_c2lHAS1HiFibroblasts.pdf"), width = 16, height = 4);p_hf;dev.off()


#clean-up
rm(se.he_crop)
rm(se.he)


#### Fig 3 #### 
##### 3A: A-NMF30-F14-high subcluster UMAP ####
f14_high_spots <- colnames(se.f14high)
c2l_subset <- c2l_all[f14_high_spots, ]
se.f14high <- AddMetaData(se.f14high, c2l_subset)

# F14hi cutoff
se_nmf_emb <- se.subset@reductions$NMF@cell.embeddings %>% as.data.frame()
f14_cutoff <- as.numeric(quantile(se_nmf_emb[,"factor_14"], c(.99)))

p_f14_box <- se_nmf_emb %>% 
  select(factor_14) %>% 
  rownames_to_column() %>% 
  mutate(cutoff = ifelse(factor_14>f14_cutoff, "high", "low"),
         factor_x = "F14")

p_f14_bp <- ggplot(p_f14_box, aes(x=factor_x, y=factor_14)) +
  geom_boxplot(outlier.size = 0.1) +
  geom_hline(yintercept = f14_cutoff, color = "#4B3B66", alpha = 0.6) +
  labs(y="factor 14 activity") +
  # scale_y_log10() +
  theme_dotplot +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(face = "plain"),
        axis.title.y = element_text());p_f14_bp

pdf(file = file.path(DIR_FIG, "hs_visium_fig_A-NMF30_F14-high_selection_boxplot.pdf"), width = 1, height = 2);p_f14_bp;dev.off()

# Umap 
se.f14high$f14_subclusters2 <- paste0("C", se.f14high$f14_subclusters) %>% as.factor()
n_clusters <- 5
col_clusters <- setNames(c(col_scale_acton[8], RColorBrewer::brewer.pal(5, "Pastel1")), nm = levels(se.f14high$f14_subclusters))
col_clusters2 <- setNames(c(col_scale_acton[8], RColorBrewer::brewer.pal(5, "Pastel1")), nm = levels(se.f14high$f14_subclusters2))

p_umap <- DimPlot(se.f14high, 
                  reduction = "umap", 
                  group.by = "f14_subclusters", 
                  cols = col_clusters, 
                  pt.size = 0.5, #0.1
                  label = T,
                  label.box = T, 
                  repel = F) & 
  theme_void() +
  theme(panel.background = element_blank(), 
        axis.text = element_blank(), 
        panel.grid = element_blank(), 
        aspect.ratio = 1, 
        legend.position = "none", 
        plot.title = element_blank());p_umap
p_umap <- ggrastr::rasterize(p_umap, layers = "Point", dpi = 600)
pdf(file = file.path(DIR_FIG, "hs_visium_fig_A-NMF30_F14-high_subcluster_umap.pdf"), width = 2, height = 2);p_umap;dev.off()


#' With markers
# se.f14high$seurat_clusters2 <- paste0("C", as.numeric(se.f14high$seurat_clusters)) %>% as.factor()
# col_clusters2 <- setNames(c(scico::scico(10, palette = "acton")[3], RColorBrewer::brewer.pal(5, "Pastel1")), nm = levels(se.f14high$seurat_clusters2))

p_umap_markers <- FeaturePlot(
  se.f14high, 
  features = c("PRSS2", "KRT17", "KRT5", 
               grep("KRT5", colnames(c2l_all), value = T), "c2l_Myofibroblasts", "c2l_Basal"),  # c2l_Myofibroblasts, c2l_Fibroblasts, c2l_AT2, c2l_Macrophages
  reduction = "umap", 
  ncol = 6,
  cols = c(col_scale_acton, "black"),
  pt.size = 0.5
) & 
  theme_void() +
  theme(panel.background = element_blank(), 
        axis.text = element_blank(), 
        panel.grid = element_blank(), 
        aspect.ratio = 1, 
        plot.title = element_text(hjust = 0.5, size=10),
        legend.position = "bottom" 
  );p_umap_markers
p_umap_markers <- ggrastr::rasterize(p_umap_markers, layers = "Point", dpi = 600)

pdf(file = file.path(DIR_FIG, "hs_visium_fig_A-NMF30_F14-high_subcluster_umap_markers.pdf"), width = 6, height = 3);
(p_umap | patchwork::plot_spacer() | p_umap_markers) + patchwork::plot_layout(widths = c(1,0.1,1.2))
dev.off()

pdf(file = file.path(DIR_FIG, "hs_visium_fig_A-NMF30_F14-high_subcluster_umap_markers2.pdf"), width = 12, height = 3);
p_umap_markers
dev.off()



##### 3B: A-NMF30-F14-hi cluster spatial ####
# HE selected crops
f14_sample_select <- list(IPF_1 = "IPF_1.TD012.B2.1", 
                          IPF_3a = "IPF_3.TD031.B1.2", 
                          IPF_3b = "IPF_3.TD032B.B3.2", 
                          IPF_4 = "IPF_4.TD042.B3.1")

se.he <- SubsetSTData(se.f14high, expression = sample_name %in% unlist(f14_sample_select))
se.he <- LoadImages(se.he, xdim = 2e3)

# cropped
crop_size <- "500x500"
crop_window_list <- list(
  # "1" = paste0(crop_size, "+180+1200"),  # IPF_1
  "2" = paste0(crop_size, "+1100+950"),  # IPF4
  "3" = paste0(crop_size, "+200+250"), # IPF_3a
  "4" = paste0(crop_size, "+300+1000"))   # IPF_3b

se.he_crop <- se.he
se.he_crop <- CropImages(se.he_crop, crop_window_list) #list("1"="500x500+500+500")) #crop_window_list
ImagePlot(se.he_crop)

p <- FeatureOverlay(se.he_crop, 
                    sampleids = 1:(se.he_crop$sample_id %>% unique() %>% length()),
                    ncols = (se.he_crop$sample_id %>% unique() %>% length()),
                    features = "f14_subclusters2", 
                    label.by = "sample_name",
                    cols = col_clusters2, 
                    min.cutoff = 1,
                    max.cutoff = 4,
                    pt.alpha = 0.8,
                    pt.size = 3.5) &
  theme(legend.position = "bottom", 
        plot.title = element_blank(),
        plot.margin = margin(4, 5, 2.5, 5))

# p <- ggrastr::rasterize(p, layers = "Point", dpi = fig_res)
pdf(file = file.path(DIR_FIG, "hs_visium_fig_A-NMF30_F14-high_subcluster_spatial_he_crops.pdf"), width = 13, height = 5);p;dev.off()


rm(se.he_crop);rm(se.he)


##### 3C: F14hi-C0 radial distance cell densities ##### 
#' See script: "hs_visium_radial_distance.R"



##### 3D: F14hi-C0 radial distance gene expression ##### 
#' See script: "hs_visium_radial_distance.R" for full plotting script
 
#' (Underlying data:)
data_path <- file.path(DIR_RES, "objects/radial_distance/excl_singletons")
p_dat_g_concat <- read.csv(file.path(data_path, "hs_visium_A-NMF30-F14hi-C0_distance_gene_top_pos_neg_data.csv"))

gene_cor_df <- read.csv(file.path(data_path, "hs_visium_A-NMF30-F14hi-C0_distance_gene_cor_500um_180923.csv"))



##### 3E: F14hi-C0 neighboring clusters ##### 
levels(se.subset$f14_nbs_clusters)
se.subset$f14_nbs_clusters <- factor(se.subset$f14_nbs_clusters, levels = c("F14_C0", 0:5, "other"))
cluster_names <- levels(se.subset$f14_nbs_clusters)
n_clusters <- length(cluster_names)
se.subset_f14hi <- SubsetSTData(se.subset, f14_nbs_clusters != "other")


#' i. Gene markers dotplot:
makers_f14_subcluster_nbs$updown <- ifelse(makers_f14_subcluster_nbs$avg_log2FC>0, "up", "down")
genes_plot_dp <- makers_f14_subcluster_nbs %>%
  filter(p_val_adj<0.05) %>% 
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group=TRUE) %>%
  top_n(5, (avg_log2FC))

p_marker_dp <- DotPlot(se.subset_f14hi, 
                       features = genes_plot_dp$gene %>% unique(), 
                       group.by = "f14_nbs_clusters", 
                       col.max = 1, col.min = -1,
                       # dot.scale = 4,
                       # cols="RdYlBu", 
                       scale = T) + 
  scale_color_gradientn(colours = c(rev(col_scale_mako[1:8]), col_scale_acton[1:8])) +
  # scale_color_gradientn(colours = col_scale_div_expr) +
  theme_dotplot +
  theme(axis.text.y = element_text(face = "plain"),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, face="italic"),
        # axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face="italic"), 
        legend.position = "bottom"); p_marker_dp


#' ii. c2l mapping dotplot:
cell_types_select <- grep("c2l_", colnames(se.subset_f14hi@meta.data), value = T)

dat <- se.subset_f14hi@meta.data %>% 
  pivot_longer(cols = all_of(cell_types_select), names_to = "cell_type", values_to = "cell_density")
dat_summary <- dat %>% 
  select(f14_nbs_clusters, cell_type, cell_density) %>% 
  group_by(f14_nbs_clusters, cell_type) %>% 
  summarise(
         avg_density = mean(cell_density),
         pct_expr = (mean(cell_density > 0.5))*100)

dat_summary <- dat_summary %>% 
  group_by(cell_type) %>% 
  mutate(avg_density_scaled = Scale01(avg_density),
         avg_density_zscore = scale(avg_density)[,1])
dat_summary$avg_density_zscore <- as.numeric(dat_summary$avg_density_zscore)
dat_summary$cell_type2 <- gsub("\\.\\.", "_", dat_summary$cell_type)
dat_summary$cell_type2 <- gsub("\\.$", "", dat_summary$cell_type2)
dat_summary <- merge(dat_summary, habermann_cell_names, by.x='cell_type2', by.y='c2l_col_names')


color_lims <- round(max(abs(dat_summary$avg_density_zscore))+0.4)
p_celltype_dp_all <- ggplot(
       dat_summary,
       aes(x=cell_name_special, y = f14_nbs_clusters, fill = avg_density_zscore, size = pct_expr)) +
  geom_point(shape=21) +
  scale_fill_gradientn(colours = c(rev(col_scale_mako[1:8]), col_scale_acton[1:8]),
                       limits = c(-color_lims, color_lims)) +
  theme_dotplot +
  theme(
    axis.text.y = element_text(face = "plain"),
    axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
    legend.position = "bottom"); p_celltype_dp_all

cell_types_plot <- dat_summary %>% 
  group_by(f14_nbs_clusters) %>% 
  # arrange(avg_density_scaled) %>% 
  slice_max(order_by = avg_density_zscore, n = 3)
cell_types_plot <- cell_types_plot$cell_name_special %>% unique()

dat_summary_subset <- subset(dat_summary, cell_name_special %in% cell_types_plot)
dat_summary_subset$cell_name_special <- factor(dat_summary_subset$cell_name_special, levels = cell_types_plot)

color_lims <- max(abs(dat_summary_subset$avg_density_zscore))
p_celltype_dp <- ggplot(
  dat_summary_subset,
  aes(x=cell_name_special, y = f14_nbs_clusters, color = avg_density_zscore, size = pct_expr)) +
  geom_point() +
  scale_color_gradientn(colours = c(rev(col_scale_mako[1:8]), col_scale_acton[1:8]),
                       limits = c(-color_lims, color_lims)) +
  theme_dotplot +
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
    legend.position = "bottom"); p_celltype_dp


pdf(file = file.path(DIR_FIG, "hs_visium_fig_A-NMF30_F14-high_C0_nbs_cluster_c2l_all.pdf"), width = 6, height = 4)
p_celltype_dp_all
dev.off()

pdf(file = file.path(DIR_FIG, "hs_visium_fig_A-NMF30_F14-high_C0_nbs_cluster_c2l_top3.pdf"), width = 4, height = 4)
p_celltype_dp
dev.off()



#' iii. donor prop bars:
# Plot cluster proportions in subjects 
d_prop <- se.subset_f14hi@meta.data %>%
  group_by(f14_nbs_clusters, subject_alias) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

p_prop <- ggplot(d_prop, aes(x=f14_nbs_clusters, y=n, fill=reorder(subject_alias, desc(subject_alias)))) +
  geom_bar(stat = 'identity', colour=NA, position = "stack", width = 0.8) +
  scale_fill_manual(values = cols_donor) +
  theme_dotplot +
  labs(y="spots") +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        # axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        axis.title.x = element_text(size=10),
        legend.position = "right", 
        legend.key.size = unit(12, "pt"),
        legend.title = element_blank()) +
  coord_flip();p_prop


#' Export grid
# (p_marker_dp | p_celltype_dp)
pdf(file = file.path(DIR_FIG, "hs_visium_fig_A-NMF30_F14-high_C0_nbs_cluster_data_new.pdf"), width = 10, height = 4)
(p_marker_dp | p_celltype_dp | p_prop) + patchwork::plot_layout(widths = c(4,3,1)) & theme(plot.margin = margin(2, 0, 2, 0))
dev.off()

##### 3F: IPA ##### 
#' IPA results - plotted elsewhere

##### 3G: NicheNet ##### 
#' Plotted in NicheNet script (hs_visium_nichenet.R)


##### 3H: APOE ##### 
#' Plotted in NicheNet script (hs_visium_nichenet.R)


#### Mixed ####
#' The rest of main figs plotted in separate scripts.

##### Extended Data: A-NMF30-F14-high subcluster dotplot ####
# Plot marker dotplot 
makers_f14_subcluster$updown <- ifelse(makers_f14_subcluster$avg_log2FC>0, "up", "down")
genes_plot_dp <- makers_f14_subcluster %>%
  filter(p_val_adj<0.05) %>% 
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group=TRUE) %>%
  top_n(10, (avg_log2FC))

p_marker_dp <- DotPlot(se.f14high, 
                       features = genes_plot_dp$gene, 
                       group.by = "f14_subclusters2", 
                       dot.scale = 4, 
                       col.max = 1.5,
                       col.min = -1.5,
                       # cols="RdYlBu", 
                       # cols=c(col_scale_acton[5], "white", col_scale_mako[5]),
                       scale = T) + 
  scale_color_gradientn(colours = c(rev(col_scale_mako[1:5]), col_scale_acton[1:8])) +
  labs(y="Cluster") +
  theme_dotplot +
  theme(axis.text.y = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain", size=11),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, face="italic"),
        legend.position = "bottom");p_marker_dp


# Plot celltype dotplot 
cell_types_plot <- c(grep("KRT5", colnames(c2l_subset), value = T),
                     grep("B.Cell", colnames(c2l_subset), value = T),
                     grep("Fibroblast", colnames(c2l_subset), value = T),
                     grep("_AT[1-2]", colnames(c2l_subset), value = T),
                     grep("_Macrophage", colnames(c2l_subset), value = T),
                     grep("Myofibroblast", colnames(c2l_subset), value = T),
                     grep("MUC5", colnames(c2l_subset), value = T),
                     grep("Basal", colnames(c2l_subset), value = T),
                     grep("Ciliated", colnames(c2l_subset), value = T))

p_celltype_dp <- DotPlot(se.f14high, 
                         features = cell_types_plot, #colnames(c2l_subset), 
                         group.by = "f14_subclusters2", 
                         dot.min = 0.1,
                         col.max = 1.5,
                         col.min = -1.5,
                         # cols = "RdYlBu",
                         scale = T) + 
  scale_color_gradientn(colours = c(rev(col_scale_mako[1:5]), col_scale_acton[1:8])) +
  theme_dotplot +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
    legend.position = "bottom");p_celltype_dp


# Plot cluster proportions in subjects 
d_prop <- se.f14high@meta.data %>%
  group_by(f14_subclusters2, subject_alias) %>% # f14_subclusters
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

p_prop <- ggplot(d_prop, aes(x=f14_subclusters2, y=n, fill=reorder(subject_alias, desc(subject_alias)))) +
  geom_bar(stat = 'identity', colour=NA, position = "stack", width = 0.8) +
  scale_fill_manual(values = cols_donor) +
  theme_dotplot +
  labs(y="spots") +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        axis.title.x = element_text(size=10),
        legend.position = "right", 
        legend.key.size = unit(12, "pt"),
        legend.title = element_blank()) +
  coord_flip();p_prop

# Save grid
pdf(file = file.path(DIR_FIG, "hs_visium_fig_A-NMF30_F14-high_subcluster_data.pdf"), width = 15, height = 4)
(p_umap2 | p_marker_dp | p_celltype_dp | p_prop) + patchwork::plot_layout(widths = c(2,7,2,1)) & theme(plot.margin = margin(2, 0, 2, 0))
dev.off()


# Cluster subject prop, export
d_prop_export <- se.f14high@meta.data %>%
  add_count(sample_name, name = "n_total_spots") %>% 
  group_by(f14_subclusters2, subject_alias, sample_name) %>%
  mutate(n = n()) %>%
  mutate(freq = n / sum(n)) %>% 
  slice(1) %>% 
  select(f14_subclusters2, subject_alias, sample_name, n_total_spots, n, freq)

write.csv(d_prop_export, file = file.path(DIR_RES, "objects", "hs_visium_A-NMF30_F14-high_subcluster_summary.csv"))

##### A-NMF30-F14-high-C0 nbs clusters #####
# levels(se.subset$f14_nbs_clusters)
# se.subset$f14_nbs_clusters <- factor(se.subset$f14_nbs_clusters, levels = c("F14_C0", 0:5, "other"))
# 
# cluster_names <- levels(se.subset$f14_nbs_clusters)
# n_clusters <- length(cluster_names)

# col_clusters <- setNames(c(scico::scico(10, palette = "acton")[2], cols_d3_20[1:6], "grey80"), 
#                          nm = cluster_names)
# se.subset_f14hi <- SubsetSTData(se.subset, f14_nbs_clusters != "other")

# Plot marker dotplot 
makers_f14_subcluster_nbs$updown <- ifelse(makers_f14_subcluster_nbs$avg_log2FC>0, "up", "down")
genes_plot_dp <- makers_f14_subcluster_nbs %>%
  filter(p_val_adj<0.05) %>% 
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group=TRUE) %>%
  top_n(5, (avg_log2FC))

p_marker_dp <- DotPlot(se.subset_f14hi, 
                       features = genes_plot_dp$gene %>% unique(), 
                       group.by = "f14_nbs_clusters", 
                       col.max = 1, col.min = -1,
                       dot.scale = 4,
                       # cols="RdYlBu", 
                       scale = T) + 
  scale_color_gradientn(colours = c(rev(col_scale_mako[1:5]), col_scale_acton[1:8])) +
  theme_dotplot +
  theme(axis.text.y = element_text(face = "plain"),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, face="italic"),
        # axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face="italic"), 
        legend.position = "bottom");


# Plot celltype dotplot 
d_cell_cluster <- se.subset_f14hi@meta.data %>% 
  select(all_of(c("f14_nbs_clusters", cell_type_col_names))) %>% 
  rownames_to_column(var = "barcode") %>% 
  pivot_longer(cols = cell_type_col_names, names_to = "cell_type", values_to = "value") %>% 
  group_by(f14_nbs_clusters, cell_type) %>% 
  summarise(cell_value_sum = sum(value),
            cell_value_avg = mean(value))

cell_types_plot <- d_cell_cluster %>% 
  group_by(f14_nbs_clusters) %>% 
  slice_max(order_by = cell_value_avg, n = 5) %>% 
  ungroup() %>% 
  select(cell_type) %>% 
  unique()
cell_types_plot <- cell_types_plot$cell_type


p_celltype_dp <- DotPlot(se.subset_f14hi, 
                         features = cell_types_plot, #cell_type_col_names,
                         group.by = "f14_nbs_clusters", 
                         dot.min = 0.1, 
                         col.max = 1, col.min = -1,
                         # cols = "RdYlBu",
                         scale = T) + 
  scale_color_gradientn(colours = c(rev(col_scale_mako[1:5]), col_scale_acton[1:8])) +
  theme_dotplot +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
    # axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
    panel.grid = element_blank(),
    legend.position = "bottom");


# Plot cluster proportions in subjects 
d_prop <- se.subset_f14hi@meta.data %>%
  group_by(f14_nbs_clusters, subject_alias) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

p_prop <- ggplot(d_prop, aes(x=f14_nbs_clusters, y=n, fill=reorder(subject_alias, desc(subject_alias)))) +
  geom_bar(stat = 'identity', colour=NA, position = "stack", width = 0.8) +
  scale_fill_manual(values = cols_donor) +
  theme_dotplot +
  labs(y="spots") +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        # axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        axis.title.x = element_text(size=10),
        legend.position = "right", 
        legend.key.size = unit(12, "pt"),
        legend.title = element_blank()) +
  coord_flip();

# Save grid
pdf(file = file.path(DIR_FIG, "hs_visium_fig_A-NMF30_F14-high_C0_nbs_cluster_data.pdf"), width = 10, height = 4)
(p_marker_dp | p_celltype_dp | p_prop) + patchwork::plot_layout(widths = c(4,2.5,1)) & theme(plot.margin = margin(2, 0, 2, 0))
dev.off()


# Plot cluster proportions radial plot 
d_prop_rad <- se.subset_f14hi@meta.data %>%
  filter(f14_nbs_clusters != "F14_C0") %>% 
  group_by(f14_nbs_clusters, subject_alias) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

p_prop_rad <- ggplot(d_prop_rad, aes(f14_nbs_clusters, n, fill = subject_alias)) +
  # geom_hline(yintercept = c(0.01, 0.25, 0.5), color = "black", linewidth=0.1) +
  geom_hline(yintercept = c(1, 1e3, 2e3), color = "black", linewidth=0.1) +
  geom_bar(width = 0.9, stat = "identity") +
  scale_fill_manual(values = cols_donor) +
  labs(fill="") +
  theme_void() +
  theme(axis.ticks = element_blank(),
        axis.text = element_text(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "right",
        axis.line = element_blank()) +
  coord_polar()

pdf(file = file.path(DIR_FIG, "hs_visium_fig_A-NMF30_F14-high_C0_nbs_cluster_prop_radial.pdf"), width = 3, height = 2)
p_prop_rad
dev.off()


#### B-NMF30-TLS ####
##### H&E panels ####
tls_factors <- list(IPF_1 = "factor_10", IPF_2 = "factor_15", IPF_3 = "factor_17")
tls_sample_select <- list(IPF_1 = "IPF_1.TD012.B2.2", IPF_2 = "IPF_2.TD022.B3.1", IPF_3 = "IPF_3.TD031.B1.2")

se.he_list <- lapply(names(tls_sample_select), function(s){
  message(s)
  se.he <- SubsetSTData(se_subset_split[[s]], expression = sample_name %in% tls_sample_select[[s]])
  se.he <- LoadImages(se.he, xdim = 2e3)
})
names(se.he_list) <- names(tls_sample_select)
crop_window_list <- list(IPF_1 = "600x600+250+1250", IPF_2 = "600x600+960+420", IPF_3 = "600x600+960+1300")


p_list <- lapply(names(tls_sample_select), function(s){
  se.he_s_crop <- se.he_list[[s]]
  se.he_s_crop <- CropImages(se.he_s_crop, crop.geometry.list = list("1" = crop_window_list[[s]]))
  # ImagePlot(se.he_s_crop, annotate = F)
  p <- FeatureOverlay(se.he_s_crop, features = tls_factors[[s]], 
                      cols = rev(col_scale_mako), 
                      # max.cutoff = 6,
                      add.alpha = T, pt.size = 2.5) & 
    labs(fill=tls_factors[[s]], title = s) &
    theme(legend.position = "bottom", plot.subtitle = element_blank(), plot.title = element_text(size=9))
  p <- ggrastr::rasterize(p, layers = "Point", dpi = fig_res)
})

pdf(file = file.path(DIR_FIG, "hs_visium_fig_B-NMF30_TLS_factors_selected_crop.pdf"), width = 10, height = 4.2)
patchwork::wrap_plots(p_list, nrow = 1) & theme(plot.margin = margin(4, 2.5, 2.5, 2.5))
dev.off()


##### Gene heatmap ####
tls_f_genes <- lapply(names(tls_factors), function(s){
  se_subset_split[[s]]@reductions$NMF@feature.loadings[,tls_factors[[s]]]
})
gene_ovelap <- Reduce(intersect, lapply(1:length(tls_f_genes), function(i){names(tls_f_genes[[i]])}))

tls_f_df <- data.frame(gene = gene_ovelap,
                       IPF_1.F10 = tls_f_genes[[1]][gene_ovelap],
                       IPF_2.F15 = tls_f_genes[[2]][gene_ovelap],
                       IPF_3.F17 = tls_f_genes[[3]][gene_ovelap])

tls_f_df[,2:ncol(tls_f_df)] <- apply(tls_f_df[,2:ncol(tls_f_df)], 2, Scale01)
tls_f_df_filt <- tls_f_df
tls_f_df_filt$gene <- NULL
rsums <- rowSums(tls_f_df_filt)
cutoff <- as.numeric(quantile(rsums, .95)[1])
tls_f_df_filt <- tls_f_df_filt[rsums>cutoff,]
head(tls_f_df_filt);dim(tls_f_df_filt)

pdf(file = file.path(DIR_FIG, "hs_visium_fig_B-NMF30_TLS_factors_gene_hm.pdf"), width = 12, height = 3, useDingbats = F)
print(pheatmap(t(tls_f_df_filt), 
               # main = "Strongest gene contributors to the TLS factors across IPF donors",
               display_numbers = F, 
               cluster_rows = F, 
               treeheight_col = 20, 
               color = col_scale_mako, 
               cellwidth = 14, 
               cellheight = 14, 
               angle_col = 90))
dev.off()










