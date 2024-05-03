#' [mm_visium_nmf.R]
#'
#' Run, analyse, and plot NMF results 
#'
#'
#' Aug 2022, L. Franzén [lovisa.franzen@scilifelab.se]

#### Set up ####
##### Define params. ####
set.seed(1)
SPECIES <- "mouse"
DIR_ROOT <- "/home/st-analysis_home/lovisa.franzen/analysis/lung/spatial-lung-fibrosis"  #getwd()
DIR_DATA <- file.path(DIR_ROOT, "data", SPECIES, "visium")
DIR_RES <- file.path(DIR_ROOT, "results", SPECIES)
fig_res <- 300

##### Load libs ####
library(STutility)
library(RColorBrewer)
library(patchwork)
library(dplyr)
library(magrittr)
library(writexl)
library(readxl)
library(gprofiler2)

##### Other ####
# source(file.path(DIR_ROOT, "scripts", "colors.R"))
source(file.path(DIR_ROOT, "scripts", "custom_functions.R"))
source(file.path(DIR_ROOT, "scripts", "custom_colors.R"))
theme_custom <- theme(axis.title.x = element_blank())
cols_group <- RColorBrewer::brewer.pal(4, "Spectral")[c(1,4, 2,3)]
names(cols_group) <- c("d21_BLM", "d21_CTRL","d7_BLM", "d7_CTRL")


##### Read objects ####
# fname <- paste0("mm_visium_preproc_se_obj.rds")
fname <- "mm_visium_preproc_se_obj_nmf30_c2l.rds" # if re-running script 
se <- readRDS(file = file.path(DIR_RES, "objects", fname))

# Process annotation labels
se$annotation2 <- ifelse(se$annotation == "Inflammation", paste0(se$annotation, "_", se$day), se$annotation)


#### NMF: All data ####
##### Run NMF ####
# n_factors <- 20  # previous setting
n_factors <- 30  # new setting to match human NMF and to obtain more refined NMF results
DIR_FIG_NMF <- file.path(DIR_RES, "figures", paste0("NMF", n_factors))
DIR_OBJ_NMF <- file.path(DIR_RES, "objects", paste0("NMF", n_factors))
dir.create(DIR_OBJ_NMF); dir.create(DIR_FIG_NMF)
file_save_prefix <- paste0("mm_visium_nmf", n_factors, "_")

factor_names <- paste0("factor_", 1:n_factors)

se <- RunNMF(se, nfactors = n_factors)


##### NMF Results ####
fib_genes <- c(grep("^Col[0-9]", rownames(se@reductions$NMF@feature.loadings), value = T), "Fn1")
bcell_genes <- c("Jchain", "Igha", "Igkc")
at1_genes <- c("Hopx", "Ager")
at2_genes <- c("Sftpc", "Sftpb")

factor_fib <- names(sort(colSums(se@reductions$NMF@feature.loadings[fib_genes, ]), decreasing = T))[1]
factor_bcell <- names(sort(colSums(se@reductions$NMF@feature.loadings[bcell_genes, ]), decreasing = T))[1]
factor_at1 <- names(sort(colSums(se@reductions$NMF@feature.loadings[at1_genes, ]), decreasing = T))[1]
factor_at2 <- names(sort(colSums(se@reductions$NMF@feature.loadings[at2_genes, ]), decreasing = T))[1]


#' Save top genes for each factor
# se@reductions$NMF@feature.loadings
factor_gene_loadings <- lapply(1:n_factors, function(factor_x){
  message(paste("Factor", factor_x))
  feat_loads <- as.data.frame(se@reductions$NMF@feature.loadings[,paste0("factor_",factor_x)])
  colnames(feat_loads) <- "gene_loading"
  feat_loads$gene <- rownames(feat_loads)
  feat_loads$factor <- factor_x
  top_n_genes <- 100
  feat_loads <- feat_loads %>%
    dplyr::slice_max(order_by = gene_loading, n = top_n_genes) %>%
    mutate(rank = dense_rank(desc(gene_loading)))
})

write_xlsx(
  factor_gene_loadings,
  # path = file.path(DIR_OBJ_NMF, "mm_visium_nmf_1-20_top100_gene_loadings.xlsx"),
  path = file.path(DIR_OBJ_NMF, paste0(file_save_prefix, "top100_gene_loadings.xlsx")),
  col_names = TRUE,
  format_headers = TRUE,
  use_zip64 = FALSE
)
factor_gene_loadings <- lapply(1:n_factors, function(factor_x){
  read_xlsx(path = file.path(DIR_OBJ_NMF, paste0(file_save_prefix, "top100_gene_loadings.xlsx")), sheet = factor_x)
})

##### gProfiler ####
library(gprofiler2)

org <- "mmusculus"
# fname <- "mm_visium_nmf_1-20_gProfiler_top25genes"
fname <- paste0(file_save_prefix, "gProfiler_top25genes")
n_genes <- 25
gene_query_list <- list()
gostres_list <- list()

for(i in 1:n_factors){
  gene_query_list[[i]] <- subset(factor_gene_loadings[[i]], rank <= n_genes)$gene
}
names(gene_query_list) <- factor_names

for (f in names(gene_query_list)){
  message(f)
  gostres_list[[f]] <- gprofiler2::gost(query = gene_query_list[[f]],
                                        organism = org, 
                                        ordered_query = T)
}

gostres_export_list <- list()
for (f in names(gostres_list)){
  gostres_export_list[[f]] <- gostres_list[[f]]$result
  gostres_export_list[[f]]$query <- f
}
write_xlsx(
  x = gostres_export_list,
  path = file.path(DIR_OBJ_NMF, paste0(fname, ".xlsx")),
  col_names = TRUE,
  format_headers = TRUE,
  use_zip64 = FALSE
)

pdf(file = file.path(DIR_FIG_NMF, paste0(fname, ".pdf")), width = 10, height = 8, useDingbats = F)
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




##### Plot NMF ####
###### Spatial heatmaps ####
selected_samples <- c(2, 8, 13, 19)
plot_list <- lapply(1:n_factors, function(i) {
  print(i)
  FactorPlot(seurat.object = se, 
             indices = selected_samples,
             factor = i, 
             col.scale = rev(RColorBrewer::brewer.pal(10, "Spectral")))
})
png(file.path(DIR_FIG_NMF, paste0(file_save_prefix, "spatial.png")), width = 24*fig_res, height = 40*fig_res, res = fig_res)
cowplot::plot_grid(plotlist = plot_list[1:n_factors], ncol = 2)
dev.off()


###### Gene loadings ####
factor_gene_loading_plot <- lapply(1:n_factors, function(i){
  n_genes <- 10
  dat_plot <- factor_gene_loadings[[i]]
  p <- ggplot(dat_plot, aes(x=rank, y=gene_loading)) +
    geom_line(linewidth=.5) +
    geom_point(data = dat_plot[1:n_genes,], mapping = aes(x=rank, y=gene_loading, color=reorder(gene, rank))) +
    scale_color_manual(values = RColorBrewer::brewer.pal(n_genes, "Spectral")) +
    labs(title = paste("Factor", i), x="gene rank", y="loading", color="") +
    theme_bw() +
    theme(legend.position = "right", legend.justification = "left", 
          legend.text = element_text(size=7), 
          legend.box.margin=margin(t = -.1, l = -.25, unit='cm'),  # pos right
          # legend.box.margin=margin(t = -.1, l = -1, unit='cm'),  # pos bottom
          # axis.text.x = element_blank(),
          plot.title = element_text(hjust=0.5, face="bold"),
          panel.grid = element_blank())
})
png(file.path(DIR_FIG_NMF, paste0(file_save_prefix, "gene_loadings_.png")), width = 17*fig_res, height = 12*fig_res, res = fig_res)
cowplot::plot_grid(plotlist = factor_gene_loading_plot[1:n_factors], ncol = 5)
dev.off()


###### NMF vs histopath annotations ####
png(file.path(DIR_FIG_NMF, paste0(file_save_prefix, "histopath.png")), width = 12*fig_res, height = 20*fig_res, res = fig_res)
VlnPlot(se, 
        features = factor_names, 
        group.by = "annotation",
        # split.by = "group", cols = cols_group, 
        cols = RColorBrewer::brewer.pal(6, "Set1"),
        ncol = 5,
        pt.size = 0) & theme_custom & NoLegend()
dev.off()

png(file.path(DIR_FIG_NMF, paste0(file_save_prefix, "histopath_group.png")), width = 14*fig_res, height = 20*fig_res, res = fig_res)
VlnPlot(se, 
        features = factor_names,
        group.by = "annotation",
        split.by = "group", cols = cols_group,
        ncol = 5,
        pt.size = 0) & theme_custom & NoLegend()
dev.off()


###### Rank plot ####
# Sample group
plot_list <- lapply(1:n_factors, function(i) {
  p <- FactorRankSpotPlot(seurat.object = se, factor = i, split.by = "sample_name", group.by = "group", color.group = cols_group)
})
png(file.path(DIR_FIG_NMF, paste0(file_save_prefix, "ranks.png")), width = 24*fig_res, height = 12*fig_res, res = fig_res)
cowplot::plot_grid(plotlist = plot_list[1:n_factors], ncol = 5)
dev.off()

# Annotation - only annotated samples included
se_annotated <- SubsetSTData(se, annotation != "NA")
plot_list <- lapply(1:n_factors, function(i) {
  p <- FactorRankSpotPlot(seurat.object = se_annotated, factor = i, split.by = "annotation2", group.by = "annotation2", color.group = cols_annotation)
})
png(file.path(DIR_FIG_NMF, paste0(file_save_prefix, "ranks_histopath.png")), width = 24*fig_res, height = 12*fig_res, res = fig_res)
cowplot::plot_grid(plotlist = plot_list[1:n_factors], ncol = 5)
dev.off()


###### Factor plot grid ####
plot_list <- lapply(1:n_factors, function(factor_plot) {
  p1 <- ST.DimPlot(se, 
                   # indices = c(1:3, 13:15, 4:6, 16:18, 7:9, 19:21),
                   indices = c(1:9, 13:21),
                   dims = factor_plot, 
                   min.cutoff = 0.2,
                   label.by = "sample_name",
                   center.zero = F, 
                   reduction = "NMF", 
                   ncol = 9, 
                   pt.size = 0.8, 
                   pt.border = F,
                   cols = rev(RColorBrewer::brewer.pal(10, "Spectral")),
                   # cols = viridis::magma(10, direction = -1),
                   show.sb = F)  &
    labs(fill = "factor\nactivity", title = paste0("Factor ", factor_plot)) & 
    theme(legend.position = "right", plot.title = element_text(face = "bold", size = 16))
  p2 <- FactorGeneLoadingPlotHorizontal(se, 
                                        factor = factor_plot, 
                                        topn = 30) + 
    labs(title = "Genes contributing to factor") & 
    theme(axis.text.x = element_text(size=10), panel.grid = element_blank(),
          plot.margin = margin(10, 10, 0, 10, unit = "pt"))
  
  p3 <- FactorRankSpotPlot(se, factor = factor_plot, split.by = "sample_name", group.by = "group", color.group = cols_group, plot.title = "Spots ranked by factor activity")
  
  p_bottom <- cowplot::plot_grid(p2, p3, nrow = 1)
  p <- cowplot::plot_grid(p1, p_bottom, ncol = 1, rel_heights = c(3, 2))
})

dir.create(path = file.path(DIR_FIG_NMF, paste0(file_save_prefix, "grid")))
for(i in 1:length(plot_list)){
  png(file = file.path(DIR_FIG_NMF, paste0(file_save_prefix, "grid"), paste0("mm_visium_nmf_grid_f", i, ".png")),
      width = 16*fig_res, height = 7*fig_res, res = fig_res)
  print(plot_list[[i]])
  dev.off()
}

# pdf(file = file.path(DIR_FIG_NMF, "mm_visium_nmf_1-20_grid.pdf"), width = 16, height = 7, pointsize = 4, fg = "white")
# print(plot_list[1])
# dev.off()



###### UMAP & clusters & NMF ####
n_clusters <- length(unique(se$seurat_clusters))
col_scale <- viridis::rocket(n_clusters, direction = -1)
plot_list <- lapply(1:n_factors, function(factor_plot) {
  factor_name <- paste0("factor_", factor_plot)
  p1 <- VlnPlot(se, features = factor_name, group.by = "seurat_clusters", cols = rev(col_scale), pt.size = 0) & 
    labs(title = paste("Factor", factor_plot), x = "Cluster") &
    coord_flip() & 
    NoLegend() & 
    theme(axis.text.x = element_text(angle=0, hjust = 0.5))
  p2 <- FeaturePlot(se, features = factor_name, reduction = "umap.harmony", cols = col_scale, min.cutoff = .2, raster = T) &
    labs(title = "", x = "UMAP 1", y="UMAP 2") & 
    NoLegend()
  p <- p1 + p2 + plot_layout(ncol = 2, widths = c(1,2))
})

png(file = file.path(DIR_FIG_NMF, paste0(file_save_prefix, "umap_clusters.png")), width = 24*fig_res, height = 18*fig_res, res = fig_res)
cowplot::plot_grid(plotlist = plot_list[1:n_factors], ncol = 5)
dev.off()


##### Factor ~ Factor correlations ##### 
factor_cor_data <- bind_cols(se@meta.data %>% select(sample_name, day, condition, group),
                             se@reductions$NMF@cell.embeddings[, factor_names])

groups_use <- c("all", "CTRL", "d7_BLM", "d21_BLM")

factor_cor_data_group <- setNames(lapply(groups_use, function(g){
  if (g=="all") {
    factor_cor_subset <- factor_cor_data
  } else if (g=="CTRL") {
    factor_cor_subset <- subset(factor_cor_data, condition == "control")
  } else {
    factor_cor_subset <- subset(factor_cor_data, group == g)
  }
  factor_cor_hm <- cor(factor_cor_subset[, factor_names])
  diag(factor_cor_hm) <- 0
  factor_cor_hm <- as.data.frame(factor_cor_hm)
}), 
nm = groups_use)

# Plot as pheatmaps
pal_length <- 11
ph_colors <- RColorBrewer::brewer.pal(pal_length, "RdBu") %>% rev()
abs_max_val <- lapply(factor_cor_data_group, function(m){max(abs(m))}) %>% unlist() %>% max()

pdf(file = file.path(DIR_FIG_NMF, paste0(file_save_prefix, "factor_cor_heatmap_per_group.pdf")), width = 8, height = 8, useDingbats = F)
for(g in groups_use){
  factor_cor_hm <- factor_cor_data_group[[g]]
  
  abs_max_val <- max(abs(factor_cor_hm))
  ph_breaks <- c(seq(-abs_max_val, 0, length.out=ceiling(pal_length/2) + 1), 
                 seq(abs_max_val/pal_length, abs_max_val, length.out=floor(pal_length/2)))
  
  print(
    pheatmap::pheatmap(factor_cor_hm, 
                       cellwidth = 10, 
                       cellheight = 10, 
                       color = ph_colors, 
                       breaks = ph_breaks, 
                       main = g)
  )
}
dev.off()


##### Factor ~ Cell correlation #####
c2l_all <- read.csv(file.path(DIR_ROOT, "data", SPECIES, "sc_deconvolution_strunz", "compiled_all_samples_cell_abundances.csv"), row.names = 1)
c2l_names <- colnames(c2l_all)
new_c2l_names <- gsub("c2l_", "", gsub("[.]$", "", gsub("..", "_", c2l_names, fixed = TRUE)))

se <- AddMetaData(se, c2l_all)


factor_cell_cor_data <- bind_cols(se@meta.data %>% select(sample_name, day, condition, group, matches(c2l_names)),
                                  se@reductions$NMF@cell.embeddings[, factor_names])

colnames(factor_cell_cor_data)[colnames(factor_cell_cor_data) %in% c2l_names] <- new_c2l_names

saveRDS(factor_cell_cor_data, file = file.path(DIR_RES, "objects", paste0(file_save_prefix, "c2l_spot_data.rds")))


# Filter cells to show
new_c2l_names_filtered <- new_c2l_names[!new_c2l_names %in% c("NA", "low.quality.cells")]

groups_use <- c("all", "CTRL", "d7_BLM", "d21_BLM")
factor_cell_cor_data_group <- setNames(lapply(groups_use, function(g){
  if (g=="all") {
    factor_cell_cor_subset <- factor_cell_cor_data
  } else if (g=="CTRL") {
    factor_cell_cor_subset <- subset(factor_cell_cor_data, condition == "control")
  } else {
    factor_cell_cor_subset <- subset(factor_cell_cor_data, group == g)
  }
  factor_cell_density_hm <- cor(factor_cell_cor_subset[,c(new_c2l_names_filtered, factor_names)])
  diag(factor_cell_density_hm) <- 0
  factor_cell_density_hm <- as.data.frame(factor_cell_density_hm)
  factor_cell_density_hm <- factor_cell_density_hm[factor_names, new_c2l_names_filtered]
}), 
nm = groups_use)


saveRDS(factor_cell_cor_data_group, file = file.path(DIR_RES, "objects", paste0(file_save_prefix, "c2l_pearson_corr_matrices.rds")))


# Plot pheatmap
pal_length <- 11
ph_colors <- RColorBrewer::brewer.pal(pal_length, "RdBu") %>% rev()
abs_max_val <- lapply(factor_cell_cor_data_group, function(m){max(abs(m))}) %>% unlist() %>% max()

pdf(file = file.path(DIR_FIG_NMF, paste0(file_save_prefix, "c2l_res_cell_factor_cor_heatmap_per_group.pdf")), 
    width = 10, height = 10, useDingbats = F)
for(g in groups_use){
  factor_cell_density_hm <- factor_cell_cor_data_group[[g]]
  
  abs_max_val <- max(abs(factor_cell_density_hm))
  ph_breaks <- c(seq(-abs_max_val, 0, length.out=ceiling(pal_length/2) + 1), 
                 seq(abs_max_val/pal_length, abs_max_val, length.out=floor(pal_length/2)))
  
  print(
    pheatmap::pheatmap(t(factor_cell_density_hm), 
                       cellwidth = 10, 
                       cellheight = 10, 
                       color = ph_colors, 
                       breaks = ph_breaks, 
                       main = g)
  )
}
dev.off()


###### Factor ~ Krt8ADI correlation ######
factor_cell_cor_data_group_krt8 <- do.call(cbind, factor_cell_cor_data_group[2:4]) %>% select(contains("Krt8"))
factor_cell_cor_data_group_krt8

pal_length <- 11
ph_colors <- RColorBrewer::brewer.pal(pal_length, "RdBu") %>% rev()
abs_max_val <- max(factor_cell_cor_data_group_krt8)

pdf(file = file.path(DIR_FIG_NMF, paste0(file_save_prefix, "c2l_res_cell_factor_cor_heatmap_krt8ADI_group.pdf")), 
    width = 6, height = 6, useDingbats = F)
pheatmap::pheatmap(factor_cell_cor_data_group_krt8, 
                   cellwidth = 20, 
                   cellheight = 10, 
                   treeheight_row = 10,  
                   treeheight_col = 4,
                   color = ph_colors, 
                   breaks = ph_breaks)
dev.off()


##### Wrap up ####
# Save new se obj with new NMF
fname <- paste0("mm_visium_preproc_se_obj_nmf30_c2l.rds")
saveRDS(se, file = file.path(DIR_RES, "objects", fname))



#### NMF: Selecte sample(s) ####
#' Run NMF on separate sections to see if we can obtain a factor that contains
#' a more Krt8+ADI specific signature

##### Data: d21 ####
###### Subset ######
subset_id <- "d21"
se_subset <- SubsetSTData(se, day == "d21")

###### Run NMF ######
#' Define params
n_factors <- 30  # new setting to match human NMF and to obtain more refined NMF results
DIR_FIG_NMF <- file.path(DIR_RES, "figures", paste0("NMF", n_factors, "_", subset_id))
DIR_OBJ_NMF <- file.path(DIR_RES, "objects", paste0("NMF", n_factors, "_", subset_id))
dir.create(DIR_OBJ_NMF); dir.create(DIR_FIG_NMF)
file_save_prefix <- paste0("mm_visium_nmf", n_factors, "_", subset_id, "_")
factor_names <- paste0("factor_", 1:n_factors)

#' Run NMF
se_subset <- RunNMF(se_subset, nfactors = n_factors)

# Save new se obj with new NMF
fname <- paste0("mm_visium_preproc_se_obj_subset_", subset_id, "_nmf30_c2l.rds")
saveRDS(se_subset, file = file.path(DIR_RES, "objects", fname))
se_subset <- readRDS(file = file.path(DIR_RES, "objects", fname))


###### NMF Results ######
#' Look at factors that contains the Krt8+ADI marker genes most highly
adi_genes <- c(grep("^Krt", rownames(se_subset@reductions$NMF@feature.loadings), value = T), "Hbegf", "Plaur", "Areg", "Sprr1a", "Tnip3", "Ckn1", "Itgb6", "Cldn4", "End1")
adi_genes <- intersect(adi_genes, rownames(se_subset@reductions$NMF@feature.loadings)); adi_genes
factors_adi <- names(sort(colSums(se_subset@reductions$NMF@feature.loadings[adi_genes, ]), decreasing = T))[1:3]; factors_adi

# FactorPlot(se_subset, factor = 14, indices = 1:6, n.columns = 6, col.scale = col_scale_mako)
f_plot <- gsub("factor_", "", factors_adi) %>% as.numeric()
f_plot_list <- lapply(f_plot, function(f){
  p <- FactorPlot(se_subset, factor = f, indices = 1:12, n.columns = 6, col.scale = col_scale_mako)
})

png(file = file.path(DIR_FIG_NMF, paste0("mm_visium_nmf_potential_Krt8ADI_factors_spatial_gene.png")),
    width = 18*fig_res, height = 16*fig_res, res = fig_res)
wrap_plots(f_plot_list, ncol = 1)
dev.off()


plot_list <- lapply(f_plot, function(factor_plot) {
  p1 <- ST.DimPlot(se_subset, 
                   indices = 1:12,
                   dims = factor_plot, 
                   min.cutoff = 0.2,
                   label.by = "sample_name",
                   center.zero = F, 
                   reduction = "NMF", 
                   ncol = 6, 
                   pt.size = 0.8,
                   pt.border = F,
                   cols = col_scale_mako,
                   show.sb = F)  &
    labs(fill = "factor\nactivity", title = paste0("Factor ", factor_plot)) & 
    theme(legend.position = "right", plot.title = element_text(face = "bold", size = 16), aspect.ratio = 1)
  p2 <- FactorGeneLoadingPlotHorizontal(se_subset, 
                                        factor = factor_plot, 
                                        topn = 30) + 
    labs(title = "Genes contributing to factor") & 
    theme(axis.text.x = element_text(size=10), panel.grid = element_blank(),
          plot.margin = margin(10, 10, 0, 10, unit = "pt"))
  
  p3 <- FactorRankSpotPlot(se_subset, factor = factor_plot, 
                           split.by = "sample_name", group.by = "group", 
                           color.group = cols_group, plot.title = "Spots ranked by factor activity")
  
  p_bottom <- cowplot::plot_grid(p2, p3, nrow = 1)
  p <- cowplot::plot_grid(p1, p_bottom, ncol = 1, rel_heights = c(3, 2))
})


png(file = file.path(DIR_FIG_NMF, paste0("mm_visium_nmf_potential_Krt8ADI_factors_grid.png")),
    width = 12*fig_res, height = 7*3*fig_res, res = fig_res)
wrap_plots(plot_list, ncol = 1)
dev.off()

#' Look for collgen containing factors
fib_genes <- c(grep("^Col", rownames(se_subset@reductions$NMF@feature.loadings), value = T), "Fn1")
fib_genes <- intersect(fib_genes, rownames(se_subset@reductions$NMF@feature.loadings)); fib_genes
factors_fib <- names(sort(colSums(se_subset@reductions$NMF@feature.loadings[fib_genes, ]), decreasing = T))[1:3]; factors_fib


#' Save top genes for each factor
# se_subset@reductions$NMF@feature.loadings
factor_gene_loadings <- lapply(1:n_factors, function(factor_x){
  message(paste("Factor", factor_x))
  feat_loads <- as.data.frame(se_subset@reductions$NMF@feature.loadings[,paste0("factor_",factor_x)])
  colnames(feat_loads) <- "gene_loading"
  feat_loads$gene <- rownames(feat_loads)
  feat_loads$factor <- factor_x
  top_n_genes <- 100
  feat_loads <- feat_loads %>%
    dplyr::slice_max(order_by = gene_loading, n = top_n_genes) %>%
    mutate(rank = dense_rank(desc(gene_loading)))
})

write_xlsx(
  factor_gene_loadings,
  path = file.path(DIR_OBJ_NMF, paste0(file_save_prefix, "top100_gene_loadings.xlsx")),
  col_names = TRUE,
  format_headers = TRUE,
  use_zip64 = FALSE
)
factor_gene_loadings <- lapply(1:n_factors, function(factor_x){
  read_xlsx(path = file.path(DIR_OBJ_NMF, paste0(file_save_prefix, "top100_gene_loadings.xlsx")), sheet = factor_x)
})

###### gProfiler #####
org <- "mmusculus"
fname <- paste0(file_save_prefix, "gProfiler_top25genes")
n_genes <- 25
gene_query_list <- list()
gostres_list <- list()

for(i in 1:n_factors){
  gene_query_list[[i]] <- subset(factor_gene_loadings[[i]], rank <= n_genes)$gene
}
names(gene_query_list) <- factor_names

for (f in names(gene_query_list)){
  message(f)
  gostres_list[[f]] <- gprofiler2::gost(query = gene_query_list[[f]],
                                        organism = org, 
                                        ordered_query = T)
}

gostres_export_list <- list()
for (f in names(gostres_list)){
  gostres_export_list[[f]] <- gostres_list[[f]]$result
  gostres_export_list[[f]]$query <- f
}
write_xlsx(
  x = gostres_export_list,
  path = file.path(DIR_OBJ_NMF, paste0(fname, ".xlsx")),
  col_names = TRUE,
  format_headers = TRUE,
  use_zip64 = FALSE
)

pdf(file = file.path(DIR_FIG_NMF, paste0(fname, ".pdf")), width = 10, height = 8, useDingbats = F)
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


###### Factor plot grid ####
plot_list <- lapply(1:n_factors, function(factor_plot) {
  p1 <- ST.DimPlot(se_subset, 
                   indices = c(1:9),
                   dims = factor_plot, 
                   min.cutoff = 0.2,
                   label.by = "sample_name",
                   center.zero = F, 
                   reduction = "NMF", 
                   ncol = 6, 
                   pt.size = 0.8, 
                   pt.border = F,
                   cols = rev(RColorBrewer::brewer.pal(10, "Spectral")),
                   show.sb = F)  &
    labs(fill = "factor\nactivity", title = paste0("Factor ", factor_plot)) & 
    theme(aspect.ratio = 1, legend.position = "right", plot.title = element_text(face = "bold", size = 16))
  p2 <- FactorGeneLoadingPlotHorizontal(se_subset, 
                                        factor = factor_plot, 
                                        topn = 30) + 
    labs(title = "Genes contributing to factor") & 
    theme(axis.text.x = element_text(size=10), panel.grid = element_blank(),
          plot.margin = margin(10, 10, 0, 10, unit = "pt"))
  
  p3 <- FactorRankSpotPlot(se_subset, factor = factor_plot, split.by = "sample_name", 
                           group.by = "group", color.group = cols_group, plot.title = "Spots ranked by factor activity")
  
  p_bottom <- cowplot::plot_grid(p2, p3, nrow = 1)
  p <- cowplot::plot_grid(p1, p_bottom, ncol = 1, rel_heights = c(3, 2))
})

dir.create(path = file.path(DIR_FIG_NMF, paste0(file_save_prefix, "grid")))
for(i in 1:length(plot_list)){
  png(file = file.path(DIR_FIG_NMF, paste0(file_save_prefix, "grid"), paste0("mm_visium_nmf_grid_f", i, ".png")),
      width = 12*fig_res, height = 7*fig_res, res = fig_res)
  print(plot_list[[i]])
  dev.off()
}


###### Factor ~ Cell correlation ######
# c2l_all <- read.csv(file.path(DIR_ROOT, "data", SPECIES, "sc_deconvolution_strunz", "compiled_all_samples_cell_abundances.csv"), row.names = 1)
# se_subset <- AddMetaData(se_subset, c2l_all)
c2l_names <- grep("c2l_", colnames(se_subset@meta.data), value = T)
new_c2l_names <- gsub("c2l_", "", gsub("[.]$", "", gsub("..", "_", c2l_names, fixed = TRUE)))

factor_cell_cor_data <- bind_cols(se_subset@meta.data %>% select(sample_name, animal, day, condition, group, matches(c2l_names)),
                                  se_subset@reductions$NMF@cell.embeddings[, factor_names])
colnames(factor_cell_cor_data)[colnames(factor_cell_cor_data) %in% c2l_names] <- new_c2l_names
# write.csv(factor_cell_cor_data, file.path(DIR_RES, "..", "mouse", "objects", "NMF30_d21", "mm_visium_nmf30_d21_factorweight_cell2location_strunz_metadata.csv"))


# Filter cells to show
new_c2l_names_filtered <- new_c2l_names[!new_c2l_names %in% c("NA", "low.quality.cells")]

animal_ids <- se_subset$animal %>% unique()
groups_use <- c("all", "CTRL", "d21_BLM", "animal_id")
factor_cell_cor_data_group <- setNames(lapply(groups_use, function(g){
  message(g)
  if (g == "all") {
    factor_cell_cor_subset <- factor_cell_cor_data
  } else if (g == "CTRL") {
    factor_cell_cor_subset <- subset(factor_cell_cor_data, condition == "control")
  } else if (g == "d21_BLM") {
    factor_cell_cor_subset <- subset(factor_cell_cor_data, group == g)
  } else {
    factor_cell_density_hm <- setNames(lapply(animal_ids, function(a){
      message(a)
      factor_cell_cor_subset <- subset(factor_cell_cor_data, animal == a)
      factor_cell_density_hm <- cor(factor_cell_cor_subset[,c(new_c2l_names_filtered, factor_names)])
      diag(factor_cell_density_hm) <- 0
      factor_cell_density_hm <- as.data.frame(factor_cell_density_hm)
      factor_cell_density_hm <- factor_cell_density_hm[factor_names, new_c2l_names_filtered]
      rownames(factor_cell_density_hm) <- paste0(subset_id, ".F", gsub("factor_", "", rownames(factor_cell_density_hm)))
      return(factor_cell_density_hm)
    }), nm = animal_ids)
    return(factor_cell_density_hm)
  }
  factor_cell_density_hm <- cor(factor_cell_cor_subset[,c(new_c2l_names_filtered, factor_names)])
  diag(factor_cell_density_hm) <- 0
  factor_cell_density_hm <- as.data.frame(factor_cell_density_hm)
  factor_cell_density_hm <- factor_cell_density_hm[factor_names, new_c2l_names_filtered]
  rownames(factor_cell_density_hm) <- paste0(subset_id, ".F", gsub("factor_", "", rownames(factor_cell_density_hm)))
  return(factor_cell_density_hm)
}), 
nm = groups_use)

saveRDS(factor_cell_cor_data_group, file.path(DIR_OBJ_NMF, "mm_visium_nmf30_d21_factor_cell_cor_data_group.rds"))


# Plot pheatmap
pal_length <- col_scale_div_custom2 %>% length()
ph_colors <- col_scale_div_custom2
abs_max_val <- lapply(factor_cell_cor_data_group, function(m){max(abs(m))}) %>% unlist() %>% max()

pdf(file = file.path(DIR_FIG_NMF, paste0(file_save_prefix, "c2l_res_cell_factor_cor_heatmap_per_group.pdf")), 
    width = 10, height = 10, useDingbats = F)
for(g in groups_use){
  factor_cell_density_hm <- factor_cell_cor_data_group[[g]]
  
  abs_max_val <- max(abs(factor_cell_density_hm))
  ph_breaks <- c(seq(-abs_max_val, 0, length.out=ceiling(pal_length/2) + 1), 
                 seq(abs_max_val/pal_length, abs_max_val, length.out=floor(pal_length/2)))
  
  print(
    pheatmap::pheatmap(t(factor_cell_density_hm), 
                       cellwidth = 10, 
                       cellheight = 10, 
                       color = ph_colors, 
                       breaks = ph_breaks, 
                       main = g)
  )
}
dev.off()

###### Factor ~ Krt8ADI correlation ######
factor_cell_cor_data_group_krt8 <- do.call(cbind, factor_cell_cor_data_group) %>% select(contains("Krt8"))
factor_cell_cor_data_group_krt8

pal_length <- col_scale_div_custom2 %>% length()
ph_colors <- col_scale_div_custom2
abs_max_val <- max(factor_cell_cor_data_group_krt8)

pdf(file = file.path(DIR_FIG_NMF, paste0(file_save_prefix, "c2l_res_cell_factor_cor_heatmap_krt8ADI_group.pdf")), 
    width = 6, height = 6, useDingbats = F)
pheatmap::pheatmap(factor_cell_cor_data_group_krt8, 
                   cellwidth = 20, 
                   cellheight = 10, 
                   treeheight_row = 10,  
                   treeheight_col = 4,
                   color = ph_colors, 
                   breaks = ph_breaks)
dev.off()


###### Factor ~ Factor correlations ###### 
factor_cor_data <- bind_cols(se_subset@meta.data %>% select(sample_name, day, condition, group),
                             se_subset@reductions$NMF@cell.embeddings[, factor_names])

groups_use <- c("all", "CTRL", "d21_BLM")

factor_cor_data_group <- setNames(lapply(groups_use, function(g){
  if (g=="all") {
    factor_cor_subset <- factor_cor_data
  } else if (g=="CTRL") {
    factor_cor_subset <- subset(factor_cor_data, condition == "control")
  } else {
    factor_cor_subset <- subset(factor_cor_data, group == g)
  }
  factor_cor_hm <- cor(factor_cor_subset[, factor_names])
  diag(factor_cor_hm) <- 0
  factor_cor_hm <- as.data.frame(factor_cor_hm)
  rownames(factor_cor_hm) <- paste0(subset_id, ".F", gsub("factor_", "", rownames(factor_cor_hm)))
  colnames(factor_cor_hm) <- paste0(subset_id, ".F", gsub("factor_", "", colnames(factor_cor_hm)))
  return(factor_cor_hm)
}), 
nm = groups_use)


# Plot as pheatmaps
pal_length <- col_scale_div_custom2 %>% length()
ph_colors <- col_scale_div_custom2
abs_max_val <- lapply(factor_cell_cor_data_group, function(m){max(abs(m))}) %>% unlist() %>% max()

pdf(file = file.path(DIR_FIG_NMF, paste0(file_save_prefix, "factor_cor_heatmap_per_group.pdf")), width = 8, height = 8, useDingbats = F)
for(g in groups_use){
  factor_cor_hm <- factor_cor_data_group[[g]]
  
  abs_max_val <- max(abs(factor_cor_hm))
  ph_breaks <- c(seq(-abs_max_val, 0, length.out=ceiling(pal_length/2) + 1), 
                 seq(abs_max_val/pal_length, abs_max_val, length.out=floor(pal_length/2)))
  
  print(
    pheatmap::pheatmap(factor_cor_hm, 
                       cellwidth = 10, 
                       cellheight = 10, 
                       color = ph_colors, 
                       breaks = ph_breaks, 
                       main = g)
  )
}
dev.off()

###### Spatial plots ###### 

# Selected sample
se_b4a <- SubsetSTData(se_subset, sample_name %in% "d21_bleo_4a")
se_b5a <- SubsetSTData(se_subset, sample_name %in% "d21_bleo_5a")
se_b3 <- SubsetSTData(se_subset, sample_name %in% "d21_bleo_3")

se_b4a <- LoadImages(se_b4a, xdim = 1e3)
se_b5a <- LoadImages(se_b5a, xdim = 1e3)
se_b3 <- LoadImages(se_b3, xdim = 1e3)

# Selected factor vs cell2location densities - comparison
FeatureOverlay(se_b4a, features = c("factor_14", "c2l_Krt8.ADI"), label.by = "sample_name", 
               add.alpha = T, cols = rev(col_scale_spec), pt.size = 1.25) & 
  theme(aspect.ratio = 1, legend.position = "bottom")
FeatureOverlay(se_b4a, features = c("factor_18", "c2l_AT1.cells"), label.by = "sample_name", 
               add.alpha = T, cols = rev(col_scale_spec), pt.size = 1.25) & 
  theme(aspect.ratio = 1, legend.position = "bottom")
FeatureOverlay(se_b4a, features = c("factor_11", "c2l_Fibroblasts"), label.by = "sample_name", 
               add.alpha = T, cols = rev(col_scale_spec), pt.size = 1.25) & 
  theme(aspect.ratio = 1, legend.position = "bottom")
FeatureOverlay(se_b4a, features = c("factor_8", "c2l_Resolution.macrophages"), label.by = "sample_name", 
               add.alpha = T, cols = rev(col_scale_spec), pt.size = 1.25) & 
  theme(aspect.ratio = 1, legend.position = "bottom")

# Select pairs and export plots
FX_c2l_pair <- list(c("factor_21", "c2l_Myofibroblasts"),
                    c("factor_14", "c2l_Krt8.ADI"),
                    c("factor_18", "c2l_AT1.cells"),
                    c("factor_11", "c2l_Fibroblasts"),
                    c("factor_8", "c2l_Resolution.macrophages"))

pdf(file = file.path(DIR_FIG_NMF, paste0(file_save_prefix, "spatial_b4a_FX_vs_c2l.pdf")), 
    width = 8, height = 6, useDingbats = F)
for(pair in FX_c2l_pair){
  print(
    FeatureOverlay(se_b4a, features = pair, label.by = "sample_name",
                   add.alpha = T, cols = rev(col_scale_spec), pt.size = 1.25) & 
      theme(aspect.ratio = 1, legend.position = "bottom")
  )
}
dev.off()

pdf(file = file.path(DIR_FIG_NMF, paste0(file_save_prefix, "spatial_b5a_FX_vs_c2l.pdf")), 
    width = 8, height = 6, useDingbats = F)
for(pair in FX_c2l_pair){
  print(
    FeatureOverlay(se_b5a, features = pair, label.by = "sample_name",
                   add.alpha = T, cols = rev(col_scale_spec), pt.size = 1.25) & 
      theme(aspect.ratio = 1, legend.position = "bottom")
  )
}
dev.off()

pdf(file = file.path(DIR_FIG_NMF, paste0(file_save_prefix, "spatial_b3_FX_vs_c2l.pdf")), 
    width = 8, height = 6, useDingbats = F)
for(pair in FX_c2l_pair){
  print(
    FeatureOverlay(se_b3, features = pair, label.by = "sample_name",
                   add.alpha = T, cols = rev(col_scale_spec), pt.size = 1.25) & 
      theme(aspect.ratio = 1, legend.position = "bottom")
  )
}
dev.off()


##### Data: d7 ####
###### Subset ######
subset_id <- "d7"
se_subset <- SubsetSTData(se, day == "d7")
n_samples <- se_subset$sample_id %>% unique() %>% length()

###### Run NMF ######
#' Define params
n_factors <- 30  # new setting to match human NMF and to obtain more refined NMF results
DIR_FIG_NMF <- file.path(DIR_RES, "figures", paste0("NMF", n_factors, "_", subset_id))
DIR_OBJ_NMF <- file.path(DIR_RES, "objects", paste0("NMF", n_factors, "_", subset_id))
dir.create(DIR_OBJ_NMF); dir.create(DIR_FIG_NMF)
file_save_prefix <- paste0("mm_visium_nmf", n_factors, "_", subset_id, "_")
factor_names <- paste0("factor_", 1:n_factors)

#' Run NMF
se_subset <- RunNMF(se_subset, nfactors = n_factors)

# Save new se obj with new NMF
fname <- paste0("mm_visium_preproc_se_obj_subset_", subset_id, "_nmf30_c2l.rds")
saveRDS(se_subset, file = file.path(DIR_RES, "objects", fname))
se_subset <- readRDS(file = file.path(DIR_RES, "objects", fname))


###### NMF Results ######
#' Look at factors that contains the Krt8+ADI marker genes most highly
adi_genes <- c(grep("^Krt", rownames(se_subset@reductions$NMF@feature.loadings), value = T), "Hbegf", "Plaur", "Areg", "Sprr1a", "Tnip3", "Ckn1", "Itgb6", "Cldn4", "End1")
adi_genes <- intersect(adi_genes, rownames(se_subset@reductions$NMF@feature.loadings)); adi_genes
factors_adi <- names(sort(colSums(se_subset@reductions$NMF@feature.loadings[adi_genes, ]), decreasing = T))[1:3]; factors_adi

# FactorPlot(se_subset, factor = 14, indices = 1:6, n.columns = 6, col.scale = col_scale_mako)
f_plot <- gsub("factor_", "", factors_adi) %>% as.numeric()
f_plot_list <- lapply(f_plot, function(f){
  p <- FactorPlot(se_subset, factor = f, indices = 1:12, n.columns = 6, col.scale = col_scale_mako)
})

png(file = file.path(DIR_FIG_NMF, paste0("mm_visium_nmf_potential_Krt8ADI_factors_spatial_gene.png")),
    width = 18*fig_res, height = 16*fig_res, res = fig_res)
wrap_plots(f_plot_list, ncol = 1)
dev.off()

plot_list <- lapply(f_plot, function(factor_plot) {
  p1 <- ST.DimPlot(se_subset, 
                   indices = 1:12,
                   dims = factor_plot, 
                   min.cutoff = 0.2,
                   label.by = "sample_name",
                   center.zero = F, 
                   reduction = "NMF", 
                   ncol = 6, 
                   pt.size = 0.8,
                   pt.border = F,
                   cols = col_scale_mako,
                   show.sb = F)  &
    labs(fill = "factor\nactivity", title = paste0("Factor ", factor_plot)) & 
    theme(legend.position = "right", plot.title = element_text(face = "bold", size = 16), aspect.ratio = 1)
  p2 <- FactorGeneLoadingPlotHorizontal(se_subset, 
                                        factor = factor_plot, 
                                        topn = 30) + 
    labs(title = "Genes contributing to factor") & 
    theme(axis.text.x = element_text(size=10), panel.grid = element_blank(),
          plot.margin = margin(10, 10, 0, 10, unit = "pt"))
  
  p3 <- FactorRankSpotPlot(se_subset, factor = factor_plot, 
                           split.by = "sample_name", group.by = "group", 
                           color.group = cols_group, plot.title = "Spots ranked by factor activity")
  
  p_bottom <- cowplot::plot_grid(p2, p3, nrow = 1)
  p <- cowplot::plot_grid(p1, p_bottom, ncol = 1, rel_heights = c(3, 2))
})

png(file = file.path(DIR_FIG_NMF, paste0("mm_visium_nmf_potential_Krt8ADI_factors_grid.png")),
    width = 12*fig_res, height = 7*3*fig_res, res = fig_res)
wrap_plots(plot_list, ncol = 1)
dev.off()

#' Look for collgen containing factors
fib_genes <- c(grep("^Col", rownames(se_subset@reductions$NMF@feature.loadings), value = T), "Fn1")
fib_genes <- intersect(fib_genes, rownames(se_subset@reductions$NMF@feature.loadings)); fib_genes
factors_fib <- names(sort(colSums(se_subset@reductions$NMF@feature.loadings[fib_genes, ]), decreasing = T))[1:3]; factors_fib


#' Save top genes for each factor
# se_subset@reductions$NMF@feature.loadings
factor_gene_loadings <- lapply(1:n_factors, function(factor_x){
  message(paste("Factor", factor_x))
  feat_loads <- as.data.frame(se_subset@reductions$NMF@feature.loadings[,paste0("factor_",factor_x)])
  colnames(feat_loads) <- "gene_loading"
  feat_loads$gene <- rownames(feat_loads)
  feat_loads$factor <- factor_x
  top_n_genes <- 100
  feat_loads <- feat_loads %>%
    dplyr::slice_max(order_by = gene_loading, n = top_n_genes) %>%
    mutate(rank = dense_rank(desc(gene_loading)))
})

write_xlsx(
  factor_gene_loadings,
  path = file.path(DIR_OBJ_NMF, paste0(file_save_prefix, "top100_gene_loadings.xlsx")),
  col_names = TRUE,
  format_headers = TRUE,
  use_zip64 = FALSE
)
factor_gene_loadings <- lapply(1:n_factors, function(factor_x){
  read_xlsx(path = file.path(DIR_OBJ_NMF, paste0(file_save_prefix, "top100_gene_loadings.xlsx")), sheet = factor_x)
})

###### gProfiler #####
org <- "mmusculus"
fname <- paste0(file_save_prefix, "gProfiler_top25genes")
n_genes <- 25
gene_query_list <- list()
gostres_list <- list()

for(i in 1:n_factors){
  gene_query_list[[i]] <- subset(factor_gene_loadings[[i]], rank <= n_genes)$gene
}
names(gene_query_list) <- factor_names

for (f in names(gene_query_list)){
  message(f)
  gostres_list[[f]] <- gprofiler2::gost(query = gene_query_list[[f]],
                                        organism = org, 
                                        ordered_query = T)
}

gostres_export_list <- list()
for (f in names(gostres_list)){
  gostres_export_list[[f]] <- gostres_list[[f]]$result
  gostres_export_list[[f]]$query <- f
}
write_xlsx(
  x = gostres_export_list,
  path = file.path(DIR_OBJ_NMF, paste0(fname, ".xlsx")),
  col_names = TRUE,
  format_headers = TRUE,
  use_zip64 = FALSE
)

pdf(file = file.path(DIR_FIG_NMF, paste0(fname, ".pdf")), width = 10, height = 8, useDingbats = F)
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


###### Factor plot grid ####
plot_list <- lapply(1:n_factors, function(factor_plot) {
  p1 <- ST.DimPlot(se_subset, 
                   indices = c(1:9),
                   dims = factor_plot, 
                   min.cutoff = 0.2,
                   label.by = "sample_name",
                   center.zero = F, 
                   reduction = "NMF", 
                   ncol = 6, 
                   pt.size = 0.8, 
                   pt.border = F,
                   cols = rev(RColorBrewer::brewer.pal(10, "Spectral")),
                   show.sb = F)  &
    labs(fill = "factor\nactivity", title = paste0("Factor ", factor_plot)) & 
    theme(aspect.ratio = 1, legend.position = "right", plot.title = element_text(face = "bold", size = 16))
  p2 <- FactorGeneLoadingPlotHorizontal(se_subset, 
                                        factor = factor_plot, 
                                        topn = 30) + 
    labs(title = "Genes contributing to factor") & 
    theme(axis.text.x = element_text(size=10), panel.grid = element_blank(),
          plot.margin = margin(10, 10, 0, 10, unit = "pt"))
  
  p3 <- FactorRankSpotPlot(se_subset, factor = factor_plot, split.by = "sample_name", 
                           group.by = "group", color.group = cols_group, plot.title = "Spots ranked by factor activity")
  
  p_bottom <- cowplot::plot_grid(p2, p3, nrow = 1)
  p <- cowplot::plot_grid(p1, p_bottom, ncol = 1, rel_heights = c(3, 2))
})

dir.create(path = file.path(DIR_FIG_NMF, paste0(file_save_prefix, "grid")))
for(i in 1:length(plot_list)){
  png(file = file.path(DIR_FIG_NMF, paste0(file_save_prefix, "grid"), paste0("mm_visium_nmf_grid_f", i, ".png")),
      width = 12*fig_res, height = 7*fig_res, res = fig_res)
  print(plot_list[[i]])
  dev.off()
}


###### Factor ~ Cell correlation ######
# c2l_all <- read.csv(file.path(DIR_ROOT, "data", SPECIES, "sc_deconvolution_strunz", "compiled_all_samples_cell_abundances.csv"), row.names = 1)
# se_subset <- AddMetaData(se_subset, c2l_all)
c2l_names <- grep("c2l_", colnames(se_subset@meta.data), value = T)
new_c2l_names <- gsub("c2l_", "", gsub("[.]$", "", gsub("..", "_", c2l_names, fixed = TRUE)))

factor_cell_cor_data <- bind_cols(se_subset@meta.data %>% select(sample_name, day, condition, animal, group, matches(c2l_names)),
                                  se_subset@reductions$NMF@cell.embeddings[, factor_names])
colnames(factor_cell_cor_data)[colnames(factor_cell_cor_data) %in% c2l_names] <- new_c2l_names


# Filter cells to show
new_c2l_names_filtered <- new_c2l_names[!new_c2l_names %in% c("NA", "low.quality.cells")]

animal_ids <- se_subset$animal %>% unique()
groups_use <- c("all", "CTRL", "d7_BLM", "animal_id")
# groups_use <- c("all", "CTRL", "d7_BLM")
factor_cell_cor_data_group <- setNames(lapply(groups_use, function(g){
  if (g=="all") {
    factor_cell_cor_subset <- factor_cell_cor_data
  } else if (g=="CTRL") {
    factor_cell_cor_subset <- subset(factor_cell_cor_data, condition == "control")
    # } else {
    #   factor_cell_cor_subset <- subset(factor_cell_cor_data, group == g)
    # }
  } else if (g == "d7_BLM") {
    factor_cell_cor_subset <- subset(factor_cell_cor_data, group == g)
  } else {
    factor_cell_density_hm <- setNames(lapply(animal_ids, function(a){
      message(a)
      factor_cell_cor_subset <- subset(factor_cell_cor_data, animal == a)
      factor_cell_density_hm <- cor(factor_cell_cor_subset[,c(new_c2l_names_filtered, factor_names)])
      diag(factor_cell_density_hm) <- 0
      factor_cell_density_hm <- as.data.frame(factor_cell_density_hm)
      factor_cell_density_hm <- factor_cell_density_hm[factor_names, new_c2l_names_filtered]
      rownames(factor_cell_density_hm) <- paste0(subset_id, ".F", gsub("factor_", "", rownames(factor_cell_density_hm)))
      return(factor_cell_density_hm)
    }), nm = animal_ids)
    return(factor_cell_density_hm)
  }
  factor_cell_density_hm <- cor(factor_cell_cor_subset[,c(new_c2l_names_filtered, factor_names)])
  diag(factor_cell_density_hm) <- 0
  factor_cell_density_hm <- as.data.frame(factor_cell_density_hm)
  factor_cell_density_hm <- factor_cell_density_hm[factor_names, new_c2l_names_filtered]
  rownames(factor_cell_density_hm) <- paste0(subset_id, ".F", gsub("factor_", "", rownames(factor_cell_density_hm)))
  return(factor_cell_density_hm)
}), 
nm = groups_use)

saveRDS(factor_cell_cor_data_group, file.path(DIR_OBJ_NMF, "mm_visium_nmf30_d7_factor_cell_cor_data_group.rds"))


# Plot pheatmap
pal_length <- col_scale_div_custom2 %>% length()
ph_colors <- col_scale_div_custom2
abs_max_val <- lapply(factor_cell_cor_data_group, function(m){max(abs(m))}) %>% unlist() %>% max()

pdf(file = file.path(DIR_FIG_NMF, paste0(file_save_prefix, "c2l_res_cell_factor_cor_heatmap_per_group.pdf")), 
    width = 10, height = 10, useDingbats = F)
for(g in groups_use[1:3]){
  factor_cell_density_hm <- factor_cell_cor_data_group[[g]]
  
  abs_max_val <- max(abs(factor_cell_density_hm))
  ph_breaks <- c(seq(-abs_max_val, 0, length.out=ceiling(pal_length/2) + 1), 
                 seq(abs_max_val/pal_length, abs_max_val, length.out=floor(pal_length/2)))
  
  print(
    pheatmap::pheatmap(t(factor_cell_density_hm), 
                       cellwidth = 10, 
                       cellheight = 10, 
                       color = ph_colors, 
                       breaks = ph_breaks, 
                       main = g)
  )
}
dev.off()

###### Factor ~ Krt8ADI correlation ######
factor_cell_cor_data_group_krt8 <- do.call(cbind, factor_cell_cor_data_group[1:3]) %>% select(contains("Krt8"))
factor_cell_cor_data_group_krt8

pal_length <- col_scale_div_custom2 %>% length()
ph_colors <- col_scale_div_custom2
abs_max_val <- max(factor_cell_cor_data_group_krt8)

pdf(file = file.path(DIR_FIG_NMF, paste0(file_save_prefix, "c2l_res_cell_factor_cor_heatmap_krt8ADI_group.pdf")), 
    width = 6, height = 6, useDingbats = F)
pheatmap::pheatmap(factor_cell_cor_data_group_krt8, 
                   cellwidth = 20, 
                   cellheight = 10, 
                   treeheight_row = 10,  
                   treeheight_col = 4,
                   color = ph_colors, 
                   breaks = ph_breaks)
dev.off()


###### Factor ~ Factor correlations ###### 
factor_cor_data <- bind_cols(se_subset@meta.data %>% select(sample_name, day, condition, group),
                             se_subset@reductions$NMF@cell.embeddings[, factor_names])

groups_use <- c("all", "CTRL", "d7_BLM")

factor_cor_data_group <- setNames(lapply(groups_use, function(g){
  if (g=="all") {
    factor_cor_subset <- factor_cor_data
  } else if (g=="CTRL") {
    factor_cor_subset <- subset(factor_cor_data, condition == "control")
  } else {
    factor_cor_subset <- subset(factor_cor_data, group == g)
  }
  factor_cor_hm <- cor(factor_cor_subset[, factor_names])
  diag(factor_cor_hm) <- 0
  factor_cor_hm <- as.data.frame(factor_cor_hm)
  rownames(factor_cor_hm) <- paste0(subset_id, ".F", gsub("factor_", "", rownames(factor_cor_hm)))
  colnames(factor_cor_hm) <- paste0(subset_id, ".F", gsub("factor_", "", colnames(factor_cor_hm)))
  return(factor_cor_hm)
}), 
nm = groups_use)


# Plot as pheatmaps
pal_length <- col_scale_div_custom2 %>% length()
ph_colors <- col_scale_div_custom2
abs_max_val <- lapply(factor_cell_cor_data_group, function(m){max(abs(m))}) %>% unlist() %>% max()

pdf(file = file.path(DIR_FIG_NMF, paste0(file_save_prefix, "factor_cor_heatmap_per_group.pdf")), width = 8, height = 8, useDingbats = F)
for(g in groups_use){
  factor_cor_hm <- factor_cor_data_group[[g]]
  
  abs_max_val <- max(abs(factor_cor_hm))
  ph_breaks <- c(seq(-abs_max_val, 0, length.out=ceiling(pal_length/2) + 1), 
                 seq(abs_max_val/pal_length, abs_max_val, length.out=floor(pal_length/2)))
  
  print(
    pheatmap::pheatmap(factor_cor_hm, 
                       cellwidth = 10, 
                       cellheight = 10, 
                       color = ph_colors, 
                       breaks = ph_breaks, 
                       main = g)
  )
}
dev.off()


###### Spatial plots ###### 

# Selected sample
se_b1 <- SubsetSTData(se_subset, sample_name %in% "d7_bleo_1")
se_b3 <- SubsetSTData(se_subset, sample_name %in% "d7_bleo_3")

se_b1 <- LoadImages(se_b1, xdim = 1e3)
se_b3 <- LoadImages(se_b3, xdim = 1e3)


# Selected factor vs cell2location densities - comparison
FeatureOverlay(se_b1, features = c("factor_12", "c2l_Myofibroblasts"), label.by = "sample_name", 
               add.alpha = T, cols = rev(col_scale_spec), pt.size = 1.25) & 
  theme(aspect.ratio = 1, legend.position = "bottom")
FeatureOverlay(se_b1, features = c("factor_16", "c2l_Krt8.ADI"), label.by = "sample_name", 
               add.alpha = T, cols = rev(col_scale_spec), pt.size = 1.25) & 
  theme(aspect.ratio = 1, legend.position = "bottom")
FeatureOverlay(se_b1, features = c("factor_5", "c2l_AT1.cells"), label.by = "sample_name", 
               add.alpha = T, cols = rev(col_scale_spec), pt.size = 1.25) & 
  theme(aspect.ratio = 1, legend.position = "bottom")
FeatureOverlay(se_b1, features = c("factor_4", "c2l_Fibroblasts"), label.by = "sample_name", 
               add.alpha = T, cols = rev(col_scale_spec), pt.size = 1.25) & 
  theme(aspect.ratio = 1, legend.position = "bottom")
FeatureOverlay(se_b1, features = c("factor_15", "c2l_Resolution.macrophages"), label.by = "sample_name", 
               add.alpha = T, cols = rev(col_scale_spec), pt.size = 1.25) & 
  theme(aspect.ratio = 1, legend.position = "bottom")
FeatureOverlay(se_b1, features = c("factor_11", "c2l_T.lymphocytes"), label.by = "sample_name", 
               add.alpha = T, cols = rev(col_scale_spec), pt.size = 1.25) & 
  theme(aspect.ratio = 1, legend.position = "bottom")

FeatureOverlay(se_b3, features = c("c2l_VECs", "c2l_Krt8.ADI"), label.by = "sample_name", 
               add.alpha = T, cols = rev(col_scale_spec), pt.size = 1.25) & 
  theme(aspect.ratio = 1, legend.position = "bottom")

# Select pairs and export plots
FX_c2l_pair <- list(c("factor_12", "c2l_Myofibroblasts"),
                    c("factor_16", "c2l_Krt8.ADI"),
                    c("factor_5", "c2l_AT1.cells"),
                    c("factor_4", "c2l_Fibroblasts"),
                    c("factor_15", "c2l_Resolution.macrophages"),
                    c("factor_11", "c2l_T.lymphocytes"),
                    c("c2l_VECs", "c2l_Krt8.ADI"))

pdf(file = file.path(DIR_FIG_NMF, paste0(file_save_prefix, "spatial_b1_FX_vs_c2l.pdf")), 
    width = 8, height = 6, useDingbats = F)
for(pair in FX_c2l_pair){
  print(
    FeatureOverlay(se_b1, features = pair, label.by = "sample_name",
                   add.alpha = T, cols = rev(col_scale_spec), pt.size = 1.25) & 
      theme(aspect.ratio = 1, legend.position = "bottom")
  )
}
dev.off()


pdf(file = file.path(DIR_FIG_NMF, paste0(file_save_prefix, "spatial_b3_FX_vs_c2l.pdf")), 
    width = 8, height = 6, useDingbats = F)
for(pair in FX_c2l_pair){
  print(
    FeatureOverlay(se_b3, features = pair, label.by = "sample_name",
                   add.alpha = T, cols = rev(col_scale_spec), pt.size = 1.25) & 
      theme(aspect.ratio = 1, legend.position = "bottom")
  )
}
dev.off()


##### d21 NMF30: AbBa factor (F14) analysis ####
#' Further analysis of spots high for AbBa-associated factor
#' in the d21 subset NMF30 analysis
#' 
# Read data
subset_id <- "d21"
n_factors <- 30
DIR_FIG_NMF <- file.path(DIR_RES, "figures", paste0("NMF", n_factors, "_", subset_id))
DIR_OBJ_NMF <- file.path(DIR_RES, "objects", paste0("NMF", n_factors, "_", subset_id))
file_save_prefix <- paste0("mm_visium_nmf", n_factors, "_", subset_id, "_")
factor_names <- paste0("factor_", 1:n_factors)

fname <- paste0("mm_visium_preproc_se_obj_subset_", subset_id, "_nmf30_c2l.rds")
se_subset <- readRDS(file = file.path(DIR_RES, "objects", fname))

# Select factor
adi_genes <- c(grep("^Krt", rownames(se_subset@reductions$NMF@feature.loadings), value = T), "Hbegf", "Plaur", "Areg", "Sprr1a", "Tnip3", "Ckn1", "Itgb6", "Cldn4", "End1")
adi_genes <- intersect(adi_genes, rownames(se_subset@reductions$NMF@feature.loadings)); adi_genes
factors_adi <- names(sort(colSums(se_subset@reductions$NMF@feature.loadings[adi_genes, ]), decreasing = T))[1:3]; factors_adi

# Cut-off: Select top 1% of spots
#' Label spots top 1% most factor enriched spots of the ADI/fibrotic factors
nmf_gene_loads <- se_subset@reductions$NMF@feature.loadings
nmf_gene_loads_scaled <- apply(nmf_gene_loads, 2, Scale01)
se_nmf_emb <- se_subset@reductions$NMF@cell.embeddings
# factor_xy <- names(sort(rowSums(nmf_gene_loads_scaled), decreasing = T))[1:4]
factor_xy <- factors_adi
se_subset <- AddMetaData(se_subset, metadata = se_nmf_emb[,factor_xy], col.name = factor_xy)

factor_cutoff <- setNames(
  lapply(factor_xy, function(x){
    as.numeric(quantile(se_nmf_emb[,x], c(.99)))
  }), 
  nm = factor_xy)
f_metadata_add <- se_subset@meta.data[,1:2]
for(f in factor_xy){
  f_metadata_add[[f]] <- paste0(f, "_low")
  spots_f_cutoff <- intersect(colnames(se_subset), rownames(se_nmf_emb[se_nmf_emb[,f]>factor_cutoff[[f]],]))
  f_metadata_add[spots_f_cutoff,f] <- paste0(f, "_high")
}
se_subset <- AddMetaData(se_subset, metadata = f_metadata_add[,factor_xy], col.name = paste0(factor_xy, "_cutoff"))


# Look at selection
summary(as.factor(se_subset$factor_14_cutoff))

p1 <- ggplot(se_subset@meta.data, aes(x="", y=factor_14)) +
  geom_boxplot(fill="grey", outlier.size = 0.5) +
  geom_hline(yintercept = factor_cutoff$factor_14, color="orange") +
  theme_bw() + theme(axis.title.x = element_blank(), panel.grid = element_blank())

p2 <- ggplot(se_subset@meta.data, aes(animal, fill = factor_14_cutoff)) +
  geom_histogram(stat = "count") +
  scale_fill_manual(values = c("orange", "grey")) +
  scale_y_log10() +
  theme_bw() + theme(panel.grid = element_blank())

p3 <- ST.FeaturePlot(
  se_subset, 
  indices = 1:12,
  features = "factor_14_cutoff", 
  label.by = "sample_name",
  ncol = 6, 
  pt.size = 0.8,
  pt.border = F,
  cols = c("grey90", "darkorange"),
  show.sb = F) +
  theme(aspect.ratio = 1, legend.position = "right")

pdf(file = file.path(DIR_FIG_NMF, "mm_visium_nmf30_d21_f14_cutoff_plots.pdf"), width = 12, height = 8, useDingbats = F)
(((p1|p2) + plot_layout(widths = c(1,3))) / p3) + plot_layout(heights = c(1,2))
dev.off()


###### Subcluster d21-NMF30-F14-high spots #####
spots_f14 <- rownames(se_subset@meta.data[se_subset$factor_14_cutoff %in% "factor_14_high",])
se_f14high <- SubsetSTData(se_subset, spots = spots_f14)

se_f14high <- RunPCA(se_f14high)
ElbowPlot(se_f14high)
DimHeatmap(se_f14high, reduction = "pca", dims = 1:14)
VlnPlot(se_f14high, features = paste0("PC_",1:14), group.by = "animal", ncol=7, pt.size = 0)

#' PCA only:
dims_select <- 1:14
se_f14high <- RunUMAP(se_f14high, reduction = "pca", dims = dims_select)
se_f14high <- FindNeighbors(se_f14high, reduction = "pca", dims = dims_select)
se_f14high <- FindClusters(se_f14high, resolution = 0.5)

#' Save clusters
se_f14high$f14_subclusters <- se_f14high$seurat_clusters
n_clusters <- length(unique(se_f14high$f14_subclusters))
cols_clusters <- setNames(object = ggsci::pal_tron()(n_clusters), #c(brewer.pal(n_clusters, "Set1")), 
                          nm = levels(se_f14high$f14_subclusters))
DimPlot(se_f14high, reduction = "umap", group.by = "f14_subclusters", cols = cols_clusters)

# Plot umap
# DimPlot(se_f14high, reduction = "umap", group.by = "animal")
# DimPlot(se_f14high, reduction = "umap", group.by = "f14_subclusters", split.by = "animal")

# Marker genes
se_f14high <- SetIdent(se_f14high, value = "f14_subclusters")
makers_f14_subcluster <- FindAllMarkers(se_f14high)

#'export
fname <- paste0("mm_visium_nmf30_d21_f14_subclusters_markers.csv")
write.csv(x = makers_f14_subcluster, file = file.path(DIR_OBJ_NMF, fname), row.names = F)

saveRDS(object = se_f14high, file = file.path(DIR_OBJ_NMF, "mm_visium_nmf30_d21_f14high_subset.rds"))
se_f14high <- readRDS(file.path(DIR_OBJ_NMF, "mm_visium_nmf30_d21_f14high_subset.rds"))


# Plot marker heatmap 
genes_plot <- subset(makers_f14_subcluster, p_val_adj<0.05)
genes_plot$updown <- "up"
genes_plot[genes_plot$avg_log2FC<0,"updown"] <- "down"
genes_plot <- genes_plot %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group=TRUE) %>%
  top_n(30, (avg_log2FC))

mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(256)
p <- DoHeatmap(se_f14high, features = c(genes_plot$gene), group.by = "f14_subclusters", group.colors =  cols_clusters, slot = "scale.data") +
  scale_fill_gradientn(colours = rev(mapal));p

png(file.path(DIR_FIG_NMF, paste0("mm_visium_nmf30_d21_f14hi_subclusters_marker_heatmap.png")), width = 6*fig_res, height = 10*fig_res, res = fig_res);p;dev.off()

# Add in metadata from 'se_f14high' to 'se_subset'
mdata_add <- setNames(paste0("F14hi_C", se_f14high$f14_subclusters), nm=names(se_f14high$f14_subclusters))
se_subset <- AddMetaData(se_subset, mdata_add, "f14_subclusters")

se_subset$f14_subclusters[is.na(se_subset$f14_subclusters)] <- "other"
f14hi_cluster_names <- se_subset$f14_subclusters[se_subset$f14_subclusters != "other"] %>% unique() %>% sort()

# PLOT
#' Fig 5b::
cluster_cols_all <- setNames(
  c("#DB8712", "#49C1AD", "#8D71A6",  # "#E7B922"
             # ggsci::pal_tron()(n_clusters),  # n_clusters=3
             "grey95"), 
             nm = c(f14hi_cluster_names, "other"))
cluster_cols_all2 <- c(cluster_cols_all, setNames("grey85", nm= "fibrosis"))

#' Spatial
se_subset$f14_subclusters2 <- ifelse(se_subset$f14_subclusters != "other", 
                                     as.character(se_subset$f14_subclusters),
                                     ifelse(se_subset$annotation %in% "Suspect Fibrosis/Fibroplasia", 
                                            "fibrosis",
                                            "other"
                                     ))
se_subset$f14_subclusters2 %>% unique()

p_spat <- ST.FeaturePlot(se_subset,
                         features = "f14_subclusters2", 
                         ncol=6, cols = cluster_cols_all2, 
                         label.by = "sample_name", 
                         pt.size = 0.5, show.sb = F) + 
  guides(fill = guide_legend(override.aes = list(size=2))) &
  theme(legend.position = "right", aspect.ratio = 1);p_spat


p_spat2 <- ST.FeaturePlot(SubsetSTData(se_subset, sample_name %in% paste0("d21_bleo_", c("1","2","3","4a","6a","5b"))),
                          features = "f14_subclusters2", 
                          ncol=3, cols = cluster_cols_all2, 
                          label.by = "animal", 
                          pt.size = 0.5, show.sb = F) + 
  guides(fill = guide_legend(override.aes = list(size=2))) &
  theme(legend.position = "right", aspect.ratio = 1);p_spat2


#' Proportions
d <- se_subset@meta.data %>%
  filter(f14_subclusters != "other") %>%
  group_by(sample_name, animal, f14_subclusters) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

p_prop <- ggplot(d, aes(x=sample_name, y=n, fill=f14_subclusters)) +
  geom_bar(stat = 'identity', colour=NA, position = "stack") +
  scale_fill_manual(values = cluster_cols_all) +
  scale_y_log10() +
  labs(x="", y="n spots (log10 scale)", fill="") +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        panel.grid = element_blank(),
        legend.position = "right", 
        text = element_text(size=10),
        plot.title = element_text(hjust=0.5));p_prop

#' Marker genes dotplot
genes_plot2 <- makers_f14_subcluster %>%
  filter(avg_log2FC>0, p_val_adj<0.01, (pct.1-pct.2)>0) %>% 
  group_by(cluster) %>%
  slice_max(n = 15, order_by = avg_log2FC)

p_dp <- DotPlot(se_subset,
                features = unique(genes_plot2$gene), 
                group.by = "f14_subclusters",
                scale = T) +
  scale_color_gradientn(colours = c(rev(col_scale_mako[1:5]), col_scale_acton[1:8])) +
  theme_dotplot +
  theme(axis.text.y = element_text(face = "plain"),
        axis.title = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, face="italic"),
        legend.position = "right");p_dp


#' Marker genes dotplot - reduced
genes_plot3 <- makers_f14_subcluster %>%
  filter(avg_log2FC>0, p_val_adj<0.01, (pct.1-pct.2)>0.1) %>% 
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

p_dp2 <- DotPlot(subset(se_subset, f14_subclusters != "other"),
                 features = unique(genes_plot3$gene) %>% rev(), 
                 group.by = "f14_subclusters",
                 scale = T) +
  scale_color_gradientn(colours = c(rev(col_scale_mako[1:5]), col_scale_acton[1:8])) +
  scale_x_discrete(position = "top") +
  coord_flip() +
  theme_dotplot +
  theme(
    # axis.text.y = element_text(face = "plain"),
    axis.title = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1, vjust=1), # , face="italic"
    legend.position = "right");p_dp2s

#' Cell density violin plots
png(file.path(DIR_FIG_NMF, paste0("mm_visium_nmf30_d21_f14hi_subclusters_c2l_violins.png")), width = 16*fig_res, height = 8*fig_res, res = fig_res)
VlnPlot(se_subset, 
        features = c2l_names,
        ncol = 10,
        pt.size = 0, 
        group.by = "f14_subclusters", cols = cluster_cols_all) &
  theme_dotplot &
  theme(legend.position = "none", 
        axis.text.y = element_text(face = "plain"))
dev.off()

c2l_names_plot <- c(grep("AT|Krt", c2l_names, value = T),
                    grep("ibrob", c2l_names, value = T))

p_vln <- VlnPlot(se_subset, 
                 features = c2l_names_plot,
                 ncol = 3,
                 pt.size = 0, 
                 group.by = "f14_subclusters", cols = cluster_cols_all) &
  theme_dotplot &
  theme(legend.position = "none", 
        axis.text.y = element_text(face = "plain"));p_vln

p_vln2 <- VlnPlot(se_subset, 
                  features = c2l_names_plot,
                  ncol = 3,
                  pt.size = 0, 
                  group.by = "f14_subclusters", cols = cluster_cols_all) &
  theme_dotplot & 
  theme(legend.position = "none", 
        axis.text.x = element_blank(),
        axis.text.y = element_text(face = "plain")
  )

se_subset$f14_subclusters_reord <- factor(se_subset$f14_subclusters,
                                          levels = sort(unique(se_subset$f14_subclusters)) %>% rev())
p_vln3 <- VlnPlot(se_subset, 
                  features = c2l_names_plot,
                  ncol = 1, y.max = 6,
                  pt.size = 0, 
                  group.by = "f14_subclusters_reord", cols = cluster_cols_all) &
  coord_flip() &
  theme_dotplot & 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle=0),
        axis.text.y = element_text(face = "plain")
  );p_vln3
pdf(file = file.path(DIR_FIG_NMF, "mm_visium_nmf30_d21_f14hi_subclusters_c2l_plots_violin.pdf"), 
    width = 3, height = 8, useDingbats = F)
p_vln3
dev.off()

# Combine plot grid
p_r <- (p_spat / p_dp) + plot_layout(heights = c(4,1))
p_l <- (p_prop / p_vln)

pdf(file = file.path(DIR_FIG_NMF, "mm_visium_nmf30_d21_f14hi_subclusters_c2l_plots.pdf"), 
    width = 18, height = 9, useDingbats = F)
(p_l|p_r)+plot_layout(widths = c(2,5))
dev.off()

# Combine plot grid-2
pdf(file = file.path(DIR_FIG_NMF, "mm_visium_nmf30_d21_f14hi_subclusters_c2l_plots2.pdf"), 
    width = 7+6.5, height = 4, useDingbats = F)
(p_vln2|plot_spacer()|p_dp2|p_spat2)+plot_layout(widths = c(5,0.25,1,7))
dev.off()

###### Select d21-NMF30-F14-high neighbors ###### 
# VlnPlot(se_subset, features = "S100a6", group.by = "f14_subclusters")

# Neighbouring spots
se_subset <- LoadImages(se_subset, xdim=50)
se_subset$f14_subclusters <- as.character(se_subset$f14_subclusters)
se_subset <- SetIdent(se_subset, value ="f14_subclusters")
# se_subset$d_cF14hi_C0 <- NULL
se_subset <- IdentifyNNeighbors(se_subset, feature.column.name = "f14_subclusters", center.feature = "F14hi_C0", n.neighbors = 4)

# se_subset <- RegionNeighbours(se_subset, id = "F14hi_C0")


unique(se_subset$d_cF14hi_C0)
cols_dist <- c("1"=col_scale_spec[2], "2"=col_scale_spec[4])

se_nbs <- SubsetSTData(se_subset, d_cF14hi_C0 == 1) #nbs_F14hi_C0 %in% "nbs_F14hi_C0"
se_nbs <- RunPCA(se_nbs, npcs = 20)

# ElbowPlot(se_nbs)
# DimHeatmap(se_nbs, reduction = "pca", dims = 1:16)
# VlnPlot(se_nbs, features = paste0("PC_", 1:15), group.by = "animal", ncol=5, pt.size = 0)

dims_select <- 1:10 # nbs 1-2
se_nbs <- RunUMAP(se_nbs, reduction = "pca", dims = dims_select)
se_nbs <- FindNeighbors(se_nbs, reduction = "pca", dims = dims_select)
se_nbs <- FindClusters(se_nbs, resolution = 0.2)
n_clusters_nbs <- length(unique(se_nbs$seurat_clusters))
cluster_cols_nbs <- setNames(object = RColorBrewer::brewer.pal(n_clusters_nbs+2, "Spectral")[-c(1,4)], 
                             nm = sort(unique(se_nbs$seurat_clusters)))

p1 <- DimPlot(se_nbs, reduction = "umap", group.by = "seurat_clusters", pt.size=0.5,  cols = cluster_cols_nbs, label = T, label.box = T, repel = T) + 
  labs(title="F14-C0 nbs clusters") +
  theme_umap
p2 <- DimPlot(se_nbs, reduction = "umap", group.by = "d_cF14hi_C0", pt.size=0.5, cols = cols_dist) + 
  labs(title="Distance to F14-C0") +
  theme_umap
p_nbs_umap <- (p1/p2) & theme(aspect.ratio = 1);p_nbs_umap

# Marker genes
se_nbs <- SetIdent(se_nbs, value = "seurat_clusters")
makers_f14_subcluster_nbs <- FindAllMarkers(se_nbs)

#'export
fname <- paste0("mm_visium_nmf30_d21_f14hi_subclusters_nbs_cluster_markers.csv")
write.csv(x = makers_f14_subcluster_nbs, file = file.path(DIR_OBJ_NMF, fname), row.names = F)
makers_f14_subcluster_nbs <- read.csv(file = file.path(DIR_OBJ_NMF, fname))

# Check top genes
makers_f14_subcluster_nbs %>% 
  filter(avg_log2FC>0, p_val_adj<0.01, (pct.1-pct.2)>0) %>% 
  group_by(cluster) %>% 
  slice_max(n = 10, order_by = avg_log2FC) %>% 
  print(n=30)


# Add nbs cluster metadata to se_subset object
se_nbs$f14_nbs_clusters <- paste0("F14hi_nbs_C", se_nbs$seurat_clusters)
summary(as.factor(se_nbs$f14_nbs_clusters));summary(as.factor(se_nbs$d_cF14hi_C0))

se_subset <- AddMetaData(se_subset, se_nbs$f14_nbs_clusters, "f14_nbs_clusters")

se_subset$f14_nbs_clusters <- as.character(se_subset$f14_nbs_clusters)
se_subset$f14_nbs_clusters[is.na(se_subset$f14_nbs_clusters)] <- "other"
se_subset$f14_nbs_clusters[se_subset$f14_subclusters=="F14hi_C0"] <- "F14hi_C0"
se_subset$f14_nbs_clusters <- factor(se_subset$f14_nbs_clusters, levels = c(unique(se_nbs$f14_nbs_clusters) %>% sort(), "F14hi_C0", "other"))
summary(as.factor(se_subset$f14_nbs_clusters));summary(as.factor(se_subset$d_cF14hi_C0))

cluster_cols_nbs_all <- setNames(object = c(cluster_cols_nbs, "grey10", "grey90"), nm = levels(as.factor(se_subset$f14_nbs_clusters)))


# Identify markers vs rest
se_subset <- SetIdent(se_subset, value = "f14_nbs_clusters")
makers_f14_subcluster_nbs_all <- FindAllMarkers(se_subset)

#'export
fname <- paste0("mm_visium_nmf30_d21_f14hi_subclusters_nbs_cluster_markers_all.csv")
write.csv(x = makers_f14_subcluster_nbs_all, file = file.path(DIR_OBJ_NMF, fname), row.names = F)
makers_f14_subcluster_nbs_all <- read.csv(file = file.path(DIR_OBJ_NMF, fname))

# Check top genes
top_markers <- makers_f14_subcluster_nbs_all %>% 
  filter(cluster != "other") %>% 
  filter(avg_log2FC>0, p_val_adj<0.01, (pct.1-pct.2)>0) %>% 
  group_by(cluster) %>% 
  slice_max(n = 20, order_by = avg_log2FC)

# top_markers %>% print(n=40)

# Plot marker dotplot 
genes_plot <- top_markers$gene
p_dp <- DotPlot(se_subset,
                features = unique(genes_plot), 
                group.by = "f14_nbs_clusters",
                scale = T) +
  scale_color_gradientn(colours = c(rev(col_scale_mako[1:5]), col_scale_acton[1:8])) +
  theme_dotplot +
  theme(axis.text.y = element_text(face = "plain"),
        axis.title = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, face="italic"),
        legend.position = "right");p_dp

# proportions plots
d <- se_subset@meta.data %>%
  filter(f14_nbs_clusters != "other") %>% 
  group_by(sample_name, animal, f14_nbs_clusters) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

p_prop <- ggplot(d, aes(x=sample_name, y=n, fill=f14_nbs_clusters)) +
  geom_bar(stat = 'identity', colour=NA, position = "stack") +
  scale_fill_manual(values = cluster_cols_nbs_all) +
  scale_y_log10() +
  labs(x="", y="n spots (log10 scale)", fill="") +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        panel.grid = element_blank(),
        legend.position = "right", 
        text = element_text(size=10),
        plot.title = element_text(hjust=0.5));p_prop


#' Plot spatial
se_subset$f14_nbs_clusters2 <- ifelse(as.character(se_subset$f14_nbs_clusters) != "other", 
                                      as.character(se_subset$f14_nbs_clusters),
                                      ifelse(se_subset$annotation %in% "Suspect Fibrosis/Fibroplasia", 
                                             "fibrosis",
                                             "other"
                                      ))
# se_subset$f14_nbs_clusters2[is.na(se_subset$f14_nbs_clusters2)] <- "other"
cluster_cols_nbs_all2 <- c(cluster_cols_nbs_all, setNames("grey", "fibrosis"))

p_spat <- ST.FeaturePlot(se_subset,
                         features = "f14_nbs_clusters2", 
                         ncol=6, cols = cluster_cols_nbs_all2, 
                         label.by = "sample_name", 
                         pt.size = 0.75, show.sb = F) + 
  guides(fill = guide_legend(override.aes = list(size=2))) &
  theme(legend.position = "right", aspect.ratio = 1);p_spat


p_up <- (p_nbs_umap|p_spat)+plot_layout(widths = c(1,6))
p_do <- (p_prop|p_dp)+plot_layout(widths = c(1,6)) 

png(file.path(DIR_FIG_NMF, paste0("mm_visium_nmf30_d21_f14hi_subclusters_nbs_clusters_plots.png")), width = 16*fig_res, height = 8*fig_res, res = fig_res)
p_up/p_do+plot_layout(heights = c(2,1)) 
dev.off()


###### Export metadata #####
mdat_export <- se_subset[[]]
mdat_export$barcode <- rownames(mdat_export)
dim(mdat_export);head(mdat_export)

write.csv(x = mdat_export, file = file.path(DIR_OBJ_NMF, "mm_visium_nmf30_d21_se_processed_metadata.csv"), row.names = T)


##### Export NMF_d7 and NMF_d21 table ####
n_factors <- 30
DIR_OBJ_NMF <- file.path(DIR_RES, "objects", paste0("NMF", n_factors))

# d7
factor_gene_loadings_d7 <- lapply(1:n_factors, function(factor_x){
  read_xlsx(path = file.path(gsub("NMF30", "NMF30_d7", DIR_OBJ_NMF),
                             "mm_visium_nmf30_d7_top100_gene_loadings.xlsx"), sheet = factor_x)
})
factor_gprofiler_d7 <- lapply(1:n_factors, function(factor_x){
  read_xlsx(path = file.path(gsub("NMF30", "NMF30_d7", DIR_OBJ_NMF),
                             "mm_visium_nmf30_d7_gProfiler_top25genes.xlsx"), sheet = factor_x)
})


# d21
factor_gene_loadings_d21 <- lapply(1:n_factors, function(factor_x){
  read_xlsx(path = file.path(gsub("NMF30", "NMF30_d21", DIR_OBJ_NMF),
                             "mm_visium_nmf30_d21_top100_gene_loadings.xlsx"), sheet = factor_x)
})
factor_gprofiler_d21 <- lapply(1:n_factors, function(factor_x){
  read_xlsx(path = file.path(gsub("NMF30", "NMF30_d21", DIR_OBJ_NMF),
                             "mm_visium_nmf30_d21_gProfiler_top25genes.xlsx"), sheet = factor_x)
})

# Prep data
prepFactorExportData <- function(factor_gene_loading_list,
                                 factor_gprofiler_list){
  # Read data
  factor_gene_loadings_all <- bind_rows(factor_gene_loading_list, .id = "factor")
  factor_gene_loadings_all$factor <- as.numeric(factor_gene_loadings_all$factor)
 
  factor_gprofiler_all <- bind_rows(factor_gprofiler_list, .id = "factor")
  factor_gprofiler_all$factor <- as.numeric(factor_gprofiler_all$factor)
  
  n_factors <- max(factor_gene_loadings_all$factor)
  
  # Select top 50 genes
  factor_gene_loadings_all <- factor_gene_loadings_all %>% 
    filter(rank <= 50)
  
  # Collapse genes per factor
  factor_gene_loadings_coll <- factor_gene_loadings_all %>% 
    group_by(factor) %>% 
    summarize(Top.50.genes = paste(gene, collapse = "; "))
  
  # Select top 10 significant pathways
  factor_gprofiler_sign <- factor_gprofiler_all %>% 
    filter(significant == TRUE) %>%
    group_by(factor) %>% 
    slice_max(n=10, order_by = desc(p_value))
  
  # Collapse pathway per factor
  factor_gprofiler_coll <- factor_gprofiler_sign %>% 
    group_by(factor) %>% 
    summarize(Top.pathway.enrichment = paste(term_name, collapse = "; "))
  
  # join output
  df_out <- merge(factor_gene_loadings_coll, factor_gprofiler_coll, by = "factor")
  return(df_out)
}

# Create summarized data frames
df_out_d7 <- prepFactorExportData(factor_gene_loading_list = factor_gene_loadings_d7,
                                  factor_gprofiler_list = factor_gprofiler_d7
                                  )

df_out_d21 <- prepFactorExportData(factor_gene_loading_list = factor_gene_loadings_d21,
                                  factor_gprofiler_list = factor_gprofiler_d21
)

# Concatenate
df_out_d7$NMF <- "NMF30_d7"
df_out_d21$NMF <- "NMF30_d21"

df_out_all <- bind_rows(df_out_d7, df_out_d21)
df_out_all <- df_out_all[, c("NMF", "factor", "Top.50.genes", "Top.pathway.enrichment")]

dim(df_out_all)

write_xlsx(
  df_out_all,
  path = file.path(DIR_OBJ_NMF, "..", paste0(file_save_prefix, "d7_d21_factor_description.xlsx")),
  col_names = TRUE,
  format_headers = TRUE,
  use_zip64 = FALSE
)



