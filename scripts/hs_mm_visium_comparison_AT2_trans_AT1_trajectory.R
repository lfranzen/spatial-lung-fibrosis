#' [hs_mm_visium_comparison_AT2_trans_AT1_trajectory.R]
#'
#' Analysis of potential trajectory of AT2 -> Kt8+ADI / AbBa -> AT1
#' in BLM mouse (d21) and human IPF
#'
#' L. Franzén [lovisa.franzen@scilifelab.se_hs]

#### Set up ####
##### Define params. ####
# SPECIES <- "mouse"
# SPECIES <- "human"
SPECIES <- "translational"
DIR_ROOT <- getwd()
DIR_DATA <- file.path(DIR_ROOT, "data")
DIR_RES <- file.path(DIR_ROOT, "results", "translational")
DIR_FIG <- file.path(DIR_RES, "figures")
DIR_FIG_OUT <- file.path(DIR_FIG, "AT2_trajectory")
DIR_OBJ_OUT <- file.path(DIR_RES, "objects", "AT2_trajectory")
dir.create(DIR_FIG_OUT); dir.create(DIR_OBJ_OUT)

##### Load libs ####
# BiocManager::install("tradeSeq")

library(tidyverse)
library(dplyr)
library(STutility)
library(patchwork)
library(scales)
library(readxl)
library(writexl)
library(slingshot)
# library(tradeSeq)


##### Other ####
source(file.path(DIR_ROOT, "scripts", "custom_functions.R"))
source(file.path(DIR_ROOT, "scripts", "custom_colors.R"))
theme_custom <- theme(axis.title.x = element_blank())
fig_res <- 500


#### Gene conversion data ####
gene_conv_df <- read.csv(file.path(DIR_DATA, "misc", "orthogene_conv_combined_hs_mm.csv"), row.names = 1)
rownames(gene_conv_df) <- gene_conv_df$symbol_hs_mm
gene_conv <- setNames(gene_conv_df$symbol_hs_mm, nm = gene_conv_df$symbol_hs)
gene_conv_mm <- setNames(gene_conv_df$symbol_hs_mm, nm = gene_conv_df$symbol_mm)
gene_conv_df$symbol_hs_mm2 <- gsub("_", "-", gene_conv_df$symbol_hs_mm)


#### Read gene annotation data ####
mm_cell_anno <- read.csv(file.path(DIR_ROOT, "data", "misc", "strunz_cell_type_groups.csv"), sep = ";")
hs_cell_anno <- read.csv(file.path(DIR_ROOT, "data", "misc", "habermann_cell_type_groups.csv"), sep = ";")

rownames(mm_cell_anno) <- ifelse(!is.na(mm_cell_anno$cell_name), mm_cell_anno$cell_name, "NA")
rownames(hs_cell_anno) <- hs_cell_anno$cell_name


#### Read se_hs obj ####
##### Human #####
se_hs <- readRDS(file.path(DIR_RES, "..", "human", "objects", "hs_visium_preproc_A_se_obj_nmf.rds"))

# Add cell2location data
c2l_all <- read.csv(file.path(DIR_DATA, "human", "sc_deconvolution_habermann", "compiled_all_samples_cell_abundances.csv"), row.names = 1)
hs_cell_anno$c2l_colnames <- colnames(c2l_all)
se_hs <- AddMetaData(se_hs, c2l_all)

# update c2l cell names
c2l_names <- colnames(c2l_all)
new_c2l_names <- gsub("c2l_", "", gsub("[.]$", "", gsub("..", "_", c2l_names, fixed = TRUE)))
new_c2l_names_filtered <- new_c2l_names


##### Mouse #####
se_mm <- readRDS(file = file.path(DIR_RES, "..", "mouse", "objects", "mm_visium_preproc_se_obj_subset_d21_nmf30_c2l.rds"))
animal_ids <- se_mm$animal %>% unique()

mm_c2l_colnames <- grep("c2l_", colnames(se_mm[[]]), value = T)
mm_cell_anno$c2l_colnames <- mm_c2l_colnames


#### Trajectory analysis ####

##### Mouse #####
#' d21 NMF30: Cell type trajectory analysis ##### 
#' Further analysis of potential trajectory of AT2->Kt8+ADI->AT1
#' in the d21 subset NMF30 analysis
#' AT2 - Trans. AT2 - F14 (ADI) - AT1
#' 
file_save_prefix_mm <- paste0("mm_visium_nmf30", "_", "d21", "_")

###### Prep input data ######
# Set cell2location cell type densities as a reduction assay in order to use it for UMAP
c2l_columns <- grep("c2l_", colnames(se_mm@meta.data), value = T)
props_assay <- CreateAssayObject(data = as(t(se_mm@meta.data[,c2l_columns]), "matrix"))
se_mm[["cell2location"]] <- props_assay

c2l_columns_key <- setNames(1:length(c2l_columns), nm = c2l_columns)
feature.loadings <- matrix(nrow = length(c2l_columns), ncol = length(c2l_columns))
rownames(x = feature.loadings) <- paste0("c2l_", c2l_columns_key)
colnames(x = feature.loadings) <- paste0("c2l_", c2l_columns_key)
cell.embeddings <- se_mm@meta.data[,c2l_columns] %>% as.matrix()
rownames(x = cell.embeddings) <- rownames(se_mm@meta.data)
colnames(x = cell.embeddings) <- paste0("c2l_", c2l_columns_key)

reduction.data <- CreateDimReducObject (
  embeddings = cell.embeddings,
  loadings = feature.loadings,
  assay = DefaultAssay(se_mm),
  key = "c2l_"
)
se_mm[["c2l"]] <- reduction.data


# Select cell types and spots for UMAP and clustering
cells_select <- c("c2l_AT2.cells", "c2l_Activated.AT2.cells", "c2l_Krt8.ADI", "c2l_AT1.cells")
quant <- 0.95
cell_density <- se_mm@meta.data[,cells_select] %>% as.data.frame()
cell_density <- cell_density %>% 
  mutate(c_sum = rowSums(across(all_of(cells_select))),
         AT2_high = ifelse(c2l_AT2.cells>quantile(cell_density$c2l_AT2.cells, probs = quant), TRUE, FALSE),
         ActAT2_high = ifelse(c2l_Activated.AT2.cells>quantile(cell_density$c2l_Activated.AT2.cells, probs = quant), TRUE, FALSE),
         ADI_high = ifelse(c2l_Krt8.ADI>quantile(cell_density$c2l_Krt8.ADI, probs = quant), TRUE, FALSE),
         AT1_high = ifelse(c2l_AT1.cells>quantile(cell_density$c2l_AT1.cells, probs = quant), TRUE, FALSE),
  )

cell_density <- cell_density %>% 
  mutate(top_half = ifelse(c_sum > mean(cell_density$c_sum), TRUE, FALSE),
         cell_hi_sum = rowSums(across(contains("_high"))))

se_mm <- AddMetaData(se_mm, cell_density[,-c(1:4)])

# Subset data
se_subset_bleo_at <- SubsetSTData(se_mm, expression = condition == "bleomycin" & cell_hi_sum > 0)
se_subset_at <- se_subset_bleo_at
dim(se_subset_at)
cor(se_subset_at@reductions$NMF@cell.embeddings[,"factor_14"], se_subset_at$c2l_Krt8.ADI)
cor(se_subset_at@reductions$NMF@cell.embeddings[,"factor_14"], se_subset_at$c2l_Krt8.ADI, method = "spearman")


###### Create UMAP ######
# c2l-based UMAP
nn <- 30
md <- 0.1
se_subset_at <- RunUMAP(se_subset_at,
                        reduction = "c2l",
                        dims = c2l_columns_key[cells_select],
                        seed.use = 1,
                        n.neighbors = nn,
                        min.dist = md,
                        reduction.name = "umap_c2l", 
                        verbose = T)

p1 <- FeaturePlot(se_subset_at, reduction = "umap_c2l", 
                  pt.size = 0.25, 
                  max.cutoff = 2,
                  features = cells_select, 
                  ncol = length(cells_select)) &
  scale_colour_gradientn(colors = viridis::rocket(n = 8, direction = -1)) &
  theme_void() &
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5),
        legend.position = "top", panel.grid = element_blank())
p1 <- p1 + plot_annotation(caption = paste0("UMAP: n neigh ", nn, "; min dist ", md), 
                           theme = theme(plot.caption = element_text(hjust=0.5, face = "italic")));p1

# Cluster (c2l-based)
se_subset_at <- FindNeighbors(se_subset_at, 
                              reduction = "c2l",
                              dims = c2l_columns_key[cells_select])
se_subset_at <- FindClusters(se_subset_at, resolution = 0.2)
se_subset_at$seurat_clusters_reordered <- factor(as.numeric(se_subset_at$seurat_clusters),
                                                 levels = unique(as.numeric(se_subset_at$seurat_clusters)))
levels(se_subset_at$seurat_clusters_reordered) <- c("C4", "C2", "C1", "C3")
se_subset_at$seurat_clusters_reordered <- factor(se_subset_at$seurat_clusters_reordered, levels = sort(c("C4", "C2", "C1", "C3")))
cluster_cols <- setNames(hcl.colors(n = 4, palette = "Sunset"), nm = paste0("C",1:4))

p2 <- DimPlot(se_subset_at, 
              reduction = "umap_c2l", 
              pt.size = 0.5, 
              group.by = "seurat_clusters_reordered", 
              cols = cluster_cols) &
  labs(title = "Clusters") &
  theme_bw() &
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5),
        legend.position = "right", panel.grid = element_blank());p2


FeaturePlot(se_subset_at, 
            reduction = "umap_c2l", 
            pt.size = 0.1,
            features = c("factor_30", "factor_14", "factor_18", "factor_5"), ncol = 4) &
  scale_colour_gradientn(colors = viridis::mako(n = 8, direction = -1)) &
  theme_void() &
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5),
        legend.position = "top", panel.grid = element_blank())

se_subset_at$seurat_clusters_reordered2 <- factor(se_subset_at$seurat_clusters_reordered,
                                                  levels = rev(levels(se_subset_at$seurat_clusters_reordered)))
VlnPlot(se_subset_at, pt.size = 0, 
        ncol=1,
        cols = cluster_cols,
        features = c("c2l_AT2.cells", "c2l_Activated.AT2.cells", "c2l_Krt8.ADI", "c2l_AT1.cells"), 
        group.by = "seurat_clusters_reordered2") &
  coord_flip() &
  theme(plot.title = element_text(size = 14, face = "plain"), axis.title = element_blank())


###### Run Slingshot ######
# Slingshot analysis
#' https://nbisweden.github.io/workshop-archive/workshop-scRNAseq/2020-01-27/labs/compiled/slingshot/slingshot.html
umap_emb_c2l <- se_subset_at@reductions$umap_c2l@cell.embeddings
umap_emb <- umap_emb_c2l

clusts <- as.character(se_subset_at$seurat_clusters_reordered)
counts <- as.matrix(se_subset_at@assays$RNA@counts[se_subset_at@assays$SCT@var.features, ])

set.seed(1)
lineages <- getLineages(data = umap_emb,
                        clusterLabels = clusts,
                        start.clus = "C1") # define where to start the trajectories
lineages

# Get curves
curves <- getCurves(lineages, approx_points = 300, 
                    thresh = 0.01, stretch = 0.8, 
                    allow.breaks = FALSE, shrink = 0.99)
curves

par(mfrow = c(1, 1))
plot(umap_emb, col = pal[clusts], asp = 1, cex = 0.5, pch = 16)
lines(curves, lwd = 3, col = "black")

# Plot ggplot
sds <- curves
cl <- slingClusterLabels(sds)[, paste0("C", 1:4)]
cl <- apply(cl, 1, which.max)
cl <- paste0("C", as.character(cl))
rd <- reducedDim(sds)
df <- data.frame(dim1 = rd[, 1], dim2 = rd[, 2], col = cl)
clust_col <- setNames(hcl.colors(n = 4, palette = "Sunset"), nm = paste0("C", 1:4))
p <- ggplot() +
  geom_point(data = df, mapping = aes(x = dim1, y = dim2, col = col), size=0.5) +
  scale_color_manual(values = clust_col) +
  labs(x="UMAP 1", y="UMAP 2", title = "Slingshot trajectory", col="Cluster") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1, 
        plot.title = element_text(hjust=0.5)) +
  guides(color = guide_legend(override.aes = list(size = 3)));p

# Adding the curves
for (i in seq_along(slingCurves(curves))) {
  curve_i <- slingCurves(curves)[[i]]
  curve_i <- curve_i$s[curve_i$ord, ]
  colnames(curve_i) <- c("dim1", "dim2")
  p <- p + geom_path(data = as.data.frame(curve_i), 
                     aes(x = dim1, y = dim2), linewidth=1)
}
p

# Plot grid
pdf(file = file.path(DIR_FIG_OUT, paste0(file_save_prefix_mm, "slingshot_AT_celldensties_c2l-UMAP.pdf")), 
    width = 12, height = 4, useDingbats = F)
(p1)|p
dev.off()

pdf(file = file.path(DIR_FIG_OUT, paste0(file_save_prefix_mm, "slingshot_AT_celldensties_c2l-UMAP_clusters.pdf")), 
    width = 5, height = 2.5, useDingbats = F)
(p2&theme(legend.position = "none"))|p
dev.off()

###### Pseudo-time #####
dim(counts)
filt_counts <- counts[rowSums(counts > 5) > ncol(counts)/100, ]
dim(filt_counts)


#' run from scratch
sce <- fitGAM(counts = as.matrix(filt_counts), sds = curves)

saveRDS(sce, file.path(DIR_OBJ_OUT, paste0(file_save_prefix_hs, "slingshot_AT_celldensties_c2l-UMAP_fitGAM_object.rds")))

#' read pre-made object
sce <- readRDS(file.path(DIR_RES, "..", "mouse", "objects", "NMF30_d21", paste0(file_save_prefix_mm, "slingshot_AT_celldensties_c2l-UMAP_fitGAM_object.rds")))

#' fetch pseudo time
pseudotime_association <- associationTest(sce)
pseudotime_association$fdr <- p.adjust(pseudotime_association$pvalue, method = "fdr")
pseudotime_association <- pseudotime_association[order(pseudotime_association$pvalue), ]
pseudotime_association$feature_id <- rownames(pseudotime_association)

se_subset_bleo_at <- AddMetaData(se_subset_bleo_at, as.data.frame(sce$slingshot)) 

df2 <- bind_cols(df, data.frame(se_subset_bleo_at$pseudotime.curve1))
colnames(df2)[4] <- "pseudotime"
df2$pseudotime <- df2$pseudotime
df2$pseudotime_scaled <- scales::rescale(df2$pseudotime)

p_pt1 <- ggplot() +
  geom_point(data = df2, mapping = aes(x = dim1, y = dim2, fill = pseudotime), 
             size=2, shape = 21, stroke = 0.1) +
  scale_fill_gradientn(colors = viridis::viridis(n=10)) +
  labs(x="UMAP 1", y="UMAP 2", title = "Slingshot trajectory", color="pseudotime") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1, 
        plot.title = element_text(hjust=0.5)) +
  guides(color = guide_legend(override.aes = list(size = 3)))

p_pt2 <- ggplot() +
  geom_point(data = df2, mapping = aes(x = dim1, y = dim2, fill = pseudotime_scaled), 
             size=2, shape = 21, stroke = 0.1) +
  scale_fill_gradientn(colors = viridis::viridis(n=10)) +
  labs(x="UMAP 1", y="UMAP 2", title = "Slingshot trajectory", color="pseudotime") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1, 
        plot.title = element_text(hjust=0.5)) +
  guides(color = guide_legend(override.aes = list(size = 3)))


p_pt <- p_pt1
for (i in seq_along(slingCurves(curves))) {
  curve_i <- slingCurves(curves)[[i]]
  curve_i <- curve_i$s[curve_i$ord, ]
  colnames(curve_i) <- c("dim1", "dim2")
  p_pt <- p_pt + geom_path(data = as.data.frame(curve_i), 
                           aes(x = dim1, y = dim2), linewidth=1)
}; 
p_pt1 <- p_pt

p_pt <- p_pt2
for (i in seq_along(slingCurves(curves))) {
  curve_i <- slingCurves(curves)[[i]]
  curve_i <- curve_i$s[curve_i$ord, ]
  colnames(curve_i) <- c("dim1", "dim2")
  p_pt <- p_pt + geom_path(data = as.data.frame(curve_i), 
                           aes(x = dim1, y = dim2), linewidth=1)
}; 
p_pt2 <- p_pt


pdf(file = file.path(DIR_FIG_OUT, paste0(file_save_prefix_mm, "slingshot_AT_celldensties_c2l-UMAP_pseudotime.pdf")), 
    width = 5*1.5, height = 2.5*1.5, useDingbats = F)
p_pt1|p_pt2
dev.off()


# save table
df_export <- bind_cols(df2, data.frame(se_subset_bleo_at@meta.data[, cells_select]))

write.csv(df_export, file.path(DIR_FIG_OUT, paste0(file_save_prefix_mm, "slingshot_AT_celldensties_c2l-UMAP_pseudotime_plot_data.csv")), row.names = T)


##### Human #####
#' IPF: Cell type trajectory analysis #### 
#'Analysis of potential trajectory of AT2->KRT5-/KR17+->AT1
#' in the IPF subset
#' AT2 - Trans. AT2 - KRT5-/KR17+ - AT1
#' 
#' 
file_save_prefix_hs <- "hs_visium_c2l_res_IPF_"

###### Prep input data ######
# Set cell2location cell type densities as a reduction assay in order to use it for UMAP
c2l_columns <- grep("c2l_", colnames(se_hs@meta.data), value = T)
props_assay <- CreateAssayObject(data = as(t(se_hs@meta.data[,c2l_columns]), "matrix"))
se_hs[["cell2location"]] <- props_assay

c2l_columns_key <- setNames(1:length(c2l_columns), nm = c2l_columns)
feature.loadings <- matrix(nrow = length(c2l_columns), ncol = length(c2l_columns))
rownames(x = feature.loadings) <- paste0("c2l_", c2l_columns_key)
colnames(x = feature.loadings) <- paste0("c2l_", c2l_columns_key)
cell.embeddings <- se_hs@meta.data[,c2l_columns] %>% as.matrix()
rownames(x = cell.embeddings) <- rownames(se_hs@meta.data)
colnames(x = cell.embeddings) <- paste0("c2l_", c2l_columns_key)

reduction.data <- CreateDimReducObject (
  embeddings = cell.embeddings,
  loadings = feature.loadings,
  assay = DefaultAssay(se_hs),
  key = "c2l_"
)
se_hs[["c2l"]] <- reduction.data


# Select cell types and spots for UMAP and clustering
cells_select <- c("c2l_AT2", "c2l_Transitional.AT2", "c2l_KRT5..KRT17.", "c2l_AT1")

quant <- 0.95
cell_density <- se_hs@meta.data[,cells_select] %>% as.data.frame()
cell_density <- cell_density %>% 
  mutate(c_sum = rowSums(across(all_of(cells_select))),
         AT2_high = ifelse(c2l_AT2>quantile(cell_density$c2l_AT2, probs = quant), TRUE, FALSE),
         TransAT2_high = ifelse(c2l_Transitional.AT2>quantile(cell_density$c2l_Transitional.AT2, probs = quant), TRUE, FALSE),
         AbBa_high = ifelse(c2l_KRT5..KRT17.>quantile(cell_density$c2l_KRT5..KRT17., probs = quant), TRUE, FALSE),
         AT1_high = ifelse(c2l_AT1>quantile(cell_density$c2l_AT1, probs = quant), TRUE, FALSE)
         # Basal_high = ifelse(c2l_Basal>quantile(cell_density$c2l_Basal, probs = quant), TRUE, FALSE)
  )

cell_density <- cell_density %>% 
  mutate(top_half = ifelse(c_sum > mean(cell_density$c_sum), TRUE, FALSE),
         cell_hi_sum = rowSums(across(contains("_high"))))

se_hs <- AddMetaData(se_hs, cell_density[,-c(1:5)])

# Subset data
se_ipf_at <- SubsetSTData(se_hs, expression = condition == "IPF" & cell_hi_sum > 0)
dim(se_ipf_at)
cor(se_ipf_at@reductions$NMF@cell.embeddings[,"factor_14"], se_ipf_at$c2l_KRT5..KRT17.)
cor(se_ipf_at@reductions$NMF@cell.embeddings[,"factor_14"], se_ipf_at$c2l_KRT5..KRT17., method = "spearman")

###### Create UMAP ######
# c2l-based UMAP
nn <- 30
md <- 0.1
se_ipf_at <- RunUMAP(se_ipf_at,
                     reduction = "c2l",
                     dims = c2l_columns_key[cells_select],
                     seed.use = 1,
                     n.neighbors = nn,
                     min.dist = md,
                     reduction.name = "umap_c2l", 
                     verbose = T)

p1 <- FeaturePlot(se_ipf_at, reduction = "umap_c2l", 
                  pt.size = 0.25, 
                  max.cutoff = 3,
                  features = cells_select, 
                  ncol = length(cells_select)) &
  scale_colour_gradientn(colors = viridis::rocket(n = 8, direction = -1)) &
  theme_void() &
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5),
        legend.position = "top", panel.grid = element_blank());p1
p1 <- p1 + plot_annotation(caption = paste0("UMAP: n neigh ", nn, "; min dist ", md), 
                           theme = theme(plot.caption = element_text(hjust=0.5, face = "italic")));p1

# Cluster (c2l-based)
se_ipf_at <- FindNeighbors(se_ipf_at, 
                           reduction = "c2l",
                           dims = c2l_columns_key[cells_select])
se_ipf_at <- FindClusters(se_ipf_at, resolution = 0.1)
se_ipf_at$seurat_clusters_reordered <- factor(as.numeric(se_ipf_at$seurat_clusters),
                                              levels = unique(as.numeric(se_ipf_at$seurat_clusters)))

DimPlot(se_ipf_at, reduction = "umap_c2l", group.by = "seurat_clusters_reordered")

levels(se_ipf_at$seurat_clusters_reordered) <- c("C1", "C2", "C4", "C3")
se_ipf_at$seurat_clusters_reordered <- factor(se_ipf_at$seurat_clusters_reordered, levels = sort(c("C1", "C2", "C4", "C3")))
cluster_cols <- setNames(hcl.colors(n = 4, palette = "Sunset"), nm = paste0("C",1:4))

p2 <- DimPlot(se_ipf_at, 
              reduction = "umap_c2l", 
              pt.size = 0.5, 
              group.by = "seurat_clusters_reordered", 
              cols = cluster_cols) &
  labs(title = "Clusters") &
  theme_bw() &
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5),
        legend.position = "right", panel.grid = element_blank());p2

FeaturePlot(se_ipf_at, 
            reduction = "umap_c2l", 
            max.cutoff = 3,
            pt.size = 0.1,
            features = c("factor_19", "factor_14", "factor_12"), ncol = 3) &
  scale_colour_gradientn(colors = viridis::mako(n = 8, direction = -1)) &
  theme_void() &
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5),
        legend.position = "top", panel.grid = element_blank())

se_ipf_at$seurat_clusters_reordered2 <- factor(se_ipf_at$seurat_clusters_reordered,
                                               levels = rev(levels(se_ipf_at$seurat_clusters_reordered)))
VlnPlot(se_ipf_at, pt.size = 0, 
        ncol=1,
        cols = cluster_cols,
        features = cells_select, 
        group.by = "seurat_clusters_reordered2") &
  coord_flip() &
  theme(plot.title = element_text(size = 14, face = "plain"), axis.title = element_blank())


###### Run Slingshot ######
umap_emb <- se_ipf_at@reductions$umap_c2l@cell.embeddings
clusts <- as.character(se_ipf_at$seurat_clusters_reordered)
counts <- as.matrix(se_ipf_at@assays$RNA@counts[se_ipf_at@assays$SCT@var.features, ])

set.seed(1)
lineages <- getLineages(data = umap_emb,
                        clusterLabels = clusts,
                        start.clus = "C1") # define where to start the trajectories
lineages

#' Get curves
curves <- getCurves(lineages, approx_points = 300, 
                    thresh = 0.01, stretch = 0.8, 
                    allow.breaks = FALSE, shrink = 0.99)
curves

pal <- cluster_cols
par(mfrow = c(1, 1))
plot(umap_emb, col = pal[clusts], asp = 1, cex = 0.5, pch = 16)
lines(curves, lwd = 3, col = "black")


#' Plot ggplot
sds <- curves
cl <- slingClusterLabels(sds)[, paste0("C", 1:4)]
# cl <- slingClusterLabels(sds)
cl <- apply(cl, 1, which.max)
cl <- paste0("C", as.character(cl))
rd <- reducedDim(sds)
df <- data.frame(dim1 = rd[, 1], dim2 = rd[, 2], col = cl)
clust_col <- setNames(hcl.colors(n = 4, palette = "Sunset"), nm = paste0("C", 1:4))
p <- ggplot() +
  geom_point(data = df, mapping = aes(x = dim1, y = dim2, col = as.factor(col)), size=0.5) +
  scale_color_manual(values = clust_col) +
  labs(x="UMAP 1", y="UMAP 2", title = "Slingshot trajectory", col="Cluster") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1, 
        plot.title = element_text(hjust=0.5)) +
  guides(color = guide_legend(override.aes = list(size = 3)));p

# Adding the curves
for (i in seq_along(slingCurves(curves))) {
  curve_i <- slingCurves(curves)[[i]]
  curve_i <- curve_i$s[curve_i$ord, ]
  colnames(curve_i) <- c("dim1", "dim2")
  p <- p + geom_path(data = as.data.frame(curve_i), 
                     aes(x = dim1, y = dim2), linewidth=1)
}
p

# Plot grid
pdf(file = file.path(DIR_FIG_OUT, paste0(file_save_prefix_hs, "slingshot_AT_celldensties_c2l-UMAP.pdf")), 
    width = 12, height = 4, useDingbats = F)
(p1)|p
dev.off()

pdf(file = file.path(DIR_FIG_OUT, paste0(file_save_prefix_hs, "slingshot_AT_celldensties_c2l-UMAP_clusters.pdf")), 
    width = 5, height = 2.5, useDingbats = F)
(p2&theme(legend.position = "none"))|p
dev.off()



###### Pseudo-time #####
dim(counts)
filt_counts <- counts[rowSums(counts > 5) > ncol(counts)/100, ]
dim(filt_counts)

#' run from scratch
sce <- fitGAM(counts = as.matrix(filt_counts), sds = curves)

saveRDS(sce, file.path(DIR_OBJ_OUT, paste0(file_save_prefix_hs, "slingshot_AT_celldensties_c2l-UMAP_fitGAM_object.rds")))

#' read pre-made object
sce <- readRDS(file.path(DIR_RES, "..", "human", "objects", "cell2location_habermann", paste0(file_save_prefix_hs, "slingshot_AT_celldensties_c2l-UMAP_fitGAM_object.rds")))

#' fetch pseudo time
pseudotime_association <- associationTest(sce)
pseudotime_association$fdr <- p.adjust(pseudotime_association$pvalue, method = "fdr")
pseudotime_association <- pseudotime_association[order(pseudotime_association$pvalue), ]
pseudotime_association$feature_id <- rownames(pseudotime_association)

se_ipf_at <- AddMetaData(se_ipf_at, as.data.frame(sce$slingshot)) 

df2 <- bind_cols(df, data.frame(se_ipf_at@meta.data[, c("pseudotime.curve1", "pseudotime.curve2")]))
colnames(df2)[4] <- "pseudotime"
colnames(df2)[5] <- "pseudotime2"
df2$pseudotime_combined <- pmax(df2$pseudotime, df2$pseudotime2)
df2$pseudotime <- df2$pseudotime
df2$pseudotime_scaled <- scales::rescale(df2$pseudotime)
df2$pseudotime_comb_scaled <- scales::rescale(df2$pseudotime_combined)

p_pt1 <- ggplot() +
  geom_point(data = df2, mapping = aes(x = dim1, y = dim2, fill = pseudotime), 
             size=2, shape = 21, stroke = 0.1) +
  scale_fill_gradientn(colors = viridis::viridis(n=10)) +
  labs(x="UMAP 1", y="UMAP 2", title = "Slingshot trajectory", color="pseudotime") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1, 
        plot.title = element_text(hjust=0.5)) +
  guides(color = guide_legend(override.aes = list(size = 3)))

p_pt2 <- ggplot() +
  geom_point(data = df2, mapping = aes(x = dim1, y = dim2, fill = pseudotime_scaled), 
             size=2, shape = 21, stroke = 0.1) +
  scale_fill_gradientn(colors = viridis::viridis(n=10)) +
  labs(x="UMAP 1", y="UMAP 2", title = "Slingshot trajectory", color="pseudotime") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1, 
        plot.title = element_text(hjust=0.5)) +
  guides(color = guide_legend(override.aes = list(size = 3)))

p_pt3 <- ggplot() +
  geom_point(data = df2, mapping = aes(x = dim1, y = dim2, fill = pseudotime_comb_scaled), 
             size=2, shape = 21, stroke = 0.1) +
  scale_fill_gradientn(colors = viridis::viridis(n=10)) +
  labs(x="UMAP 1", y="UMAP 2", title = "Slingshot trajectory", color="pseudotime") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1, 
        plot.title = element_text(hjust=0.5)) +
  guides(color = guide_legend(override.aes = list(size = 3)))


p_pt <- p_pt1
for (i in seq_along(slingCurves(curves))) {
  curve_i <- slingCurves(curves)[[i]]
  curve_i <- curve_i$s[curve_i$ord, ]
  colnames(curve_i) <- c("dim1", "dim2")
  p_pt <- p_pt + geom_path(data = as.data.frame(curve_i), 
                           aes(x = dim1, y = dim2), linewidth=1)
}; 
p_pt1 <- p_pt

p_pt <- p_pt2
for (i in seq_along(slingCurves(curves))) {
  curve_i <- slingCurves(curves)[[i]]
  curve_i <- curve_i$s[curve_i$ord, ]
  colnames(curve_i) <- c("dim1", "dim2")
  p_pt <- p_pt + geom_path(data = as.data.frame(curve_i), 
                           aes(x = dim1, y = dim2), linewidth=1)
}; 
p_pt2 <- p_pt


pdf(file = file.path(DIR_FIG_OUT, paste0(file_save_prefix_hs, "slingshot_AT_celldensties_c2l-UMAP_pseudotime.pdf")), 
    width = 5*1.5, height = 2.5*1.5, useDingbats = F)
p_pt1|p_pt2
dev.off()

p_pt <- p_pt3
for (i in seq_along(slingCurves(curves))) {
  curve_i <- slingCurves(curves)[[i]]
  curve_i <- curve_i$s[curve_i$ord, ]
  colnames(curve_i) <- c("dim1", "dim2")
  p_pt <- p_pt + geom_path(data = as.data.frame(curve_i), 
                           aes(x = dim1, y = dim2), linewidth=1)
}; 
p_pt3 <- p_pt

pdf(file = file.path(DIR_FIG_OUT, paste0(file_save_prefix_hs, "slingshot_AT_celldensties_c2l-UMAP_pseudotime_2.pdf")), 
    width = 7.5*1.5, height = 2.5*1.5, useDingbats = F)
p_pt1|p_pt2|p_pt3
dev.off()


# save table
df_export <- bind_cols(df2, data.frame(se_ipf_at@meta.data[, cells_select]))

write.csv(df_export, file.path(DIR_FIG_OUT, paste0(file_save_prefix_hs, "slingshot_AT_celldensties_c2l-UMAP_pseudotime_plot_data.csv")), row.names = T)
write.csv(df_export, file.path(DIR_FIG_OUT, paste0(file_save_prefix_hs, "slingshot_AT_celldensties_c2l-UMAP_pseudotime_plot_data_2023-12-04.csv")), row.names = T)



#### Spatial trajectory cell co-localization plot ####
#' Continue by looking at spatial localization
#' of act/trans AT2 and KRT5-/KRT17+ or Krt8+ADI
#' and AT1 cells

##### Mouse #####
se_subset_bleo_at$AT2act_ADI_product <- rescale(se_subset_bleo_at$c2l_Activated.AT2.cells) * rescale(se_subset_bleo_at$c2l_Krt8.ADI)
se_subset_bleo_at$ADI_AT1_product <- rescale(se_subset_bleo_at$c2l_Krt8.ADI) * rescale(se_subset_bleo_at$c2l_AT1.cells)
se_subset_bleo_at$AD2_AT1_product <- rescale(se_subset_bleo_at$c2l_AT2.cells) * rescale(se_subset_bleo_at$c2l_AT1.cells)

ST.FeaturePlot(se_subset_bleo_at, features = c("c2l_Activated.AT2.cells", "c2l_Krt8.ADI", "c2l_AT1.cells"), 
               ncol = 4, show.sb = F, blend = T, dark.theme = T)

ST.FeaturePlot(se_subset_bleo_at, features = "AT2act_ADI_product", cols = c("black", "magenta"),
               ncol = 4, show.sb = F, dark.theme = T, label.by = "sample_name") & theme(aspect.ratio = 1)

ST.FeaturePlot(se_subset_bleo_at, features = "AD2_AT1_product", cols = c("black", "yellow"),
               ncol = 4, show.sb = F, dark.theme = T, label.by = "sample_name") & theme(aspect.ratio = 1)



p1 <- ST.FeaturePlot(se_subset_bleo_at, features = "annotation", 
                     cols = cols_annotation,
                     ncol = 3, show.sb = F, indices = c(4,6),
                     dark.theme = T, label.by = "sample_name") & theme(aspect.ratio = 1);p1

p2 <- ST.FeaturePlot(se_subset_bleo_at, features = "AT2act_ADI_product", cols = c("grey10", "magenta"),
                     max.cutoff = 0.15,
                     ncol = 3, show.sb = F, dark.theme = T, label.by = "sample_name", indices = c(4,6)) & theme(aspect.ratio = 1)

p3 <- ST.FeaturePlot(se_subset_bleo_at, features = "ADI_AT1_product", cols = c("grey10", "cyan"),
                     max.cutoff = 0.15,
                     ncol = 3, show.sb = F, dark.theme = T, label.by = "sample_name", indices = c(4,6)) & theme(aspect.ratio = 1)

p4 <- ST.FeaturePlot(se_subset_bleo_at, features = c("AT2act_ADI_product", "ADI_AT1_product"), channels.use = c("red", "blue"), 
                     ncol = 3, show.sb = F, blend = T, dark.theme = T, indices = c(4,6), max.cutoff = 0.15) & theme(aspect.ratio = 1)

p5 <- ST.FeaturePlot(se_subset_bleo_at, features = "ADI_AT1_product", cols = c("grey10", "yellow"),
                     max.cutoff = 0.15,
                     ncol = 3, show.sb = F, dark.theme = T, label.by = "sample_name", indices = c(4,6)) & theme(aspect.ratio = 1)


#' Final plot used for Fig. 5e:
se_subset_bleo_at$plot_pt_alpha <- pmax(
  scales::rescale(se_subset_bleo_at$AT2act_ADI_product), 
  scales::rescale(se_subset_bleo_at$ADI_AT1_product))

# Decide an alpha limit
se_subset_bleo_at$plot_pt_alpha_lim <- ifelse(se_subset_bleo_at$plot_pt_alpha >= 0.4, 0.4, se_subset_bleo_at$plot_pt_alpha)
se_subset_bleo_at$plot_pt_alpha_lim <- scales::rescale(se_subset_bleo_at$plot_pt_alpha_lim)

p6 <- ST.FeaturePlot(se_subset_bleo_at, features = c("AT2act_ADI_product", "ADI_AT1_product"), 
                     channels.use = c("red", "blue"), blend = T, 
                     pt.alpha = subset(se_subset_bleo_at@meta.data, sample_name %in% c("d21_bleo_4a", "d21_bleo_6a"))$plot_pt_alpha_lim,
                     ncol = 3, show.sb = F,
                     dark.theme = F, 
                     indices = c(4,6), max.cutoff = 0.15) & theme(aspect.ratio = 1);p6


pdf(file = file.path(DIR_FIG_OUT, paste0(file_save_prefix_mm, "AT_celldensity_products_selected_spatial.pdf")), 
    width = 5*1.5, height = 7*1.5, useDingbats = F)
p1/p2/p3/p4
dev.off()

pdf(file = file.path(DIR_FIG_OUT, paste0(file_save_prefix_mm, "AT_celldensity_products_selected_spatial_revised.pdf")), 
    width = 5*1.5, height = 7*1.5, useDingbats = F)
p1/p2/p3/p6
dev.off()

pdf(file = file.path(DIR_FIG_OUT, paste0(file_save_prefix_mm, "AT_celldensity_products_selected_spatial_revised_p6.pdf")), 
    width = 5*1.5, height = 2*1.5, useDingbats = F)
p6
dev.off()


##### Human #####
se_ipf_at$AT2act_AbBa_product <- rescale(se_ipf_at$c2l_Transitional.AT2) * rescale(se_ipf_at$c2l_KRT5..KRT17.)
se_ipf_at$AT2_AT1_product <- rescale(se_ipf_at$c2l_AT2) * rescale(se_ipf_at$c2l_AT1)
se_ipf_at$AT2act_AT1_product <- rescale(se_ipf_at$c2l_Transitional.AT2) * rescale(se_ipf_at$c2l_AT1)

ST.FeaturePlot(se_ipf_at, features = c("c2l_Transitional.AT2", "AT2act_AbBa_product", "c2l_KRT5..KRT17."), 
               ncol = 4, show.sb = F, blend = T, dark.theme = T)

ST.FeaturePlot(se_ipf_at, features = c("AT2act_AbBa_product", "AT2_AT1_product"), channels.use = c("red", "blue"),
               ncol = 4, show.sb = F, blend = T, dark.theme = T)

ST.FeaturePlot(se_ipf_at, features = "AT2act_AbBa_product", cols = c("black", "magenta"), 
               max.cutoff = 0.1,
               ncol = 4, show.sb = F, dark.theme = T, label.by = "sample_name")
ST.FeaturePlot(se_ipf_at, features = "AT2_AT1_product", cols = c("black", "cyan"), 
               max.cutoff = 0.1,
               ncol = 4, show.sb = F, dark.theme = T, label.by = "sample_name")


p1 <- ST.FeaturePlot(se_ipf_at, features = "annotation", cols = cols_annotation,
               ncol = 3, show.sb = F, dark.theme = T, label.by = "sample_name", indices = c(16:18)) & theme(aspect.ratio = 1)

p2 <- ST.FeaturePlot(se_ipf_at, features = "AT2act_AbBa_product", cols = c("grey10", "magenta"),
                     max.cutoff = 0.15,
               ncol = 3, show.sb = F, dark.theme = T, label.by = "sample_name", indices = c(16:18)) & theme(aspect.ratio = 1)

p3 <- ST.FeaturePlot(se_ipf_at, features = "AT2_AT1_product", cols = c("grey10", "cyan"),
                     max.cutoff = 0.15,
               ncol = 3, show.sb = F, dark.theme = T, label.by = "sample_name", indices = c(16:18)) & theme(aspect.ratio = 1)

p3.2 <- ST.FeaturePlot(se_ipf_at, features = "AT2act_AT1_product", cols = c("grey10", "blue"),
                     max.cutoff = 0.15,
                     ncol = 3, show.sb = F, dark.theme = T, label.by = "sample_name", indices = c(16:18)) & theme(aspect.ratio = 1)

p4 <- ST.FeaturePlot(se_ipf_at, features = c("AT2act_AbBa_product", "AT2_AT1_product"), channels.use = c("red", "blue"), 
                     ncol = 3, show.sb = F, blend = T, dark.theme = T, indices = c(16:18), max.cutoff = 0.15) & theme(aspect.ratio = 1)

p5 <- ST.FeaturePlot(se_ipf_at, features = c("AT2act_AbBa_product", "AT2act_AT1_product"), channels.use = c("red", "blue"), 
                     ncol = 3, show.sb = F, blend = T, dark.theme = T, indices = c(16:18), max.cutoff = 0.15) & theme(aspect.ratio = 1)


#' Final plot used for Fig. 5e:
# Add alpha to points
se_ipf_at$plot_pt_alpha <- pmax(
  scales::rescale(se_ipf_at$AT2act_AbBa_product), 
  scales::rescale(se_ipf_at$AT2_AT1_product))

se_ipf_at$plot_pt_alpha_lim <- ifelse(se_ipf_at$plot_pt_alpha >= 0.4, 0.4, se_ipf_at$plot_pt_alpha)
se_ipf_at$plot_pt_alpha_lim <- scales::rescale(se_ipf_at$plot_pt_alpha_lim)

p6 <- ST.FeaturePlot(se_ipf_at, features = c("AT2act_AbBa_product", "AT2_AT1_product"), 
                     channels.use = c("red", "blue"), blend = T, 
                     pt.alpha = subset(se_ipf_at@meta.data, sample_name %in% c("IPF_3.TD031.G1.2", "IPF_3.TD032A.G2.2", "IPF_3.TD032B.G3.2"))$plot_pt_alpha_lim,
                     ncol = 3, show.sb = F, 
                     dark.theme = F,
                     indices = c(16:18), max.cutoff = 0.15) & theme(aspect.ratio = 1)
p6

#' Export
pdf(file = file.path(DIR_FIG_OUT, paste0(file_save_prefix_hs, "AT_celldensity_products_selected_spatial.pdf")), 
    width = 7*1.5, height = 12*1.5, useDingbats = F)
p1/p2/p3/p3.2/p4/p5
dev.off()


pdf(file = file.path(DIR_FIG_OUT, paste0(file_save_prefix_hs, "AT_celldensity_products_selected_spatial_revised.pdf")), 
    width = 7*1.5, height = 12*1.5, useDingbats = F)
p1/p2/p3/p3.2/p4/p6
dev.off()

pdf(file = file.path(DIR_FIG_OUT, paste0(file_save_prefix_hs, "AT_celldensity_products_selected_spatial_revised_p6.pdf")), 
    width = 7*1.5, height = 3*1.5, useDingbats = F)
p6
dev.off()


