#' [mm_visium_radial_distance.R]
#'
#' Analyses related to the radial distance from selected regions
#'
#'
#' Mar 2023, L. Franzén [lovisa.franzen@scilifelab.se]

#### Set up ####
##### Define params. ####
set.seed(1)
SPECIES <- "mouse"
DIR_ROOT <- getwd()
DIR_DATA <- file.path(DIR_ROOT, "data", SPECIES, "visium")
DIR_RES <- file.path(DIR_ROOT, "results", SPECIES)
DIR_FIG <- file.path(DIR_RES, "figures")
DIR_FIG_OUT <- file.path(DIR_FIG, "radial_distance")
DIR_OBJ_OUT <- file.path(DIR_RES, "objects", "radial_distance")
dir.create(DIR_FIG_OUT); dir.create(DIR_OBJ_OUT)

##### Load libs ####
library(tidyverse)
library(dplyr)
library(STutility)
library(patchwork)

##### Other ####
source(file.path(DIR_ROOT, "scripts", "custom_functions.R"))
source(file.path(DIR_ROOT, "scripts", "custom_colors.R"))
theme_custom <- theme(axis.title.x = element_blank())
fig_res <- 300

##### Read objects ####
fname <- paste0("mm_visium_preproc_se_obj_nmf30_c2l.rds")
se <- readRDS(file = file.path(DIR_RES, "objects", fname))

#'Strunz cell2location data
# c2l_all <- read.csv(file.path(DIR_ROOT, "data", SPECIES, "sc_deconvolution_strunz", "compiled_all_samples_cell_abundances.csv"), 
#                     row.names = 1)
# 
# cell_type_col_names <- colnames(c2l_all)
# cell_type_names <- gsub("c2l_", "", gsub("[.]$", "", gsub("..", "_", cell_type_col_names, fixed = TRUE)))
# 
# se <- AddMetaData(se, c2l_all)


#### F14hiC0 r-dist #### 

##### Add data ####
subset_id <- "d21"
n_factors <- 30
DIR_OBJ_NMF <- file.path(DIR_RES, "objects", paste0("NMF", n_factors, "_", subset_id))
mdat <- read.csv(file = file.path(DIR_OBJ_NMF, "mm_visium_nmf30_d21_se_processed_metadata.csv"), row.names = 1)

rdist_mdata_f14hiC0 <- read.csv(file.path(DIR_OBJ_OUT, "mm_visium_NMF30_d21-F14hi-C0_distance_meta_data.csv"), row.names = 1)
rownames(rdist_mdata_f14hiC0) <- rdist_mdata_f14hiC0$bcs_stutility
rdist_mdata_f14hiC0 <- rdist_mdata_f14hiC0 %>% select(r_dist_F14hi_C0)

# se <- AddMetaData(se, metadata = mdat %>% select(f14_subclusters, f14_subclusters2, f14_nbs_clusters, f14_nbs_clusters2))
se_d21 <- SubsetSTData(se, day %in% "d21")
se_d21 <- AddMetaData(se_d21, metadata = mdat %>% select(f14_subclusters, f14_subclusters2, d_cF14hi_C0, f14_nbs_clusters, f14_nbs_clusters2))
se_d21 <- AddMetaData(se_d21, metadata = rdist_mdata_f14hiC0, col.name = "r_dist_F14hi_C0")


# Add rdistance meta data to se object and subset
samples_include <- subset(mdat, f14_subclusters2 == "F14hi_C0") %>% pull(sample_name) %>% unique()
se.rdist <- SubsetSTData(se_d21, sample_name %in% samples_include)


##### d21: Compute correlations #####
dist_cutoff <- 250 #500

###### Genes ######
# se2 <- SubsetSTData(se.rdist, d_cF14hi_C0 < 3)
se2 <- SubsetSTData(se.rdist, r_dist_F14hi_C0 < dist_cutoff)
se2 <- FindVariableFeatures(se2, nfeatures = 1000)
sel_genes <- se2@assays$SCT@var.features
rm(se2)

gene_data <- se.rdist[[]] %>% 
  bind_cols(FetchData(se.rdist, vars = sel_genes)) %>% 
  filter(r_dist_F14hi_C0 < dist_cutoff)
  # filter(d_cF14hi_C0 < 3)

gene_cor <- setNames(lapply(sel_genes, function(g){
  # cor_res <- cor.test(x = gene_data[["d_cF14hi_C0"]], y = gene_data[[g]])
  cor_res <- cor.test(x = gene_data[["r_dist_F14hi_C0"]], y = gene_data[[g]])
  cor_res <- cbind(cor_res$estimate, cor_res$p.value)
  colnames(cor_res) <- c("cor", "pval")
  rownames(cor_res) <- g
  return(cor_res)
}), nm = sel_genes)
gene_cor_df <- do.call("rbind", gene_cor) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  as_tibble() %>%
  filter(!is.na(cor)) %>% 
  arrange(-cor) %>% 
  mutate(
    pval_FDR = p.adjust(gene_cor_df$pval, method = "BH"),
    sign_005 = ifelse(gene_cor_df$pval < 0.05, TRUE, FALSE),
    sign_001 = ifelse(gene_cor_df$pval < 0.01, TRUE, FALSE),
    sign_padj_005 = ifelse(gene_cor_df$pval_FDR < 0.05, TRUE, FALSE),
    sign_padj_001 = ifelse(gene_cor_df$pval_FDR < 0.01, TRUE, FALSE)
  )

nrow(gene_cor_df)
sum(gene_cor_df$sign_padj_005)
sum(gene_cor_df$sign_padj_001)

write.csv(gene_cor_df, 
          file.path(DIR_OBJ_OUT, paste0("mm_visium_NMF30_d21-F14hi-C0_distance_gene_cor_", dist_cutoff, "um.csv")), 
          row.names = F)

###### Cell types ######
cell_types_remove <- grep("_NA|_Low.quality.cells", colnames(se.rdist[[]]), value = T)
cell_colnames <- grep("c2l_", colnames(se.rdist[[]]), value = T)
cell_colnames <- cell_colnames[-match(cell_types_remove,cell_colnames)]

# cell_data <- se.rdist[[]][, c("d_cF14hi_C0", cell_colnames)] %>% 
#   filter(d_cF14hi_C0 < 3)
cell_data <- se.rdist[[]][, c("r_dist_F14hi_C0", cell_colnames)] %>% 
  filter(r_dist_F14hi_C0 < dist_cutoff)

cell_cor <- setNames(lapply(cell_colnames, function(g){
  # cor_res <- cor.test(x = cell_data[["d_cF14hi_C0"]], y = cell_data[[g]])
  cor_res <- cor.test(x = cell_data[["r_dist_F14hi_C0"]], y = cell_data[[g]])
  cor_res <- cbind(cor_res$estimate, cor_res$p.value)
  colnames(cor_res) <- c("cor", "pval")
  rownames(cor_res) <- g
  return(cor_res)
}), nm = cell_colnames)

cell_cor_df <- do.call("rbind", cell_cor) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "cell_type") %>% 
  as_tibble() %>%
  arrange(-cor) %>% 
  mutate(
    pval_FDR = p.adjust(cell_cor_df$pval, method = "BH"),
    sign_005 = ifelse(cell_cor_df$pval < 0.05, TRUE, FALSE),
    sign_001 = ifelse(cell_cor_df$pval < 0.01, TRUE, FALSE)
  )
cell_cor_df$sign_padj_005 <- ifelse(cell_cor_df$pval_FDR < 0.05, TRUE, FALSE)
cell_cor_df$sign_padj_001 <- ifelse(cell_cor_df$pval_FDR < 0.01, TRUE, FALSE)
cell_cor_df$cell_type <- gsub("c2l_", "", cell_cor_df$cell_type)

write.csv(cell_cor_df,
          file.path(DIR_OBJ_OUT, paste0("mm_visium_NMF30_d21-F14hi-C0_distance_cell_cor_", dist_cutoff, "um.csv")), 
          row.names = F)

cell_cor_df <- read.csv(
  file.path(DIR_OBJ_OUT,
            paste0("mm_visium_NMF30_d21-F14hi-C0_distance_cell_cor_", dist_cutoff, "um.csv"))) %>% 
  as_tibble()


##### Plot #####
###### Spatial radial distance ###### 
p_spat <- ST.FeaturePlot(se.rdist, 
                         # features = "d_cF14hi_C0", 
                         features = "r_dist_F14hi_C0", 
                         cols = c("grey10", col_scale_earth), 
                         pt.size = 0.8, 
                         ncol = 3, 
                         max.cutoff = 2e3, 
                         show.sb = F, 
                         label.by = "sample_name") + 
  labs(title = "Radial distance (µm) from F14hi C0") &
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5), 
        legend.position = "bottom", 
        legend.text = element_text(angle=45, vjust = 1, hjust=1, color="black", size=10))

png(file = file.path(DIR_FIG_OUT, "mm_visium_F14hiC0_distance_spatial_plot.png"), width = 5*fig_res, height = 10*fig_res, res = fig_res);p_spat;dev.off()


###### Fig 5f: Cell type cor - one panel per cell type ###### 
# Selected cells
dist_cutoff_plot <- 1e3
cell_selected <- c("Krt8.ADI", 
                   "AT2.cells",
                   "Activated.AT2.cells", 
                   "AT1.cells")
# cell_selected <- cell_cor_df[cell_cor_df$cor>0,]$cell_type %>% unique()

# Plot selected cell types - all samples
cell_cor_df_p2 <- cell_cor_df %>% filter(cell_type %in% cell_selected)

plot_df <- se.rdist[[]] %>% 
  # filter(d_cF14hi_C0 < 5) %>% 
  filter(r_dist_F14hi_C0 < dist_cutoff_plot) %>% 
  # filter(r_dist_krt8adi_hi < dist_cutoff_plot) %>% 
  pivot_longer(all_of(paste0("c2l_", cell_selected)), 
               names_to = "cell_type", values_to = "density") %>% 
  mutate_at(.vars="cell_type", ~ factor(., levels = paste0("c2l_", cell_selected)))
plot_df$cell_type <- gsub("c2l_", "", plot_df$cell_type)
plot_df$cell_type <- factor(plot_df$cell_type, levels = cell_selected)

p_cell <- ggplot(plot_df, aes(x=r_dist_F14hi_C0, y=density, color = cell_type)) +
  # geom_smooth(method = "loess",
  #             linewidth=0.8, fill="grey70") +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
              linewidth=0.8, fill="grey70") +
  geom_vline(aes(xintercept = 0, color = "border"), 
             linetype = "dashed", linewidth=0.5, color="black") +
  facet_wrap(~cell_type, scales = "free_y", ncol = 1, strip.position = "right") +
  scale_color_manual(values = c("#830051", col_scale_spec)) +
  labs(x=bquote("Distance (µm) from" ~F14^hi~C0), y="Cell type density") +
  theme_classic() +
  theme(strip.text.y.right = element_text(angle=0, hjust=0, color = "black"), 
        strip.background = element_rect(color = NA), 
        legend.position = "none",
        axis.text.x = element_text(angle=45, vjust = 1, hjust=1),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color="black"),
        panel.grid.major.x = element_line(linewidth=0.5, color="grey90")
  );p_cell

pdf(file.path(DIR_FIG_OUT, "mm_visium_F14hiC0_distance_cell_selected.pdf"), width=3, height=3);p_cell;dev.off()  # !Main Figure
pdf(file.path(DIR_FIG_OUT, "mm_visium_F14hiC0_distance_cell_selected_GAM.pdf"), width=3, height=3);p_cell;dev.off()  # !Main Figure


# All cells
dist_cutoff_plot <- 1e3
cell_selected <- cell_cor_df$cell_type %>% sort()

plot_df2 <- se.rdist[[]] %>% 
  filter(r_dist_F14hi_C0 < dist_cutoff_plot) %>% 
  pivot_longer(all_of(paste0("c2l_", cell_selected)), 
               names_to = "cell_type", values_to = "density") %>% 
  mutate_at(.vars="cell_type", ~ factor(., levels = paste0("c2l_", cell_selected)))
plot_df2$cell_type <- gsub("c2l_", "", plot_df2$cell_type)
plot_df2$cell_type <- factor(plot_df2$cell_type, levels = cell_selected)

p_cell_all <- ggplot(plot_df2, aes(x=r_dist_F14hi_C0, y=density, color = cell_type)) +
  geom_smooth(method = "loess",
              linewidth=0.8, fill="grey70", color = "#830051") +
  geom_vline(aes(xintercept = 0, color = "border"), 
             linetype = "dashed", linewidth=0.5, color="black") +
  facet_wrap(~cell_type, scales = "free_y", ncol = 6) +
  labs(x=bquote("Distance (µm) from" ~F14^hi~C0), y="Cell type density") +
  theme_classic() +
  theme(
    strip.background = element_rect(color = NA), 
    legend.position = "none",
    axis.text.x = element_text(angle=45, vjust = 1, hjust=1),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color="black"),
    panel.grid.major.x = element_line(linewidth=0.5, color="grey90")
  );p_cell_all

pdf(file.path(DIR_FIG_OUT, "mm_visium_F14hiC0_distance_cell_all.pdf"), width=10, height=10);p_cell_all;dev.off()



#### TLS r-dist: all samples #### 
##### Add data ####
rdist_mdata_tls <- read.csv(file.path(DIR_OBJ_OUT, "ms_visium_se_blm2_TLS_distance_meta_data.csv"), row.names = 1)
rownames(rdist_mdata_tls) <- rdist_mdata_tls$bcs_stutility
rdist_mdata_tls$r_dist_TLS <- rdist_mdata_tls$r_dist_Inflammation
samples_include <- rdist_mdata_tls[!is.na(rdist_mdata_tls$r_dist_Inflammation), "sample_name"] %>% unique()

# Add rdistance meta data to se object and subset
# se.subset <- AddMetaData(se.subset, metadata = rdist_mdata_tls %>% select(r_dist_Inflammation))
# se.rdist <- SubsetSTData(se.subset, sample_name %in% samples_include)
# se.rdist$r_dist_TLS <- se.rdist$r_dist_Inflammation


# Add rdistance meta data to se object and subset
se <- AddMetaData(se, metadata = rdist_mdata_tls %>% select(r_dist_TLS))
se.rdist <- SubsetSTData(se, sample_name %in% samples_include)

##### Compute correlations #####
dist_cutoff <- 250

###### Genes ######
se2 <- SubsetSTData(se.rdist, r_dist_TLS < dist_cutoff) %>% 
  FindVariableFeatures(nfeatures = 1000)
sel_genes <- se2@assays$SCT@var.features
rm(se2)

gene_data <- se.rdist[[]] %>% 
  bind_cols(FetchData(se.rdist, vars = sel_genes)) %>% 
  filter(r_dist_TLS < dist_cutoff)

gene_cor <- setNames(lapply(sel_genes, function(g){
  cor_res <- cor.test(x = gene_data[["r_dist_TLS"]], y = gene_data[[g]])
  cor_res <- cbind(cor_res$estimate, cor_res$p.value)
  colnames(cor_res) <- c("cor", "pval")
  rownames(cor_res) <- g
  return(cor_res)
}), nm = sel_genes)

gene_cor_df <- do.call("rbind", gene_cor) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  as_tibble() %>%
  filter(!is.na(cor)) %>% 
  arrange(-cor)

gene_cor_df$pval_FDR <- p.adjust(gene_cor_df$pval, method = "BH")
gene_cor_df$sign_005 <- ifelse(gene_cor_df$pval < 0.05, TRUE, FALSE)
gene_cor_df$sign_001 <- ifelse(gene_cor_df$pval < 0.01, TRUE, FALSE)
gene_cor_df$sign_padj_005 <- ifelse(gene_cor_df$pval_FDR < 0.05, TRUE, FALSE)
gene_cor_df$sign_padj_001 <- ifelse(gene_cor_df$pval_FDR < 0.01, TRUE, FALSE)

nrow(gene_cor_df)
sum(gene_cor_df$sign_padj_005)
sum(gene_cor_df$sign_padj_001)

write.csv(gene_cor_df, 
          file.path(DIR_OBJ_OUT, paste0("mm_visium_TLS_rdistance_gene_cor_d21_", dist_cutoff, "um.csv")), 
          row.names = F)


###### Cell types ######
cell_colnames <- grep("c2l_", colnames(se.rdist[[]]), value = T)
cell_data <- se.rdist[[]][, c("r_dist_TLS", cell_colnames)] %>% 
  filter(r_dist_TLS < dist_cutoff)

cell_cor <- setNames(lapply(cell_colnames, function(g){
  cor_res <- cor.test(x = cell_data[["r_dist_TLS"]], y = cell_data[[g]])
  cor_res <- cbind(cor_res$estimate, cor_res$p.value)
  colnames(cor_res) <- c("cor", "pval")
  rownames(cor_res) <- g
  return(cor_res)
}), nm = cell_colnames)

cell_cor_df <- do.call("rbind", cell_cor) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "cell_type") %>% 
  as_tibble() %>%
  arrange(-cor)
cell_cor_df$cell_type <- gsub("c2l_", "", cell_cor_df$cell_type)

cell_cor_df$pval_FDR <- p.adjust(cell_cor_df$pval, method = "BH")
cell_cor_df$sign_005 <- ifelse(cell_cor_df$pval < 0.05, TRUE, FALSE)
cell_cor_df$sign_001 <- ifelse(cell_cor_df$pval < 0.01, TRUE, FALSE)
cell_cor_df$sign_padj_005 <- ifelse(cell_cor_df$pval_FDR < 0.05, TRUE, FALSE)
cell_cor_df$sign_padj_001 <- ifelse(cell_cor_df$pval_FDR < 0.01, TRUE, FALSE)

nrow(cell_cor_df)
sum(cell_cor_df$sign_padj_005)
sum(cell_cor_df$sign_padj_001)

write.csv(cell_cor_df,
          file.path(DIR_OBJ_OUT, paste0("mm_visium_TLS_rdistance_cell_cor_d21_", dist_cutoff, "um.csv")), 
          row.names = F)


##### Plot #####
###### Spatial radial distance ###### 
p_spat <- ST.FeaturePlot(se.rdist, features = "r_dist_TLS", 
                         cols = c("grey10", col_scale_temps), 
                         pt.size = 0.7, 
                         ncol = 2, 
                         max.cutoff = 2e3, 
                         show.sb = F, 
                         label.by = "sample_name") + 
  labs(title = "Radial distance TLS") &
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5), 
        legend.position = "bottom", 
        legend.text = element_text(angle=45, vjust = 1, hjust=1, color="black", size=10));p_spat

pdf(file = file.path(DIR_FIG_OUT, "mm_visium_TLS_rdistance_spatial_plot.pdf"), width = 4, height = 6);p_spat;dev.off()
png(file = file.path(DIR_FIG_OUT, "mm_visium_TLS_rdistance_spatial_plot.png"), width = 4*fig_res, height = 5.5*fig_res, res = fig_res);p_spat;dev.off()


###### Gene cor ###### 
dist_cutoff_plot <- 1500

# Plot top n pos. cor
n_genes <- 30
genes_pos <- gene_cor_df %>% arrange(desc(cor)) %>% head(n_genes) %>% select(gene) %>% unlist() %>% as.character()

p_dat_g_pos <- se.rdist[[]] %>% 
  bind_cols(FetchData(se.rdist, vars = genes_pos)) %>% 
  filter(r_dist_TLS < dist_cutoff_plot) %>% 
  pivot_longer(all_of(genes_pos), names_to = "gene", values_to = "expression") %>% 
  mutate_at(.vars="gene", ~ factor(., levels = rev(genes_pos)))
p_dat_g_pos$r_dist <- round(p_dat_g_pos$r_dist_TLS)

p_dat_g_pos2 <- p_dat_g_pos %>% 
  group_by(gene, r_dist) %>% 
  summarise(avg_expression = mean(expression)) %>% 
  ungroup() %>% 
  group_by(gene) %>% 
  mutate(avg_expr_scaled = scales::rescale(avg_expression, to = c(0, 1)))

p_pos_hm <- ggplot(p_dat_g_pos2, aes(x=r_dist, y=gene, fill=avg_expr_scaled, width=(max(r_dist)+10-r_dist)*0.1)) +
  geom_tile(alpha=0.9) +
  geom_vline(aes(xintercept = 0, color = "border"), 
             linetype = "dashed", linewidth=0.5, color="black") +
  labs(x="Distance (µm)", fill="Scaled avg.\nexpression", title="Positively correlated genes") +
  scale_fill_gradientn(colours = col_scale_lajolla) +
  theme_bw() +
  theme(text = element_text(color="black", size=10),
        plot.title = element_text(hjust=0.5, size=10),
        axis.title.y = element_blank(), 
        panel.grid = element_blank(),
        axis.text.x = element_text(angle=45, vjust = 1, hjust=1),
        legend.position = "bottom",
        legend.text = element_text(angle=45, vjust = 1, hjust=1, color="black", size=10))

# Plot top n neg. cor
n_genes <- 30
genes_neg <- gene_cor_df %>% arrange(cor) %>% head(n_genes) %>% select(gene) %>% unlist() %>% as.character()

p_dat_g_neg <- se.rdist[[]] %>% 
  bind_cols(FetchData(se.rdist, vars = genes_neg)) %>% 
  filter(r_dist_TLS < dist_cutoff_plot) %>% 
  pivot_longer(all_of(genes_neg), names_to = "gene", values_to = "expression") %>% 
  mutate_at(.vars="gene", ~ factor(., levels = rev(genes_neg)))
p_dat_g_neg$r_dist <- round(p_dat_g_neg$r_dist_TLS)

p_dat_g_neg2 <- p_dat_g_neg %>% 
  group_by(gene, r_dist) %>% 
  summarise(avg_expression = mean(expression)) %>% 
  ungroup() %>% 
  group_by(gene) %>% 
  mutate(avg_expr_scaled = scales::rescale(avg_expression, to = c(0, 1)))

p_neg_hm <- ggplot(p_dat_g_neg2, aes(x=r_dist, y=gene, fill=avg_expr_scaled, width=(max(r_dist)+10-r_dist)*0.1)) +
  geom_tile(alpha=0.9) +
  geom_vline(aes(xintercept = 0, color = "border"), 
             linetype = "dashed", linewidth=0.5, color="black") +
  labs(x="Distance (µm)", fill="Scaled avg.\nexpression", title="Negatively correlated genes") +
  scale_fill_gradientn(colours = col_scale_lajolla) +
  theme_bw() +
  theme(text = element_text(color="black", size=10),
        plot.title = element_text(hjust=0.5, size=10),
        axis.title.y = element_blank(), 
        panel.grid = element_blank(),
        axis.text.x = element_text(angle=45, vjust = 1, hjust=1),
        legend.position = "bottom",
        legend.text = element_text(angle=45, vjust = 1, hjust=1, color="black", size=10))

# Export
pdf(file.path(DIR_FIG_OUT, "mm_visium_TLS_rdistance_gene_top_pos_neg_heatmap.pdf"), width=6, height=6)
(p_pos_hm | p_neg_hm) + plot_layout(guides = 'collect') & theme(axis.text = element_text(color="black", size=10), axis.text.y = element_text(face = "italic"), legend.position = "bottom")
dev.off()



###### Cell type cor ###### 
dist_cutoff_plot <- 2e3

cell_cor_df_comp <- cell_cor_df %>% select(cell_type, cor)
colnames(cell_cor_df_comp)[2] <- "cor_TLS_d21"
cell_cor_df_comp <- pivot_longer(cell_cor_df_comp, cols = c("cor_TLS_d21"), names_to = "group", values_to = "cor")

cell_pval_df_comp <- cell_cor_df %>% select(cell_type, pval_FDR)
colnames(cell_pval_df_comp)[2] <- "cor_TLS_d21"
cell_pval_df_comp <- pivot_longer(cell_pval_df_comp, cols = c("cor_TLS_d21"), names_to = "group", values_to = "pval")
cell_pval_df_comp$sign_005 <- ifelse(cell_pval_df_comp$pval < 0.05, TRUE, FALSE)
cell_pval_df_comp$sign_001 <- ifelse(cell_pval_df_comp$pval < 0.01, TRUE, FALSE)

cell_df_comp <- bind_cols(cell_cor_df_comp, cell_pval_df_comp[,3:5])
cell_df_comp$cor_p <- ifelse(cell_df_comp$pval<0.05, cell_df_comp$cor, NA)

cor_lims <- cell_cor_df_comp$cor %>% abs() %>% max() %>% round(digits = 1) + 0.1

p1 <- ggplot(cell_df_comp, aes(x=group, y=reorder(cell_type, desc(cor)), fill=cor_p)) +
  geom_tile(width=1, height=0.8) +
  scale_fill_gradientn(colours = col_scale_div_custom2, limits=c(-cor_lims, cor_lims), na.value = "grey80") +
  theme_bw() +
  theme(text = element_text(color="black", size=10),
        plot.title = element_text(hjust=0.5, size=10),
        axis.title.y = element_blank(), 
        panel.grid = element_blank(),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(angle=0, vjust = 0.5, hjust=0.5, size=10),
        axis.title.x = element_blank(),
        legend.position = "right");p1

pdf(file.path(DIR_FIG_OUT, paste0("mm_visium_TLS_rdistance_cell_cor_", dist_cutoff, "um_heatmap.pdf")), width=3.6, height=6);p1;dev.off()


# Selected cells
cell_selected <- c("B.lymphocytes", 
                   "Themis..T.lymphocytes", 
                   "T.lymphocytes", 
                   "T.cell.subset",
                   "AT2.cells",
                   "Fibroblasts")

# Plot selected cell types
cell_cor_df_p2 <- cell_cor_df %>% filter(cell_type %in% cell_selected)

plot_df <- se.rdist[[]] %>% 
  filter(r_dist_TLS < dist_cutoff_plot) %>% 
  pivot_longer(all_of(paste0("c2l_", cell_selected)), 
               names_to = "cell_type", values_to = "density") %>% 
  mutate_at(.vars="cell_type", ~ factor(., levels = paste0("c2l_", cell_selected)))
plot_df$cell_type <- gsub("c2l_", "", plot_df$cell_type)
plot_df$cell_type <- factor(plot_df$cell_type, levels = cell_selected)

p_cell <- ggplot(plot_df, aes(x=r_dist_TLS, y=density, color = cell_type)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), 
              linewidth=0.8, fill="grey70") + # color="#830051"
  geom_vline(aes(xintercept = 0, color = "border"), 
             linetype = "dashed", linewidth=0.5, color="black") +
  facet_wrap(~cell_type, scales = "free_y", ncol = 1, strip.position = "right") +
  scale_color_manual(values = rev(col_scale_mako)) +
  labs(x="Distance (µm) from TLS", y="Cell type density") +
  theme_classic() +
  theme(strip.text.y.right = element_text(angle=0, hjust=0), 
        strip.background = element_rect(color = NA), 
        legend.position = "none",
        axis.text.x = element_text(angle=45, vjust = 1, hjust=1), 
        panel.grid.major.x = element_line(linewidth=0.5, color="grey90")
  );p_cell

pdf(file.path(DIR_FIG_OUT, "mm_visium_TLS_distance_cell_selected.pdf"), width=3, height=4);p_cell;dev.off()

# Plot per donor
plot_df$animal <- factor(plot_df$animal, levels = plot_df$animal %>% unique() %>% sort())
p_cell_animal <- ggplot(plot_df, aes(x=r_dist_TLS, y=density, color = cell_type)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), 
              linewidth=0.8, fill="grey70") + # color="#830051"
  geom_vline(aes(xintercept = 0, color = "border"), 
             linetype = "dashed", linewidth=0.5, color="black") +
  facet_grid(vars(cell_type), vars(animal), scales = "free_y") +
  # facet_wrap(~cell_type + animal, scales = "free_y", ncol = 5, strip.position = "top") +
  scale_color_manual(values = rev(col_scale_mako)) +
  labs(x="Distance (µm) from TLS", y="Cell type density") +
  theme_classic() +
  theme(
    # strip.text.y.right = element_text(angle=0, hjust=0), 
    strip.background = element_rect(color = NA), 
    strip.text.x = element_text(margin = margin(2, 0, 2, 0)),
    legend.position = "none",
    axis.text.x = element_text(angle=45, vjust = 1, hjust=1), 
    panel.grid.major.x = element_line(linewidth=0.5, color="grey90")
  );p_cell_animal

pdf(file.path(DIR_FIG_OUT, "mm_visium_TLS_distance_cell_selected_per_animal.pdf"), width=7, height=7);p_cell_animal;dev.off()


