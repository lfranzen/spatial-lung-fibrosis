#' [hs_visium_radial_distance.R]
#'
#' Analyses related to the radial distance from selected regions
#'
#'
#' Mar 2023, L. Franzén [lovisa.franzen@scilifelab.se]

#### Set up ####
##### Define params. ####
set.seed(1)
SPECIES <- "human"
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
library(writexl)
library(pheatmap)

##### Other ####
source(file.path(DIR_ROOT, "scripts", "custom_functions.R"))
source(file.path(DIR_ROOT, "scripts", "custom_colors.R"))
theme_custom <- theme(axis.title.x = element_blank())
fig_res <- 300

##### Read objects ####
fname <- paste0("hs_visium_preproc_A_se_obj_nmf.rds")
se.subset <- readRDS(file = file.path(DIR_RES, "objects", fname))

#'Habermann cell2location data
c2l_all <- read.csv(file.path(DIR_ROOT, "data", SPECIES, "sc_deconvolution_habermann", "compiled_all_samples_cell_abundances.csv"), 
                    row.names = 1)
cell_type_col_names <- colnames(c2l_all)
cell_type_names <- gsub("c2l_", "", gsub("[.]$", "", gsub("..", "_", cell_type_col_names, fixed = TRUE)))

se.subset <- AddMetaData(se.subset, c2l_all)


##### Custom Functions ####
RunGeneDistanceCor <- function(
    object, 
    dist_colname,
    genes_select, 
    dist_cutoff = NULL,
    dist_lower_cutoff = 0){
  
  if (is.null(dist_cutoff)) dist_cutoff <- 500
  
  gene_data <- object[[]] %>% 
    bind_cols(FetchData(object, vars = genes_select)) %>% 
    filter(.[[dist_colname]] < dist_cutoff & .[[dist_colname]] >= dist_lower_cutoff)
  
  message(paste0("Selecting ", length(genes_select), " genes"))
  message(paste0("Selecting ", nrow(gene_data), " spots"))
  
  message()
  
  gene_cor <- setNames(lapply(genes_select, function(g){
    cor_res <- cor.test(x = gene_data[[dist_colname]], y = gene_data[[g]])
    cor_res <- cbind(cor_res$estimate, cor_res$p.value)
    colnames(cor_res) <- c("cor", "pval")
    rownames(cor_res) <- g
    return(cor_res)
  }), nm = genes_select)
  
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
  
  return(gene_cor_df)
}

RunCellDistanceCor <- function(
    object,
    dist_colname,
    cell_colnames, 
    dist_cutoff = NULL,
    dist_lower_cutoff = 0){
  
  if (is.null(dist_cutoff)) dist_cutoff <- 500
  
  cell_data <- object[[]][, c(dist_colname, cell_colnames)] %>% 
    filter(.[[dist_colname]] < dist_cutoff & .[[dist_colname]] >= dist_lower_cutoff)
  
  cell_cor <- setNames(lapply(cell_colnames, function(g){
    cor_res <- cor.test(x = cell_data[[dist_colname]], y = cell_data[[g]])
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
  
  cell_cor_df$pval_FDR <- p.adjust(cell_cor_df$pval, method = "BH")
  cell_cor_df$sign_005 <- ifelse(cell_cor_df$pval < 0.05, TRUE, FALSE)
  cell_cor_df$sign_001 <- ifelse(cell_cor_df$pval < 0.01, TRUE, FALSE)
  cell_cor_df$sign_padj_005 <- ifelse(cell_cor_df$pval_FDR < 0.05, TRUE, FALSE)
  cell_cor_df$sign_padj_001 <- ifelse(cell_cor_df$pval_FDR < 0.01, TRUE, FALSE)
  
  return(cell_cor_df)
}

#### F14hiC0 r-dist: all samples #### 
##### Add data ####
rdist_mdata_f14hiC0 <- read.csv(file.path(DIR_OBJ_OUT, "hs_visium_A-NMF30-F14hi-C0_distance_meta_data.csv"), row.names = 1)
rdist_mdata_f14hiC0_nosingle <- read.csv(file.path(DIR_OBJ_OUT, "hs_visium_A-NMF30-F14hi-C0_distance_nosingle_meta_data.csv"), row.names = 1)
rdist_mdata_f14hiC0 <- bind_cols(rdist_mdata_f14hiC0, rdist_mdata_f14hiC0_nosingle %>% select(r_dist_F14_C0_nosingle))

rdist_mdata_f14hiC0$r_dist_F14_C0 <- rdist_mdata_f14hiC0$r_dist_F14_C0_nosingle

rownames(rdist_mdata_f14hiC0) <- rdist_mdata_f14hiC0$bcs_stutility
samples_include <- rdist_mdata_f14hiC0[!is.na(rdist_mdata_f14hiC0$r_dist_F14_C0), "sample_name"] %>% unique()


# Update output paths depending on singleton status
DIR_OBJ_OUT <- file.path(DIR_OBJ_OUT, "excl_singletons")
DIR_FIG_OUT <- file.path(DIR_FIG_OUT, "excl_singletons")


# Add rdistance meta data to se object and subset
se.subset <- AddMetaData(se.subset, metadata = rdist_mdata_f14hiC0 %>% select(r_dist_F14_C0))
se.rdist <- SubsetSTData(se.subset, sample_name %in% samples_include)


##### Compute correlations #####
dist_cutoff <- 500 

###### Genes ######
se2 <- SubsetSTData(se.rdist, r_dist_F14_C0 < dist_cutoff) %>% 
  FindVariableFeatures(nfeatures = 1000)
sel_genes <- se2@assays$SCT@var.features
rm(se2)

gene_data <- se.rdist[[]] %>% 
  bind_cols(FetchData(se.rdist, vars = sel_genes)) %>% 
  filter(r_dist_F14_C0 < dist_cutoff)

gene_cor <- setNames(lapply(sel_genes, function(g){
  cor_res <- cor.test(x = gene_data[["r_dist_F14_C0"]], y = gene_data[[g]])
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
          file.path(DIR_OBJ_OUT, paste0("hs_visium_A-NMF30-F14hi-C0_distance_gene_cor_", dist_cutoff, "um_180923.csv")), 
          row.names = F)

gene_cor_df <- read.csv(file.path(DIR_OBJ_OUT, paste0("hs_visium_A-NMF30-F14hi-C0_distance_gene_cor_", 500, "um_180923.csv")))

###### Cell types ######
cell_colnames <- grep("c2l_", colnames(se.rdist[[]]), value = T)
cell_data <- se.rdist[[]][, c("r_dist_F14_C0", cell_colnames)] %>% 
  filter(r_dist_F14_C0 < dist_cutoff)

cell_cor <- setNames(lapply(cell_colnames, function(g){
  cor_res <- cor.test(x = cell_data[["r_dist_F14_C0"]], y = cell_data[[g]])
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
          file.path(DIR_OBJ_OUT, paste0("hs_visium_A-NMF30-F14hi-C0_distance_cell_cor_", dist_cutoff, "um.csv")), 
          row.names = F)

cell_cor_df <- read.csv(
  file.path(DIR_OBJ_OUT, "excl_singletons",
            paste0("hs_visium_A-NMF30-F14hi-C0_distance_cell_cor_", dist_cutoff, "um.csv")))

###### Factors ######
sel_factors <- paste0("factor_", 1:30)
nmf_data <- se.rdist[[]] %>% 
  bind_cols(FetchData(se.rdist, vars = sel_factors)) %>% 
  filter(r_dist_F14_C0 < dist_cutoff)

nmf_cor <- setNames(lapply(sel_factors, function(g){
  cor_res <- cor.test(x = nmf_data[["r_dist_F14_C0"]], y = nmf_data[[g]])
  cor_res <- cbind(cor_res$estimate, cor_res$p.value)
  colnames(cor_res) <- c("cor", "pval")
  rownames(cor_res) <- g
  return(cor_res)
}), nm = sel_factors)

nmf_cor_df <- do.call("rbind", nmf_cor) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "factor_n") %>% 
  as_tibble() %>%
  filter(!is.na(cor)) %>% 
  arrange(-cor)

nmf_cor_df$pval_FDR <- p.adjust(nmf_cor_df$pval, method = "BH")
nmf_cor_df$sign_005 <- ifelse(nmf_cor_df$pval < 0.05, TRUE, FALSE)
nmf_cor_df$sign_001 <- ifelse(nmf_cor_df$pval < 0.01, TRUE, FALSE)
nmf_cor_df$sign_padj_005 <- ifelse(nmf_cor_df$pval_FDR < 0.05, TRUE, FALSE)
nmf_cor_df$sign_padj_001 <- ifelse(nmf_cor_df$pval_FDR < 0.01, TRUE, FALSE)

nrow(nmf_cor_df)
sum(nmf_cor_df$sign_padj_005)
sum(nmf_cor_df$sign_padj_001)

write.csv(nmf_cor_df, 
          file.path(DIR_OBJ_OUT, paste0("hs_visium_A-NMF30-F14hi-C0_distance_nmf_cor_", dist_cutoff, "um.csv")), 
          row.names = F)


##### Plot #####
###### Spatial radial distance ###### 
p_spat <- ST.FeaturePlot(se.rdist, features = "r_dist_F14_C0", 
               cols = c("grey10", col_scale_earth), 
               pt.size = 0.7, 
               ncol = 3, 
               max.cutoff = 2e3, 
               show.sb = F, 
               label.by = "sample_name") + 
  labs(title = "Radial dist F14hi-C0") &
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5), 
        legend.position = "bottom", 
        legend.text = element_text(angle=45, vjust = 1, hjust=1, color="black", size=10));p_spat

pdf(file = file.path(DIR_FIG_OUT, "hs_visium_A_nmf_1-30_f14high_C0_rdistance_spatial_plot.pdf"), width = 5, height = 10);p_spat;dev.off()
png(file = file.path(DIR_FIG_OUT, "hs_visium_A_nmf_1-30_f14high_C0_rdistance_spatial_plot.png"), width = 5*fig_res, height = 10*fig_res, res = fig_res);p_spat;dev.off()



###### Gene cor expr heatmap ###### 
dist_cutoff_plot <- 1500

# Plot top n pos. cor
n_genes <- 30
genes_pos <- gene_cor_df %>% arrange(desc(cor)) %>% head(n_genes) %>% select(gene) %>% unlist() %>% as.character()

p_dat_g_pos <- se.rdist[[]] %>% 
  bind_cols(FetchData(se.rdist, vars = genes_pos)) %>% 
  filter(r_dist_F14_C0 < dist_cutoff_plot) %>% 
  pivot_longer(all_of(genes_pos), names_to = "gene", values_to = "expression") %>% 
  mutate_at(.vars="gene", ~ factor(., levels = rev(genes_pos)))
p_dat_g_pos$r_dist <- round(p_dat_g_pos$r_dist_F14_C0)

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
  filter(r_dist_F14_C0 < dist_cutoff_plot) %>% 
  pivot_longer(all_of(genes_neg), names_to = "gene", values_to = "expression") %>% 
  mutate_at(.vars="gene", ~ factor(., levels = rev(genes_neg)))
p_dat_g_neg$r_dist <- round(p_dat_g_neg$r_dist_F14_C0)

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
pdf(file.path(DIR_FIG_OUT, "hs_visium_A-NMF30-F14hi-C0_distance_gene_top_pos_neg_heatmap.pdf"), width=6, height=6)
(p_pos_hm | p_neg_hm) + plot_layout(guides = 'collect') & theme(axis.text = element_text(color="black", size=10), axis.text.y = element_text(face = "italic"), legend.position = "bottom")
dev.off()


p_dat_g_concat <- bind_rows(p_dat_g_pos, p_dat_g_neg)

write.csv(p_dat_g_concat, 
          file.path(DIR_OBJ_OUT, paste0("hs_visium_A-NMF30-F14hi-C0_distance_gene_top_pos_neg_data.csv")), 
          row.names = F)


###### Fig 3d: Gene cor expr line ###### 
dist_cutoff_plot <- 500

# Plot top n pos. cor
n_genes <- 20
genes_pos <- gene_cor_df %>% arrange(desc(cor)) %>% head(n_genes) %>% filter(sign_001 == TRUE) %>% pull(gene)

p_dat_g_pos_line <- se.rdist[[]] %>% 
  bind_cols(FetchData(se.rdist, vars = c(genes_pos))) %>% 
  filter(r_dist_F14_C0 < dist_cutoff_plot) %>% 
  pivot_longer(all_of(c(genes_pos)), names_to = "gene", values_to = "expression") %>% 
  mutate_at(.vars="gene", ~ factor(., levels = c(genes_pos)))
p_dat_g_pos_line$r_dist <- round(p_dat_g_pos_line$r_dist_F14_C0)

p_dat_g_pos_cor <- gene_cor_df %>% arrange(desc(cor)) %>% filter(gene %in% genes_pos)
p_dat_g_pos_cor$gene <- factor(p_dat_g_pos_cor$gene, levels = genes_pos)

# Plot top n neg. cor
n_genes <- 20
genes_neg <- gene_cor_df %>% arrange(cor) %>% head(n_genes) %>% filter(sign_001 == TRUE) %>% pull(gene)

p_dat_g_neg_line <- se.rdist[[]] %>% 
  bind_cols(FetchData(se.rdist, vars = c(genes_neg))) %>% 
  filter(r_dist_F14_C0 < dist_cutoff_plot) %>% 
  pivot_longer(all_of(c(genes_neg)), names_to = "gene", values_to = "expression") %>% 
  mutate_at(.vars="gene", ~ factor(., levels = c(genes_neg)))
p_dat_g_neg_line$r_dist <- round(p_dat_g_neg_line$r_dist_F14_C0)

p_dat_g_neg_cor <- gene_cor_df %>% arrange((cor)) %>% filter(gene %in% genes_neg)
p_dat_g_neg_cor$gene <- factor(p_dat_g_neg_cor$gene, levels = genes_neg)


# plot
# pos:
p_g_pos_line <- ggplot(p_dat_g_pos_line, aes(r_dist_F14_C0, expression)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), linewidth=0.5, fill="grey80", color = "black") +
  geom_vline(aes(xintercept = 0, color = "border"), 
             linetype = "dashed", linewidth=0.5, color="black") +
  facet_wrap(~gene, scales = "free_y", ncol = 1, strip.position = "right") +
  theme_classic() +
  theme(strip.text.y.right = element_text(angle=0, hjust=0), 
        strip.background = element_rect(color = NA), 
        legend.position = "none",
        axis.text.x = element_text(angle=45, vjust = 1, hjust=1), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.title.y = element_blank(),
        panel.grid.major.x = element_line(linewidth=0.5, color="grey90")
  ); p_g_pos_line

cor_lims <- max(abs(p_dat_g_pos_cor$cor))
p_g_pos_cor <- ggplot(p_dat_g_pos_cor, aes(x="cor", y=reorder(gene, cor), fill=cor)) +
  geom_tile(height = 0.75) +
  scale_fill_gradientn(colours = col_scale_div_custom2, limits = c(-cor_lims, cor_lims)) +
  theme_classic() +
  theme(text = element_text(color="black", size=10),
        axis.title = element_blank(), 
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "left",
        legend.text = element_text(color="black", size=10)); p_g_pos_cor

p1 <- (p_g_pos_cor | p_g_pos_line) + plot_layout(widths = c(1,3))

# neg:
p_g_neg_line <- ggplot(p_dat_g_neg_line, aes(r_dist_F14_C0, expression)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), linewidth=0.5, fill="grey80", color = "black") +
  geom_vline(aes(xintercept = 0, color = "border"), 
             linetype = "dashed", linewidth=0.5, color="black") +
  facet_wrap(~gene, scales = "free_y", ncol = 1, strip.position = "right") +
  theme_classic() +
  theme(strip.text.y.right = element_text(angle=0, hjust=0), 
        strip.background = element_rect(color = NA), 
        legend.position = "none",
        axis.text.x = element_text(angle=45, vjust = 1, hjust=1), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.title.y = element_blank(),
        panel.grid.major.x = element_line(linewidth=0.5, color="grey90")
  ); p_g_neg_line

cor_lims <- max(abs(p_dat_g_neg_cor$cor))
p_g_neg_cor <- ggplot(p_dat_g_neg_cor, aes(x="cor", y=reorder(gene, desc(cor)), fill=cor)) +
  geom_tile(height = 0.75) +
  scale_fill_gradientn(colours = col_scale_div_custom2, limits = c(-cor_lims, cor_lims)) +
  theme_classic() +
  theme(text = element_text(color="black", size=10),
        axis.title = element_blank(), 
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "left",
        legend.text = element_text(color="black", size=10)); p_g_neg_cor

p2 <- (p_g_neg_cor | p_g_neg_line) + plot_layout(widths = c(1,3))


# Export
pdf(file.path(DIR_FIG_OUT, "hs_visium_A-NMF30-F14hi-C0_distance_gene_top_pos_neg_line.pdf"), width=3.5, height=6)  # !Main Figure
p1
p2
dev.off()

###### Selected genes ##### 

# Plot Keratins
genes_krts <- grep("^KRT", gene_cor_df$gene, value = T) %>% rev()

p_gene_ktr <- se.rdist[[]] %>% 
  bind_cols(FetchData(se.rdist, vars = genes_krts)) %>% 
  filter(r_dist_F14_C0 < dist_cutoff_plot) %>% 
  pivot_longer(all_of(genes_krts), names_to = "gene", values_to = "expression") %>% 
  mutate_at(.vars="gene", ~ factor(., levels = genes_krts)) %>% 
  ggplot(aes(r_dist_F14_C0, expression, color = gene)) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), linewidth=1, fill="grey80") +
    geom_vline(aes(xintercept = 0, color = "border"), 
               linetype = "dashed", linewidth=0.5, color="black") +
    facet_wrap(~gene, scales = "free_y", ncol = 1, strip.position = "right") +
    scale_color_manual(values = c(col_scale_lajolla[2:9] %>% rev())) +
    theme_classic() +
    theme(strip.text.y.right = element_text(angle=0, hjust=0), 
          strip.background = element_rect(color = NA), 
          legend.position = "none",
          axis.text.x = element_text(angle=45, vjust = 1, hjust=1), 
          panel.grid.major.x = element_line(linewidth=0.5, color="grey90")
    );p_gene_ktr
pdf(file.path(DIR_FIG_OUT, "hs_visium_A-NMF30-F14hi-C0_distance_gene_keratins.pdf"), width=2.5, height=4);p_gene_ktr;dev.off()


# Plot APOE and selected receptors
genes_apoe <- c("APOE", "GPC1", "LRP1", "LRP5", "SDC1", "SDC2", "SDC3", "SDC4")
cols_apoe <- c("#830051", 
               col_scale_mako[8], 
               col_scale_mako[5], col_scale_mako[5],
               col_scale_mako[3], col_scale_mako[3],col_scale_mako[3], col_scale_mako[3])
p_gene_apoe <- se.rdist[[]] %>% 
  bind_cols(FetchData(se.rdist, vars = genes_apoe)) %>% 
  filter(r_dist_F14_C0 < dist_cutoff_plot) %>% 
  pivot_longer(all_of(genes_apoe), names_to = "gene", values_to = "expression") %>% 
  mutate_at(.vars="gene", ~ factor(., levels = genes_apoe)) %>% 
  ggplot(aes(r_dist_F14_C0, expression, color = gene)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), linewidth=1, fill="grey80") +
  geom_vline(aes(xintercept = 0, color = "border"), 
             linetype = "dashed", linewidth=0.5, color="black") +
  facet_wrap(~gene, scales = "free_y", ncol = 1, strip.position = "right") +
  scale_color_manual(values = cols_apoe) +
  theme_classic() +
  theme(strip.text.y.right = element_text(angle=0, hjust=0), 
        strip.background = element_rect(color = NA), 
        legend.position = "none",
        axis.text.x = element_text(angle=45, vjust = 1, hjust=1), 
        panel.grid.major.x = element_line(linewidth=0.5, color="grey90")
  );p_gene_apoe
pdf(file.path(DIR_FIG_OUT, "hs_visium_A-NMF30-F14hi-C0_distance_gene_apoe.pdf"), width=2.5, height=5);p_gene_apoe;dev.off()


###### Cell type cor ###### 
dist_cutoff_plot <- 2e3
cell_selected <- c("KRT5..KRT17.", 
                   "Transitional.AT2", "AT2", "SCGB3A2.", 
                   "PLIN2..Fibroblasts", "Mesothelial.Cells")
# cell_selected <- cell_cor_df[cell_cor_df$cor>0,]$cell_type %>% unique()
cell_cor_df$cell_type2 <- cell_cor_df$cell_type %>% 
  as.character() %>% 
  gsub(pattern = "\\.\\.", replacement = "+") %>% 
  gsub(pattern = "\\.$", replacement = "+") %>% 
  gsub(pattern = "\\.", replacement = " ")
cell_cor_df$cell_type2[cell_cor_df$cell_type2 == "KRT5+KRT17+"] <- "KRT5-/KRT17+"

# Plot selected cell types
cell_cor_df_p2 <- cell_cor_df %>% filter(cell_type %in% cell_selected)

plot_df <- se.rdist[[]] %>% 
  filter(r_dist_F14_C0 < dist_cutoff_plot) %>% 
  pivot_longer(all_of(paste0("c2l_", cell_selected)), 
               names_to = "cell_type", values_to = "density") %>% 
  mutate_at(.vars="cell_type", ~ factor(., levels = paste0("c2l_", cell_selected)))
plot_df$cell_type <- gsub("c2l_", "", plot_df$cell_type)

# Retrieve distance peaks to define plot order
gam_dist_max <- plot_df %>% 
  group_by(cell_type) %>% 
  mutate(gam_fit_vals = mgcv::gam(formula = density ~ s(r_dist_F14_C0), method="REML")$fitted.values)

gam_dist_max_cells <- gam_dist_max %>% 
  group_by(cell_type) %>% 
  summarise(max_dens_dist = r_dist_F14_C0[which.max(gam_fit_vals)] %>% round()) %>% 
  arrange(max_dens_dist)

cell_order <- c("KRT5..KRT17.", gam_dist_max_cells$cell_type[-1])
plot_df$cell_type <- factor(plot_df$cell_type, levels=cell_order)
cell_selected2 <- cell_selected[-c(5,6)]
plot_df2 <- subset(plot_df, cell_type %in% cell_selected2)


p_cell <- ggplot(plot_df2, aes(x=r_dist_F14_C0, y=density, color = cell_type)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), 
              linewidth=0.8, fill="grey70") + # color="#830051"
  geom_vline(aes(xintercept = 0, color = "border"), 
             linetype = "dashed", linewidth=0.5, color="black") +
  facet_wrap(~cell_type, scales = "free_y", ncol = 1, strip.position = "right") +
  scale_color_manual(values = c("#830051", col_scale_spec)) +
  labs(x="Distance (µm) from F14hiC0", y="Cell type density") +
  theme_classic() +
  theme(strip.text.y.right = element_text(angle=0, hjust=0), 
        strip.background = element_rect(color = NA), 
        legend.position = "none",
        axis.text.x = element_text(angle=45, vjust = 1, hjust=1), 
        panel.grid.major.x = element_line(linewidth=0.5, color="grey90")
  );p_cell

pdf(file.path(DIR_FIG_OUT, "hs_visium_A-NMF30-F14hi-C0_distance_cell_selected.pdf"), width=3, height=3);p_cell;dev.off()  # !Main Figure


# Plot per donor
plot_df2$subject_alias <- factor(plot_df2$subject_alias, levels = paste0("IPF_", 1:4))
p_cell_donor <- ggplot(plot_df2, aes(x=r_dist_F14_C0, y=density, color = cell_type)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), 
              linewidth=0.8, fill="grey70") + # color="#830051"
  geom_vline(aes(xintercept = 0, color = "border"), 
             linetype = "dashed", linewidth=0.5, color="black") +
  # facet_grid(vars(cell_type), vars(subject_alias), scales = "free_y") +
  facet_wrap(~cell_type + subject_alias, scales = "free_y", ncol = 3, strip.position = "top") +
  scale_color_manual(values = c("#830051", col_scale_spec)) +
  labs(x="Distance (µm) from F14hiC0", y="Cell type density") +
  theme_classic() +
  theme(
    # strip.text.y.right = element_text(angle=0, hjust=0), 
    strip.background = element_rect(color = NA), 
    strip.text.x = element_text(margin = margin(2, 0, 2, 0)),
    legend.position = "none",
    axis.text.x = element_text(angle=45, vjust = 1, hjust=1), 
    panel.grid.major.x = element_line(linewidth=0.5, color="grey90")
  );p_cell_donor

pdf(file.path(DIR_FIG_OUT, "hs_visium_A-NMF30-F14hi-C0_distance_cell_selected_per_donor.pdf"), width=6, height=5.5);p_cell_donor;dev.off()


# Plot all cell types
cell_selected <- cell_cor_df$cell_type %>% unique()

plot_df <- se.rdist[[]] %>% 
  filter(r_dist_F14_C0 < dist_cutoff_plot) %>% 
  pivot_longer(all_of(paste0("c2l_", cell_selected)), 
               names_to = "cell_type", values_to = "density") %>% 
  mutate_at(.vars="cell_type", ~ factor(., levels = paste0("c2l_", cell_selected)))
plot_df$cell_type <- gsub("c2l_", "", plot_df$cell_type)

p_cell_all <- ggplot(plot_df, aes(x=r_dist_F14_C0, y=density)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), 
              linewidth=0.8, fill="grey70", color="#830051") + # color="#830051"
  geom_vline(aes(xintercept = 0, color = "border"), 
             linetype = "dashed", linewidth=0.5, color="black") +
  facet_wrap(~cell_type, scales = "free_y", ncol = 4, strip.position = "top") +
  scale_color_manual(values = c("#830051", col_scale_spec)) +
  labs(x="Distance (µm) from F14hiC0", y="Cell type density") +
  theme_classic() +
  theme(strip.text.y.right = element_text(angle=0, hjust=0), 
        strip.background = element_rect(color = NA), 
        legend.position = "none",
        # axis.text.x = element_text(angle=45, vjust = 1, hjust=1), 
        panel.grid.major.x = element_line(linewidth=0.5, color="grey90")
  );p_cell_all

pdf(file.path(DIR_FIG_OUT, "hs_visium_A-NMF30-F14hi-C0_distance_cell_all.pdf"), width=8, height=8);p_cell_all;dev.off()




