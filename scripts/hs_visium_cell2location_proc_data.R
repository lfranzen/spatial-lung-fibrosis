#' [hs_visium_cell2location_proc_data.R]
#'
#' Run, analyse, and plot NMF results 
#'
#'
#' Feb 2023, L. Franzén [lovisa.franzen@scilifelab.se]

#### Set up ####
##### Define params. ####
set.seed(1)
SPECIES <- "human"
DIR_ROOT <- getwd()
DIR_DATA <- file.path(DIR_ROOT, "data", SPECIES, "sc_deconvolution_habermann")
DIR_RES <- file.path(DIR_ROOT, "results", SPECIES)
DIR_FIG_C2L <- file.path(DIR_RES, "figures", "cell2location_habermann")
DIR_OBJ_C2L <- file.path(DIR_RES, "objects", "cell2location_habermann")
dir.create(DIR_FIG_C2L);dir.create(DIR_OBJ_C2L)
fig_res <- 500

##### Load libs ####
library(STutility)
library(harmony)
library(RColorBrewer)
library(patchwork)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(ggpmisc)


##### Other ####
source(file.path(DIR_ROOT, "scripts", "custom_functions.R"))
source(file.path(DIR_ROOT, "scripts", "custom_colors.R"))
theme_custom <- theme(axis.title.x = element_blank())


##### Read se object ####
fname <- paste0("hs_visium_preproc_A_se_obj_nmf.rds")
se <- readRDS(file = file.path(DIR_RES, "objects", fname))


##### Read tables ####
metadata <- read.table(file.path(DIR_DATA, "../visium/hs_visium_metadata.tsv"), sep = "\t", header = T)
rownames(metadata) <- metadata$sample_id
metadata$sample_n <- 1:nrow(metadata)


##### Read cell annotatoins ####
cell_anno <- read.csv(file.path(DIR_ROOT, "data", "misc", "habermann_cell_type_groups.csv"), sep = ";")


##### Read cell2location csv ####
c2l_list <- setNames(
  lapply(metadata$sample_id, function(id){
    message(id)
    c2l_df <- read.csv(file = file.path(DIR_DATA, paste0(id, "_spot_cell_abundances_5pc.csv")), header = T)
    c2l_df$sample_id <- id
    rownames(c2l_df) <- paste0(gsub(paste0(id, "_"), "", c2l_df$spot_id), "_", subset(metadata, sample_id == id)$sample_n)
    colnames(c2l_df) <- gsub("q05cell_abundance_w_sf_", "", colnames(c2l_df))
    return(c2l_df)
  }),
  nm = metadata$sample_id
)

c2l_all <- bind_rows(c2l_list)
c2l_all <- c2l_all %>% select(c(-sample_id, -spot_id))
colnames(c2l_all) <- paste0("c2l_", colnames(c2l_all))

write.csv(c2l_all, file = file.path(DIR_DATA, "compiled_all_samples_cell_abundances.csv"), row.names = T)
c2l_all <- read.csv(file.path(DIR_DATA, "compiled_all_samples_cell_abundances.csv"), row.names = 1)

cell_anno$c2l_colnames <- colnames(c2l_all)


##### Add metadata to se obj ####
se <- AddMetaData(se, c2l_all)


##### Define variables ####
factor_names <- paste0("factor_", 1:30)
c2l_names <- colnames(c2l_all)
new_c2l_names <- gsub("c2l_", "", gsub("[.]$", "", gsub("..", "_", c2l_names, fixed = TRUE)))
new_c2l_names_filtered <- new_c2l_names


##### Violin & Spatial plots ####
png(file.path(DIR_FIG_C2L, paste0("hs_visium_c2l_res_violin_group.png")), width = 30*fig_res, height = 16*fig_res, res = fig_res)
VlnPlot(se, features = c2l_names, pt.size = 0, ncol = 8, group.by = "subject_alias", cols = cols_donor) & theme_custom
dev.off()

png(file.path(DIR_FIG_C2L, paste0("hs_visium_c2l_res_violin_histpath.png")), width = 30*fig_res, height = 16*fig_res, res = fig_res)
VlnPlot(se, features = c2l_names, ncol = 8, group.by = "annotation", pt.size = 0, cols = cols_d3_20) & theme_custom
dev.off()

png(file.path(DIR_FIG_C2L, paste0("hs_visium_c2l_res_spatial_KRT5KRT17.png")), width = 14*fig_res, height = 14*fig_res, res = fig_res)
ST.FeaturePlot(se, features = grep("KRT5", c2l_names, value = T),
               ncol = 5, cols = col_scale_mako, pt.border = F, 
               min.cutoff = 1, show.sb = F,
               label.by = "sample_name")
dev.off()

png(file.path(DIR_FIG_C2L, paste0("hs_visium_c2l_res_spatial_HAS1Fibroblasts.png")), width = 14*fig_res, height = 14*fig_res, res = fig_res)
ST.FeaturePlot(se, features = grep("HAS1", c2l_names, value = T),
               ncol = 5, cols = col_scale_mako, pt.border = F, 
               min.cutoff = 1, show.sb = F,
               label.by = "sample_name")
dev.off()



##### Cell~X correlation #####
factor_cell_cor_data <- bind_cols(se@meta.data %>% select(sample_name, subject_alias, condition, matches(c2l_names)),
                                  se@reductions$NMF@cell.embeddings[, factor_names])

colnames(factor_cell_cor_data)[colnames(factor_cell_cor_data) %in% c2l_names] <- new_c2l_names
# write.csv(factor_cell_cor_data, file.path(DIR_RES, "objects", "A_NMF30", "hs_visium_A_nmf_1-30_factorweight_cell2location_habermann_metadata.csv"))


###### Cell~Cell correlation ######
cell_cor_data <- se@meta.data %>% select(sample_name, subject_alias, condition, matches(c2l_names))
colnames(cell_cor_data)[colnames(cell_cor_data) %in% c2l_names] <- new_c2l_names

#' Filter low values data?
cell_cor_data <- cell_cor_data %>% mutate_at(vars(new_c2l_names), ~if_else(.<0.1, 0, .))

groups_use <- c("All", "Healthy", "IPF_1", "IPF_2", "IPF_3", "IPF_4")
cell_cor_data_group <- setNames(lapply(groups_use, function(g){
  if(g=="All"){
    cell_cor_data_subset <- cell_cor_data
  } else if(g=="Healthy"){
    cell_cor_data_subset <- subset(cell_cor_data, condition == "control")
  } else {
    cell_cor_data_subset <- subset(cell_cor_data, subject_alias == g)
  }
  cell_density_hm <- cor(cell_cor_data_subset[,new_c2l_names_filtered])
  diag(cell_density_hm) <- 0
  cell_density_hm <- as.data.frame(cell_density_hm)
}), 
nm = groups_use)

saveRDS(cell_cor_data_group, file.path(DIR_OBJ_C2L, "hs_visium_A_nmf_1-30_cell_cor_data_group.rds"))


# Plot pheatmap
pal_length <- 11
ph_colors <- RColorBrewer::brewer.pal(pal_length, "RdBu") %>% rev()
abs_max_val <- lapply(cell_cor_data_group, max) %>% unlist() %>% max()

pdf(file = file.path(DIR_FIG_C2L, "hs_visium_c2l_res_cell_type_cor_heatmap_per_subject.pdf"), width = 10, height = 10, useDingbats = F)
for(g in groups_use){
  cell_density_hm <- cell_cor_data_group[[g]]
  
  if(g=="all"){abs_max_val <- max(cell_density_hm)}
  
  ph_breaks <- c(seq(-abs_max_val, 0, length.out=ceiling(pal_length/2) + 1), 
                 seq(abs_max_val/pal_length, abs_max_val, length.out=floor(pal_length/2)))
  
  print(
    pheatmap::pheatmap(cell_density_hm[, new_c2l_names_filtered], 
                       cellwidth = 10, 
                       cellheight = 10, 
                       color = ph_colors, 
                       breaks = ph_breaks, 
                       main = g)
  )
}
dev.off()


#' Plot only KRT5-/KRT17+ group
cell_cor_data_group_abba <- do.call(cbind, cell_cor_data_group[2:6]) %>% select(contains("KRT17"))
cell_cor_data_group_abba %>% head()

cell_cor_data_group_abba <- cell_cor_data_group_abba[rownames(cell_cor_data_group_abba)!="KRT5_KRT17",-1]
# pal_length <- 11
# ph_colors <- RColorBrewer::brewer.pal(pal_length, "RdBu") %>% rev()

pal_length <- length(col_scale_div_custom2)
ph_colors <- col_scale_div_custom2
abs_max_val <- max(cell_cor_data_group_abba)
ph_breaks <- c(seq(-abs_max_val, 0, length.out=ceiling(pal_length/2) + 1), 
               seq(abs_max_val/pal_length, abs_max_val, length.out=floor(pal_length/2)))

pdf(file = file.path(DIR_FIG_C2L, paste0("hs_visium_c2l_res_cell_type_cor_heatmap_per_subject_KRT5KRT17_onlyIPF.pdf")), 
    width = 6, height = 6, useDingbats = F)
pheatmap::pheatmap(cell_cor_data_group_abba, 
                   cellwidth = 20, 
                   cellheight = 10, 
                   treeheight_row = 10,  
                   treeheight_col = 4,
                   # color = col_scale_div_custom2, 
                   color = ph_colors,
                   breaks = ph_breaks
                   )
dev.off()


###### Cell~Factor correlation ######
groups_use <- c("All", "Healthy", "IPF_1", "IPF_2", "IPF_3", "IPF_4")

factor_cell_cor_data_group <- setNames(lapply(groups_use, function(g){
  if (g=="All") {
    factor_cell_cor_subset <- factor_cell_cor_data
  } else if (g=="Healthy") {
    factor_cell_cor_subset <- subset(factor_cell_cor_data, condition == "control")
  } else {
    factor_cell_cor_subset <- subset(factor_cell_cor_data, subject_alias == g)
  }
  factor_cell_density_hm <- cor(factor_cell_cor_subset[,c(new_c2l_names_filtered, factor_names)])
  diag(factor_cell_density_hm) <- 0
  factor_cell_density_hm <- as.data.frame(factor_cell_density_hm)
  factor_cell_density_hm <- factor_cell_density_hm[factor_names, new_c2l_names_filtered]
  }),
  nm = groups_use)

saveRDS(factor_cell_cor_data_group, file.path(DIR_OBJ_C2L, "hs_visium_A_nmf_1-30_factor_cell_cor_data_group.rds"))


# Plot pheatmap
pal_length <- 11
ph_colors <- RColorBrewer::brewer.pal(pal_length, "RdBu") %>% rev()
abs_max_val <- lapply(factor_cell_cor_data_group, max) %>% unlist() %>% max()

pdf(file = file.path(DIR_FIG_C2L, "hs_visium_c2l_res_cell_type_A-NMF30_factor_cor_heatmap_per_subject.pdf"), width = 10, height = 10, useDingbats = F)
for(g in groups_use){
  factor_cell_density_hm <- factor_cell_cor_data_group[[g]]

  if(g=="all"){abs_max_val <- max(cell_density_hm)}
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


#' Plot only F14
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

#' Plot only KRT5-/KRT17+ group
factor_cell_cor_data_group_abba <- do.call(cbind, factor_cell_cor_data_group[2:6]) %>% select(contains("KRT17"))
factor_cell_cor_data_group_abba %>% head()

pal_length <- 11
ph_colors <- RColorBrewer::brewer.pal(pal_length, "RdBu") %>% rev()
abs_max_val <- max(factor_cell_cor_data_group_abba)

pdf(file = file.path(DIR_FIG_C2L, paste0("hs_visium_c2l_res_cell_type_A-NMF30_factor_cor_heatmap_per_subject_KRT5KRT17.pdf")), 
    width = 6, height = 6, useDingbats = F)
pheatmap::pheatmap(factor_cell_cor_data_group_abba, 
                   cellwidth = 20, 
                   cellheight = 10, 
                   treeheight_row = 10,  
                   treeheight_col = 4,
                   color = ph_colors, 
                   breaks = ph_breaks)
dev.off()



#### Count cells per sample ####
c2l_mdata <- se@meta.data
c2l_mdata$group_subject <- ifelse(c2l_mdata$condition == "control", "HC", c2l_mdata$subject_alias)

cell_select_list <- c2l_names %>% sort()

cell_stats_list <- lapply(cell_select_list, function(cell_select){
  c2l_mdata$cell_oi <- c2l_mdata[, cell_select]
  c2l_cell_stats <- c2l_mdata %>% 
    group_by(sample_name, group_subject, condition) %>% 
    summarise(.groups = "keep", 
              sum_density_cells = sum(cell_oi),
              avg_density_cells = mean(cell_oi)
    )
  c2l_cell_stats$cell_oi <- cell_select
  return(c2l_cell_stats)
}) %>% setNames(nm = cell_select_list)
cell_stats_all <- bind_rows(cell_stats_list)
cell_stats_all$cell_oi <- gsub("c2l_", "", cell_stats_all$cell_oi)

p_stats_list <- lapply(cell_stats_list, function(c2l_cell_stats){
  cell_select <- c2l_cell_stats$cell_oi %>% unique()
  p1 <- ggplot(c2l_cell_stats, aes(x=group_subject, y=sum_density_cells, color=condition)) +
    geom_point(size=1) +
    geom_boxplot(width=0.4) +
    scale_color_manual(values = cols_cond) +
    labs(y="Sum of cell density", x="", title="Cell density sum per sample") +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5), 
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  p2 <- ggplot(c2l_cell_stats, aes(x=group_subject, y=avg_density_cells, color=condition)) +
    geom_point(size=1) +
    geom_boxplot(width=0.4) +
    scale_color_manual(values = cols_cond) +
    labs(y="Mean of cell density", x="", title="Cell density mean per sample") +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5), 
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  (p1|p2) + plot_annotation(title = gsub("c2l_", "", cell_select), 
                            theme = theme(plot.title = element_text(hjust=0.5, face="bold")))
})

pdf(file = file.path(DIR_FIG_C2L, "hs_visium_c2l_res_cell_type_count_per_subject.pdf"), width = 6, height = 4, useDingbats = F)
for(p in p_stats_list){print(p)}
dev.off()

pdf(file = file.path(DIR_FIG_C2L, "hs_visium_c2l_res_cell_type_count_per_subject_grid.pdf"), width = 8, height = 14, useDingbats = F)
ggplot(cell_stats_all, aes(x=group_subject, y=avg_density_cells, color=condition)) +
  geom_point(size=1) +
  geom_boxplot(width=0.7) +
  scale_color_manual(values = cols_cond) +
  facet_wrap(~cell_oi, ncol = 4, scales = "free_y") +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5), 
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

#### HE Spatial #### 
se_he <- SubsetSTData(se, sample_name %in% c("HC_4.TH040.B0.2", 
                                             "IPF_1.TD012.B2.2", 
                                             "IPF_3.TD031.B1.2", 
                                             "IPF_4.TD041B.B2.1"))
se_he$sample_name %>% unique()
se_he <- LoadImages(se_hes, xdim = 1e3)

# ImagePlot(se_he)

##### KRT5-/KR17+ #####
p <- FeatureOverlay(se_he, sampleids = 1:4, 
                    features = grep("KRT5", c2l_names, value = T), 
                    label.by = "sample_name",
                    # value.scale = "all", 
                    min.cutoff = 0, 
                    max.cutoff = 6,
                    pt.size = 1.8,
                    ncols = 2, cols = rev(col_scale_mako), 
                    add.alpha = T) &
  theme(aspect.ratio = 1)
p <- ggrastr::rasterize(p, layers = "Point", dpi = fig_res)

png(file.path(DIR_FIG_C2L, paste0("hs_visium_c2l_res_spatial_KRT5KRT17_he.png")), 
    width = 14*fig_res, height = 14*fig_res, res = fig_res);p;dev.off()


##### HAS1+ Fibroblasts #####
p <- FeatureOverlay(se_he, sampleids = 1:4, 
                    features = grep("HAS1", c2l_names, value = T), 
                    label.by = "sample_name", 
                    value.scale = "all", 
                    min.cutoff = 0, 
                    max.cutoff = 5,
                    pt.size = 1.8,
                    ncols = 2, cols = rev(col_scale_mako), 
                    add.alpha = T) &
  theme(aspect.ratio = 1)
p <- ggrastr::rasterize(p, layers = "Point", dpi = fig_res)

png(file.path(DIR_FIG_C2L, paste0("hs_visium_c2l_res_spatial_HAS1Fibroblasts_he.png")), 
    width = 14*fig_res, height = 14*fig_res, res = fig_res);p;dev.off()

