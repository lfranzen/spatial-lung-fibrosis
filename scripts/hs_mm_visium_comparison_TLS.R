#' [hs_mm_visium_comparison_TLS.R]
#'
#' Comparison of regions of TLS in both IPF and BLM Visium data 
#'
#'
#' Mar 2023, L. Franzén [lovisa.franzen@scilifelab.se]

#### Set up ####
##### Define params. ####
SPECIES <- "mouse"
SPECIES <- "human"
DIR_ROOT <- getwd()
DIR_DATA <- file.path(DIR_ROOT, "data")
DIR_RES <- file.path(DIR_ROOT, "results", "translational")
DIR_FIG <- file.path(DIR_RES, "figures")
DIR_FIG_OUT <- file.path(DIR_FIG, "TLS")
DIR_OBJ_OUT <- file.path(DIR_RES, "objects", "TLS")
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


#### Read se obj ####
# Human NMF-B object list
hs_donor_select <- c("IPF_1", "IPF_2", "IPF_3")
tls_sample_select <- list(IPF_1 = "IPF_1.TD012.G2.2", IPF_2 = "IPF_2.TD022.G3.1", IPF_3 = "IPF_3.TD031.G1.2")

fname <- paste0("hs_visium_preproc_B_se_obj_list.rds")
hs_se_list <- readRDS(file = file.path(DIR_RES, "..", "human", "objects", fname))
hs_se_list <- hs_se_list[hs_donor_select]


# Mouse d7
# fname <- paste0("hs_visium_preproc_B_se_obj_list.rds")
# mm_se_d7 <- readRDS(file = file.path(DIR_RES, "..", "human", "objects", fname))
se_mmd7 <- readRDS(file = file.path(DIR_RES, "..", "mouse", "objects", "mm_visium_preproc_se_obj_subset_d7_nmf30_c2l.rds"))

# Mouse d21
se_mmd21 <- readRDS(file = file.path(DIR_RES, "..", "mouse", "objects", "mm_visium_preproc_se_obj_subset_d21_nmf30_c2l.rds"))


#### TLS-associated NMF factors ####
tls_factors <- list(IPF_1 = "factor_10", IPF_2 = "factor_15", IPF_3 = "factor_17", mm_d7 = "factor_11", mm_d21 = "factor_12")

##### IPF ##### 
#' Select based on NMF30-B
#' IPF_1: F10 (?); IPF_2: F15; IPF_3: F17
hs_tls_factor_genes_list <- list(
  IPF_1.NMF_10 = read_excel(path = file.path(DIR_RES, "..", "human", "objects", "B_NMF30", "hs_visium_B_nmf_1-30_top100_gene_loadings.xlsx"), sheet = "IPF_1") %>% filter(factor == 10),
  IPF_2.NMF_15 = read_excel(path = file.path(DIR_RES, "..", "human", "objects", "B_NMF30", "hs_visium_B_nmf_1-30_top100_gene_loadings.xlsx"), sheet = "IPF_2") %>% filter(factor == 15),
  IPF_3.NMF_17 = read_excel(path = file.path(DIR_RES, "..", "human", "objects", "B_NMF30", "hs_visium_B_nmf_1-30_top100_gene_loadings.xlsx"), sheet = "IPF_3") %>% filter(factor == 17)
)
hs_tls_factor_genes <- bind_rows(hs_tls_factor_genes_list, .id = "group")
# hs_tls_factor_genes$gene_loading_scaled_original <- hs_tls_factor_genes$gene_loading_scaled
hs_tls_factor_genes <- hs_tls_factor_genes %>% 
  group_by(group) %>% 
  mutate(gene_loading_scaled = Scale01(gene_loading))


##### BLM d7 ##### 
#' Select based on d7-NMF30
#' F11
mm_d7_tls_factor_genes <- read_excel(path = file.path(DIR_RES, "..", "mouse", "objects", "NMF30_d7", "mm_visium_nmf30_d7_top100_gene_loadings.xlsx"), 
                                     sheet = 11)
mm_d7_tls_factor_genes$gene_loading_scaled <- Scale01(mm_d7_tls_factor_genes$gene_loading)
mm_d7_tls_factor_genes$group <- "mm_d7.NMF_F11"

##### BLM d21 ##### 
#' Select based on d21-NMF30
#' F12
mm_d21_tls_factor_genes <- read_excel(path = file.path(DIR_RES, "..", "mouse", "objects", "NMF30_d21", "mm_visium_nmf30_d21_top100_gene_loadings.xlsx"), 
                                      sheet = 12)
mm_d21_tls_factor_genes$gene_loading_scaled <- Scale01(mm_d21_tls_factor_genes$gene_loading)
mm_d21_tls_factor_genes$group <- "mm_d21.NMF_F12"


##### Join top gene contributions ##### 
hs_mm_tls_factor_genes <- bind_rows(list(hs_tls_factor_genes, mm_d7_tls_factor_genes, mm_d21_tls_factor_genes))

# shared genes:
hs_mm_tls_factor_genes$shared_gene_hs <- ifelse(hs_mm_tls_factor_genes$gene %in% gene_conv_df$symbol_hs, 
                                             hs_mm_tls_factor_genes$gene,
                                             ifelse(hs_mm_tls_factor_genes$gene %in% gene_conv_df$symbol_mm,
                                                    gene_conv_mm[hs_mm_tls_factor_genes$gene],
                                                    NA)
                                             )
hs_mm_tls_factor_genes$shared_gene_mm <- gene_conv[hs_mm_tls_factor_genes$shared_gene_hs]

hs_mm_tls_factor_genes$shared_gene <- paste0(hs_mm_tls_factor_genes$shared_gene_hs, ":", hs_mm_tls_factor_genes$shared_gene_hs_mm)

##### Fig 6b: Plot factor gene contribution heatmap ##### 
phm_df <- hs_mm_tls_factor_genes %>% 
  filter(!is.na(shared_gene_hs)) %>% 
  select(group, shared_gene, gene_loading_scaled) %>% 
  pivot_wider(names_from = group, values_from = gene_loading_scaled) %>% 
  mutate(across(where(is.numeric), replace_na, 0)) %>% 
  as.data.frame()
rownames(phm_df) <- phm_df$shared_gene
phm_df <- phm_df[,-1]


#' Filt based on gene expr in at least two IPF subjects
#' and total row sum cutoff
phm_df_filt <- phm_df[rowSums(phm_df)>0.5,] # filter


#' View and plot
dim(phm_df_filt);head(phm_df_filt)


pdf(file.path(DIR_FIG_OUT, "hs_mm_visium_tls_factor_gene_heatmap_cutoff05.pdf"),  # !Manuscript Main figure 6b
    width = 4, height = 10)
pheatmap((phm_df_filt), 
         treeheight_row = 10, 
         treeheight_col = 4, 
         cutree_cols = 2,
         # cutree_rows = 2,
         color = c("grey85", col_scale_mako), 
         border_color = NA,
         cellwidth = 10,
         cellheight = 10,
         angle_col = 90
)
dev.off()


write.csv(x = phm_df, file = file.path(DIR_OBJ_OUT, "hs_mm_visium_tls_factor_gene_heatmap_unfiltered.csv"), row.names = T)



##### Fig 6a: Plot spatial factor overlay ##### 
# Plot factor enrichment in all samples
p_list_hs <- lapply(hs_donor_select, function(s){
  p <- ST.FeaturePlot(hs_se_list[[s]], 
                           features = tls_factors[[s]], 
                           label.by = "sample_name",
                           ncol = 3, cols = c("grey85", col_scale_mako),
                           min.cutoff = 0.25
  ) &
    theme(aspect.ratio = 1)
  return(p)
})

p_list_mm <- list(
  ST.FeaturePlot(se_mmd7,
               features = tls_factors[["mm_d7"]], 
               label.by = "sample_name",
               ncol = 3, cols = c("grey85", col_scale_mako),
               min.cutoff = 0.25) &
    theme(aspect.ratio = 1),
  ST.FeaturePlot(se_mmd21,
                 features = tls_factors[["mm_d21"]], 
                 label.by = "sample_name",
                 ncol = 3, cols = c("grey85", col_scale_mako),
                 min.cutoff = 0.25) &
    theme(aspect.ratio = 1)
  )

p_list <- c(p_list_hs, p_list_mm)

pdf(file.path(DIR_FIG_OUT, "hs_mm_visium_tls_factor_activity_spatial.pdf"), 
    width = 8, height = 12)
for(i in 1:length(p_list)){
  print(p_list[[i]])
}
dev.off()

# Plot factor enrichment in selected samples
se_selected_list <- list(SubsetSTData(hs_se_list$IPF_1, sample_name == "IPF_1.TD011.G1.1"),
                         SubsetSTData(hs_se_list$IPF_2, sample_name == "IPF_2.TD022.G3.1"),
                         SubsetSTData(hs_se_list$IPF_3, sample_name == "IPF_3.TD032B.G3.2"),
                         SubsetSTData(se_mmd7, sample_name == "d7_bleo_1"),
                         SubsetSTData(se_mmd21, sample_name == "d21_bleo_4a")
                         )
names(se_selected_list) <- c("IPF_1", "IPF_2", "IPF_3", "mm_d7", "mm_d21")

p_list <- lapply(names(se_selected_list), function(s){
  f_plot <- tls_factors[[s]]
  p <- ST.FeaturePlot(se_selected_list[[s]], 
                      features = f_plot, 
                      min.cutoff = 0.25, 
                      label.by = "sample_name",
                      cols = c("grey85", col_scale_mako)) &
    theme(aspect.ratio = 1, 
          plot.title = element_text(hjust = 0.5), 
          legend.position = "bottom")
  return(p)
})

pdf(file.path(DIR_FIG_OUT, "hs_mm_visium_tls_factor_activity_selected_spatial.pdf"), 
    width = 10, height = 4)
wrap_plots(p_list, nrow = 1)
dev.off()

# Plot factor enrichment in selected samples, H&E overlay
for(s in names(se_selected_list)){
  se_selected_list[[s]] <- LoadImages(se_selected_list[[s]], xdim = 2e3)
}

pdf(file.path(DIR_FIG_OUT, "hs_mm_visium_tls_factor_activity_selected_spatial_HE.pdf"),   # !Manuscript Main figure 6a
    width = 8, height = 8, useDingbats = F)
for(s in names(se_selected_list)){
  p <- FeatureOverlay(se_selected_list[[s]], features = tls_factors[[s]], 
                 label.by = "sample_name",
                 min.cutoff = 0.25, 
                 max.cutoff = 'q99',
                 cols=rev(col_scale_mako),
                 add.alpha = T) &
    theme(aspect.ratio = 1, 
          plot.title = element_text(hjust = 0.5), 
          legend.position = "bottom")
  print(p)
}
dev.off()


#' TLS-proximal factor in ms d21 (spatial)
pdf(file.path(DIR_FIG_OUT, "hs_mm_visium_tls_factor_activity_mmNMFd21_F12_F13_spatial_HE.pdf"),   # !Manuscript Main figure 6a
    width = 8, height = 8, useDingbats = F)
FeatureOverlay(se_selected_list[["mm_d21"]], features = "factor_12",
               label.by = "sample_name",
               min.cutoff = 0.25,
               max.cutoff = 'q99',
               cols=rev(col_scale_mako),
               add.alpha = T) &
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")
FeatureOverlay(se_selected_list[["mm_d21"]], features = "factor_13",
               label.by = "sample_name",
               min.cutoff = 0.25,
               max.cutoff = 'q99',
               cols=rev(col_scale_mako),
               add.alpha = T) &
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")
ST.FeaturePlot(se_selected_list[["mm_d21"]], 
               features = c("factor_12", "factor_13"), 
               pt.size = 2, min.cutoff = 0.25,
               blend = T, dark.theme = T, channels.use = c("red", "blue")) &
  theme(aspect.ratio = 1, 
        plot.title = element_text(hjust = 0.5), 
        legend.position = "bottom")
dev.off()


# F13 gene loading plot
mmd21_f13_loading <- data.frame(loading = se_mmd21@reductions$NMF@feature.loadings[, 13])
mmd21_f13_loading <- mmd21_f13_loading %>% arrange(desc(loading))
mmd21_f13_loading$gene <- factor(rownames(mmd21_f13_loading), levels = rownames(mmd21_f13_loading))


pdf(file.path(DIR_FIG_OUT, "hs_mm_visium_tls_factor_activity_mmNMFd21_F13_gene_loading_top15.pdf"),   # !Manuscript Extended Data figure 5c
    width = 2, height = 3, useDingbats = F)
ggplot(head(mmd21_f13_loading, n=15), aes(x=reorder(gene, loading), y=loading)) +
  geom_point(size=2, color="black") +
  geom_segment(aes(x=reorder(gene, loading), xend=reorder(gene, loading), y=0, yend=loading), 
               linewidth=0.25, linetype="dashed", color="black") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text.y = element_text(face = "italic", color="black"),
        axis.text.x = element_text(color="black"),
        axis.title.y = element_blank())
dev.off()



##### Factor gene pathway enrichment (gProfiler) ##### 
p1 <- ST.FeaturePlot(se_mmd7,
               features = tls_factors[["mm_d7"]], 
               label.by = "sample_name",
               ncol = 3, cols = c("grey", col_scale_mako[6:10]),
               min.cutoff = 'q99') &
  theme(aspect.ratio = 1)
p2 <- ST.FeaturePlot(se_mmd21,
                     features = tls_factors[["mm_d21"]], 
                     label.by = "sample_name",
                     ncol = 3, cols = c("grey", col_scale_mako[6:10]),
                     min.cutoff = 'q99') &
  theme(aspect.ratio = 1)

p1|p2

##### TLS-factor cut-off #####
hs_se_list <- setNames(
  lapply(hs_donor_select, function(s){
    f <- tls_factors[[s]]
    se <- hs_se_list[[s]]
    se_nmf_emb <- se@reductions$NMF@cell.embeddings
    f_cutoff <- as.numeric(quantile(se_nmf_emb[, f], c(.99)))
    f_metadata_add <- ifelse(se_nmf_emb[,f]>f_cutoff,"high", "low")
    se <- AddMetaData(se, metadata = f_metadata_add, col.name = "tls_factor")
    return(se)
  }), 
  nm = hs_donor_select)

mm_se_list <- setNames(
  lapply(c("mm_d7", "mm_d21"), function(mm_d){
    if (mm_d == "mm_d7") {
      se <- se_mmd7
    } else if(mm_d == "mm_d21"){
      se <- se_mmd21
    }
    f <- tls_factors[[mm_d]]
    se_nmf_emb <- se@reductions$NMF@cell.embeddings
    f_cutoff <- as.numeric(quantile(se_nmf_emb[, f], c(.99)))
    f_metadata_add <- ifelse(se_nmf_emb[,f]>f_cutoff,"high", "low")
    se <- AddMetaData(se, metadata = f_metadata_add, col.name = "tls_factor")
    return(se)
  }), 
  nm=c("mm_d7", "mm_d21"))


##### Factor - cell density enrichment ##### 

# hs_c2l_all <- read.csv(file.path(DIR_DATA, "human", "sc_deconvolution_habermann", "compiled_all_samples_cell_abundances.csv"), row.names = 1)
# hs_c2l_colnames <- colnames(hs_c2l_all)
# 
# factor_cell_cor_data <- lapply(names(hs_se_list), function(s){
#   message(s)
#   se <- hs_se_list[[s]]
#   se <- AddMetaData(se, metadata = hs_c2l_all[colnames(se), ])
#   f_dat <- se@reductions$NMF@cell.embeddings
#   colnames(f_dat) <- paste0(s, ".", gsub("factor_", "F", colnames(f_dat)))
#   cell_dat <- se@meta.data %>% select(sample_name, tissue_alias, subject_alias, condition, species, colnames(hs_c2l_all))
#   nmf_cell_dat <- bind_cols(cell_dat, f_dat)
#   return(nmf_cell_dat)
# }) %>% setNames(nm=names(hs_se_list))
#
# saveRDS(factor_cell_cor_data, file.path(DIR_RES, "..", "human", "objects", "hs_visium_preproc_B_IPF1-3_metadata.rds"))

factor_cell_cor_data <- readRDS(file.path(DIR_RES, "..", "human", "objects", "hs_visium_preproc_B_IPF1-3_metadata.rds"))

hs_c2l_colnames <- grep("c2l_", colnames(factor_cell_cor_data$IPF_1), value = T)
hs_c2l_selected <- c("c2l_B.Cells",
                     "c2l_Plasma.Cells",
                     "c2l_T.Cells",
                     "c2l_Proliferating.T.Cells",
                     "c2l_cDCs",
                     "c2l_pDCs",
                     "c2l_Monocytes",
                     "c2l_Macrophages",
                     "c2l_NK.Cells")


mm_c2l_colnames <- grep("c2l_", colnames(se_mmd21[[]]), value = T)
mm_c2l_selected <- c("c2l_B.lymphocytes",
                     "c2l_Plasma.cells",
                     "c2l_T.lymphocytes",
                     "c2l_T.cell.subset",
                     "c2l_Themis..T.lymphocytes",
                     grep(x = mm_c2l_colnames, pattern="DC", value = T),
                     grep(x = mm_c2l_colnames, pattern="mono", value = T),
                     grep(x = mm_c2l_colnames, pattern="AM|IM|macr", value = T),
                     "c2l_NK.cells")

#' Add metadata to hs se list
hs_se_list <- setNames(
  lapply(hs_donor_select, function(s){
    se <- hs_se_list[[s]]
    se <- AddMetaData(se, metadata = factor_cell_cor_data[[s]][colnames(se), hs_c2l_colnames])
    return(se)
  }), 
  nm = hs_donor_select)


#' Calculate pct present and avg cell density
for(c in hs_c2l_selected){
  dat <- subset(hs_se_list$IPF_1@meta.data, tls_factor == "high")
  pct_expr <- (sum(dat[[c]]>0.5) / nrow(dat)) * 100
  avg_expr <- mean(dat[[c]])
  message(paste(c, ":", pct_expr %>% round(digits = 2), ":", avg_expr %>% round(digits = 2)))
}

#' Human
tls_spot_c2l_stats_empty <- data.frame(cell_type = hs_c2l_selected,
                                       subject = rep("", length(hs_c2l_selected)),
                                       cell_avg_expr = rep(0, length(hs_c2l_selected)),
                                       cell_pct_expr = rep(0, length(hs_c2l_selected)))
tls_spot_c2l_stats_list_hs <- setNames(
  lapply(hs_donor_select, function(s){
    tls_spot_c2l_stats <- tls_spot_c2l_stats_empty
    tls_spot_c2l_stats$subject <- s
    
    dat <- subset(hs_se_list[[s]]@meta.data, tls_factor == "high")
    for(c in hs_c2l_selected){
      avg_expr <- mean(dat[[c]])
      pct_expr <- (sum(dat[[c]]>0.5) / nrow(dat)) * 100
      tls_spot_c2l_stats[tls_spot_c2l_stats$cell_type == c, c("cell_avg_expr", "cell_pct_expr")] <- c(avg_expr, pct_expr)
    }
    return(tls_spot_c2l_stats)
  }), 
  nm=hs_donor_select)

#' Mouse
mm_cell_types_rm <- c(grep("AM", mm_c2l_selected, value = T),
                      grep("Non.classical.monocytes", mm_c2l_selected, value = T),
                      grep("Fn1", mm_c2l_selected, value = T),
                      grep("M2", mm_c2l_selected, value = T),
                      grep("Themis", mm_c2l_selected, value = T),
                      grep("T.cell", mm_c2l_selected, value = T))
paste0(mm_cell_anno$cell_name_special[match(mm_cell_types_rm, mm_cell_anno$c2l_col_name)], collapse = ", ")
       
mm_c2l_selected_filtered <- mm_c2l_selected[!mm_c2l_selected %in% mm_cell_types_rm]
mm_c2l_selected_use <- mm_c2l_selected_filtered

tls_spot_c2l_stats_empty <- data.frame(cell_type = mm_c2l_selected_use,
                                       subject = rep("", length(mm_c2l_selected_use)),
                                       cell_avg_expr = rep(0, length(mm_c2l_selected_use)),
                                       cell_pct_expr = rep(0, length(mm_c2l_selected_use)))
tls_spot_c2l_stats_list_mm <- setNames(
  lapply(c("mm_d7", "mm_d21"), function(s){
    tls_spot_c2l_stats <- tls_spot_c2l_stats_empty
    tls_spot_c2l_stats$subject <- s
    
    dat <- subset(mm_se_list[[s]]@meta.data, tls_factor == "high")
    for(c in mm_c2l_selected_use){
      # message(c)
      avg_expr <- mean(dat[[c]])
      pct_expr <- (sum(dat[[c]]>0.5) / nrow(dat)) * 100
      tls_spot_c2l_stats[tls_spot_c2l_stats$cell_type == c, c("cell_avg_expr", "cell_pct_expr")] <- c(avg_expr, pct_expr)
    }
    return(tls_spot_c2l_stats)
  }), 
  nm=c("mm_d7", "mm_d21"))

#' Plot human
tls_spot_c2l_stats_hs <- bind_rows(tls_spot_c2l_stats_list_hs)
tls_spot_c2l_stats_hs <- tls_spot_c2l_stats_hs %>% arrange(cell_type)

tls_spot_c2l_stats_hs$cell_type <- hs_cell_anno$cell_name_special[match(gsub("c2l_", "", tls_spot_c2l_stats_hs$cell_type), hs_cell_anno$cell_name)]
hs_cell_levels <- hs_cell_anno$cell_name_special[match(gsub("c2l_", "", hs_c2l_selected), hs_cell_anno$cell_name)]

tls_spot_c2l_stats_hs$cell_type <- factor(tls_spot_c2l_stats_hs$cell_type, 
                                          levels = hs_cell_levels %>% rev())


p1 <- ggplot(tls_spot_c2l_stats_hs, aes(x=subject, y=cell_type, size=cell_pct_expr, color=cell_avg_expr)) +
  geom_point() +
  scale_color_gradientn(colours = col_scale_acton, limits=c(0,max(tls_spot_c2l_stats_hs$cell_avg_expr))) +
  labs(color = "Avg. cell density", size = "Detection rate (%)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "left",
        axis.title = element_blank(),
        axis.text = element_text(size=10, color="black"),
        axis.text.x = element_text(angle=90)) +
  scale_size(limits = c(0,100));p1


#' Plot mouse
tls_spot_c2l_stats_mm <- bind_rows(tls_spot_c2l_stats_list_mm)
tls_spot_c2l_stats_mm <- tls_spot_c2l_stats_mm %>% arrange(cell_type)

tls_spot_c2l_stats_mm$cell_type <- mm_cell_anno$cell_name_special[match(tls_spot_c2l_stats_mm$cell_type, mm_cell_anno$c2l_col_name)]
mm_cell_levels <- mm_cell_anno$cell_name_special[match(mm_c2l_selected_use, mm_cell_anno$c2l_col_name)]

tls_spot_c2l_stats_mm$cell_type <- factor(tls_spot_c2l_stats_mm$cell_type, 
                                          levels = mm_cell_levels %>% rev())

tls_spot_c2l_stats_mm$subject <- factor(tls_spot_c2l_stats_mm$subject, levels = c("mm_d7", "mm_d21"))


p2 <- ggplot(tls_spot_c2l_stats_mm, aes(x=subject, y=cell_type, size=cell_pct_expr, color=cell_avg_expr)) +
  geom_point() +
  # scale_color_gradientn(colours = col_scale_rocket) +
  scale_color_gradientn(colours = col_scale_acton, limits=c(0,max(tls_spot_c2l_stats_mm$cell_avg_expr))) +
  scale_y_discrete(position = "right") +
  labs(color = "Avg. cell density", size = "Detection rate (%)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        axis.title = element_blank(),
        axis.text = element_text(size=10, color="black"),
        axis.text.x = element_text(angle=90)) +
  scale_size(limits = c(0,100));p2

# pdf(file.path(DIR_FIG_OUT, "hs_mm_visium_tls_factor_high_c2l_enrichment.pdf"), 
#     width = 5, height = 5)
pdf(file.path(DIR_FIG_OUT, "hs_mm_visium_tls_factor_high_c2l_enrichment_filt.pdf"),  # !Manuscript Main figure 6c
    width = 7.5, height = 3)
(p1|p2) # & theme(legend.position = "bottom")
dev.off()


# Export data
tls_spot_c2l_stats_export <- bind_rows(tls_spot_c2l_stats_hs, tls_spot_c2l_stats_mm)
write.csv(x = tls_spot_c2l_stats_export, file = file.path(DIR_OBJ_OUT, "hs_mm_visium_tls_factor_high_c2l_enrichment.csv"), row.names = F)

# Export metadata of TLS selection
mdat_tls_selection_mm <- bind_rows(mm_se_list$mm_d7@meta.data, mm_se_list$mm_d21@meta.data)
mdat_tls_selection_hs <- bind_rows(hs_se_list$IPF_1@meta.data, hs_se_list$IPF_2@meta.data, hs_se_list$IPF_3@meta.data)

mdat_tls_selection_mm$barcode <- rownames(mdat_tls_selection_mm)
mdat_tls_selection_hs$barcode <- rownames(mdat_tls_selection_hs)

write.csv(x = mdat_tls_selection_mm, file = file.path(DIR_OBJ_OUT, "hs_mm_visium_tls_factor_se_metadata_mm.csv"), row.names = T)
write.csv(x = mdat_tls_selection_hs, file = file.path(DIR_OBJ_OUT, "hs_mm_visium_tls_factor_se_metadata_hs.csv"), row.names = T)



