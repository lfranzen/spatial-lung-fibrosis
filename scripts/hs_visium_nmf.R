#' [hs_visium_nmf.R]
#'
#' Run, analyse, and plot NMF results 
#'
#' L. Franzén [lovisa.franzen@scilifelab.se]

#### Set up ####
##### Define params. ####
set.seed(1)
SPECIES <- "human"
DIR_ROOT <- getwd()
DIR_DATA <- file.path(DIR_ROOT, "data", SPECIES, "visium")
DIR_RES <- file.path(DIR_ROOT, "results", SPECIES)
DIR_FIG <- file.path(DIR_RES, "figures")
fig_res <- 300
n_factors <- 30
DIR_FIG_NMF <- file.path(DIR_RES, "figures", paste0("nmf_", n_factors, "_factors"))
dir.create(path = DIR_FIG_NMF)

##### Load libs ####
library(STutility)
library(RColorBrewer)
library(patchwork)
library(dplyr)
library(magrittr)
library(writexl)
library(gprofiler2)
library(pheatmap)

##### Other ####
# source(file.path(DIR_ROOT, "scripts", "colors.R"))
source(file.path(DIR_ROOT, "scripts", "custom_functions.R"))
source(file.path(DIR_ROOT, "scripts", "custom_colors.R"))
theme_custom <- theme(axis.title.x = element_blank())

##### Read objects ####
fname <- paste0("hs_visium_stutility_obj.rds")
se.subset <- readRDS(file = file.path(DIR_RES, "objects", fname))

#### NMF ####
##### Run NMF ####
#' ...if it hasn't been run before
n_factors <- 30
se.subset <- RunNMF(se.subset, nfactors = n_factors)


##### NMF Results ####

#' Save top genes for each factor
factor_gene_loadings <- lapply(1:n_factors, function(factor_x){
  message(paste("Factor", factor_x))
  feat_loads <- as.data.frame(se.subset@reductions$NMF@feature.loadings[,paste0("factor_",factor_x)])
  colnames(feat_loads) <- "gene_loading"
  feat_loads$gene_loading_scaled <- Scale01(feat_loads$gene_loading)
  feat_loads$gene <- rownames(feat_loads)
  feat_loads$factor <- factor_x
  top_n_genes <- 100
  feat_loads <- feat_loads %>%
    dplyr::slice_max(order_by = gene_loading, n = top_n_genes) %>%
    mutate(rank = dense_rank(desc(gene_loading)))
})

write_xlsx(
  factor_gene_loadings,
  path = file.path(DIR_RES, "objects", "A_NMF30", paste0("hs_visium_A_nmf_1-", n_factors, "_top100_gene_loadings.xlsx")),
  col_names = TRUE,
  format_headers = TRUE,
  use_zip64 = FALSE
)


##### gProfiler ####
library(gprofiler2)

org <- "hsapiens"
fname <- paste0("hs_visium_A_nmf_1-", n_factors,"_gProfiler_top25genes")
n_genes <- 25
gene_query_list <- list()
gostres_list <- list()

for(i in 1:n_factors){
  gene_query_list[[i]] <- subset(factor_gene_loadings[[i]], rank <= n_genes)$gene
  }
names(gene_query_list) <- paste0("factor_", 1:n_factors)

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
  path = file.path(DIR_RES, "objects", "A_NMF30", paste0(fname, ".xlsx")),
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
selected_samples <- c(1, 2, 3, 4, 16, 21, 22, 19)
plot_list <- lapply(1:n_factors, function(i) {
  print(i)
  FactorPlot(seurat.object = se.subset, 
             indices = selected_samples,
             n.columns = 4,
             factor = i, 
             col.scale = rev(RColorBrewer::brewer.pal(10, "Spectral")))
})
png(file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors, "_a.png")), width = 24*fig_res, height = 26*fig_res, res = fig_res)
cowplot::plot_grid(plotlist = plot_list[1:10], ncol = 2)
dev.off()
png(file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors, "_b.png")), width = 24*fig_res, height = 26*fig_res, res = fig_res)
cowplot::plot_grid(plotlist = plot_list[11:20], ncol = 2)
dev.off()
if(n_factors==30){
  png(file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors, "_c.png")), width = 24*fig_res, height = 26*fig_res, res = fig_res)
  print(cowplot::plot_grid(plotlist = plot_list[21:30], ncol = 2))
  dev.off()
}


###### Gene loadings ####
factor_gene_loading_plot <- lapply(1:n_factors, function(i){
  n_genes <- 10
  dat_plot <- factor_gene_loadings[[i]]
  p <- ggplot(dat_plot, aes(x=rank, y=gene_loading)) +
    geom_line(size=.5) +
    geom_point(data = dat_plot[1:n_genes,], mapping = aes(x=rank, y=gene_loading, color=reorder(gene, rank))) +
    scale_color_manual(values = RColorBrewer::brewer.pal(n_genes, "Spectral")) +
    labs(title = paste("Factor", i), x="gene rank", y="loading", color="") +
    theme_bw() +
    theme(legend.position = "right", legend.justification = "left", 
          legend.text = element_text(size=6), 
          legend.box.margin=margin(t = -.1, l = -.25, unit='cm'),  # pos right
          # legend.box.margin=margin(t = -.1, l = -1, unit='cm'),  # pos bottom
          legend.key.size = unit(3, "mm"),
          # axis.text.x = element_blank(),
          plot.title = element_text(hjust=0.5, face="bold"),
          panel.grid = element_blank())
})
png(file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors, "_gene_loadings.png")), width = 17*fig_res, height = 12*fig_res, res = fig_res)
cowplot::plot_grid(plotlist = factor_gene_loading_plot[1:n_factors], ncol = 6)
dev.off()



###### Rank plot ####
# Grade
plot_list <- lapply(1:n_factors, function(i) {
  p <- FactorRankSpotPlot(seurat.object = se.subset, factor = i, split.by = "sample_name", group.by = "fibrotic_extent_score_by_pathologist_0.3", color.group = cols_grade)
  })
png(file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors, "_ranks.png")), width = 20*fig_res, height = 18*fig_res, res = fig_res)
cowplot::plot_grid(plotlist = plot_list[1:n_factors], ncol = 4)
dev.off()

# Donor
plot_list <- lapply(1:n_factors, function(i) { # n_factors
  p <- FactorRankSpotPlot(seurat.object = se.subset, factor = i, split.by = "sample_name", group.by = "subject_alias", color.group = cols_donor)
})
# plot_list[[14]] & theme(panel.background = element_rect(fill="grey10"))
FactorRankSpotPlot(seurat.object = se.subset, factor = 14, split.by = "sample_name", group.by = "subject_alias", color.group = cols_donor, include.zoom = F, top.n.spots = 300)

png(file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors, "_ranks_donor.png")), width = 20*fig_res, height = 18*fig_res, res = fig_res)
cowplot::plot_grid(plotlist = plot_list[1:n_factors], ncol = 4)
dev.off()


# Annotation
se_annotated <- SubsetSTData(se.subset, expression = (annotation != "NA"))
plot_list <- lapply(1:n_factors, function(i) {
  p <- FactorRankSpotPlot(seurat.object = se_annotated, factor = i, split.by = "annotation", group.by = "annotation", color.group = cols_annotation)
})
png(file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors, "_ranks_histopath.png")), width = 20*fig_res, height = 18*fig_res, res = fig_res)
cowplot::plot_grid(plotlist = plot_list[1:n_factors], ncol = 4)
dev.off()

# rm(se_annotated)


###### Factor plot grid ####
plot_list <- lapply(1:n_factors, function(factor_plot) {
  p1 <- ST.DimPlot(se.subset, 
                   # indices = c(1:3, 13:15, 4:6, 16:18, 7:9, 19:21),
                   indices = c(1:12, 20:23),
                   dims = factor_plot, 
                   min.cutoff = 0.2,
                   label.by = "sample_name",
                   center.zero = F, 
                   reduction = "NMF", 
                   ncol = 8, 
                   pt.size = 0.8, 
                   pt.border = F,
                   cols = rev(RColorBrewer::brewer.pal(10, "Spectral")),
                   # cols = viridis::magma(10, direction = -1),
                   show.sb = F)  &
    labs(fill = "factor\nactivity", title = paste0("Factor ", factor_plot)) & 
    theme(legend.position = "right", plot.title = element_text(face = "bold", size = 16))
  p2 <- FactorGeneLoadingPlotHorizontal(se.subset, 
                                        factor = factor_plot, 
                                        topn = 30) + 
    labs(title = "Genes contributing to factor") & 
    theme(axis.text.x = element_text(size=10), panel.grid = element_blank(),
          plot.margin = margin(10, 10, 0, 10, unit = "pt"))
  
  p3 <- FactorRankSpotPlot(se.subset, factor = factor_plot, split.by = "sample_name", 
                           group.by = "fibrotic_extent_score_by_pathologist_0.3", color.group = cols_grade, 
                           plot.title = "Spots ranked by factor activity")
  
  p_bottom <- cowplot::plot_grid(p2, p3, nrow = 1)
  p <- cowplot::plot_grid(p1, p_bottom, ncol = 1, rel_heights = c(3, 2))
})

out_path <- file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors, "_grid"))
dir.create(path = out_path)
for(i in 1:length(plot_list)){
  png(file = file.path(out_path, paste0("hs_visium_A_nmf_grid_f", i, ".png")),
      width = 14*fig_res, height = 7*fig_res, res = fig_res)
  print(plot_list[[i]])
  dev.off()
  }

# pdf(file = file.path(DIR_FIG, paste0("hs_visium_A_nmf_1-", n_factors, "_grid.pdf")), width = 16, height = 7, pointsize = 4, fg = "white")
# print(plot_list[1])
# dev.off()

###### Top rank box plot plot ####
factor_select <- "factor_14"
rank_dat <- cbind(se.subset@reductions$NMF@cell.embeddings[,factor_select], se.subset@meta.data)
colnames(rank_dat)[1] <- "nmf_factor"
rank_dat$group_data_by <- rank_dat[ ,"sample_name"]
rank_dat$subject_group <- ifelse(grepl("HC", rank_dat$subject_alias), "HC", rank_dat$subject_alias)

rank_dat_plot <- rank_dat %>%
  dplyr::group_by(group_data_by, subject_group) %>%
  dplyr::slice_max(order_by = nmf_factor, n = 10) %>%
  mutate(rank = dense_rank(desc(nmf_factor)))

p <- ggplot(rank_dat_plot, aes(x=subject_group, y=nmf_factor, fill=subject_group)) +
  geom_boxplot(color="black", outlier.colour = "black", outlier.size = 0.8) +
  # geom_point(shape = 21) +
  labs(x="", y="factor activity", title="Top 10 F14 enriched spots per sample") +
  scale_fill_manual(values = cols_donor[c(1,5:8)]) +
  theme_linedraw() +
  theme(panel.grid = element_blank(), 
        plot.title = element_text(size=12),
        legend.position = "none");p

p2 <- FactorRankSpotPlot(seurat.object = se.subset, factor = 14, split.by = "sample_name", group.by = "subject_alias", color.group = cols_donor, 
                         include.zoom = T, top.n.spots = 300, plot.title = "Top 300 ranked F14 enriched spots per sample") & 
  theme(plot.title = element_text(size=10), legend.position = "none", text = element_text(color="black"));p2

png(file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors, "_ranks_box_donor.png")), width = 3.8*fig_res, height = 3.8*fig_res, res = fig_res)
p2/p
dev.off()



###### UMAP & clusters & NMF ####
n_clusters <- length(unique(se.subset$seurat_clusters))
col_scale <- viridis::rocket(n_clusters, direction = -1)
plot_list <- lapply(1:n_factors, function(factor_plot) {
  factor_name <- paste0("factor_", factor_plot)
  p1 <- VlnPlot(se.subset, features = factor_name, group.by = "seurat_clusters", cols = rev(col_scale), pt.size = 0) & 
    labs(title = paste("Factor", factor_plot), x = "Cluster") &
    coord_flip() & 
    NoLegend() & 
    theme(axis.text.x = element_text(angle=0, hjust = 0.5))
  p2 <- FeaturePlot(se.subset, features = factor_name, reduction = "umap.harmony", cols = col_scale, min.cutoff = .2, raster = T) &
    labs(title = "", x = "UMAP 1", y="UMAP 2") & 
    NoLegend()
  p <- p1 + p2 + plot_layout(ncol = 2, widths = c(1,2))
})

png(file = file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors, "_umap_clusters.png")), width = 20*fig_res, height = 18*fig_res, res = fig_res)
cowplot::plot_grid(plotlist = plot_list[1:n_factors], ncol = 4, nrow = 5)
dev.off()


###### NMF on H&E image of selected samples ####
sample_name_select <- c("IPF_1.TD012.B2.1", "IPF_2.TD022.B3.2", "IPF_3.TD031.B1.2", "IPF_3.TD032A.B2.2")
se.he <- SubsetSTData(se.subset, expression = sample_name %in% sample_name_select)
se.he <- LoadImages(se.he, xdim = 2e3)

factor_select <- c(paste0("factor_", c(4, 5, 6, 8, 9, 14)), "CXCL13")
fname <- paste0("hs_visium_A_nmf_1-", n_factors,"_HE_selected_factors")
pdf(file = file.path(DIR_FIG, paste0(fname, ".pdf")), width = 10, height = 10, useDingbats = F)
for(f in factor_select){
  message(f)
  p <- FeatureOverlay(se.he, features = f,
                      sampleids = 1:4, label.by = "sample_name",
                      ncols = 2,
                      add.alpha = T, palette = "Spectral", pt.size = 1, show.sb = F)
  print(p)
}
dev.off()


###### Fibrotic Factor HE plots ####
fib_genes <- sort(c(grep("^COL[0-9]", rownames(se.subset@reductions$NMF@feature.loadings), value = T), "FN1", "FBN1", "THY1", "TNC", "VIM"))
nmf_gene_loads <- se.subset@reductions$NMF@feature.loadings
nmf_gene_loads <- nmf_gene_loads[rownames(nmf_gene_loads)%in%fib_genes,]
fib_factor_sum <- sort(colSums(nmf_gene_loads), decreasing = T)
factors_plot <- names(head(fib_factor_sum, 3)) # not really reflecting ECM factors? 
# factors_plot <- paste0("factor_", c(4,7,14,23,30))
factors_plot <- paste0("factor_", c(4,7,14))

subject_ipf_names <- grep("IPF", unique(se.subset$subject_alias), value = T)
sample_name_select <- list(HC = c("HC_1.TH010.B0.1", "HC_3.TH030.B0.2", "HC_4.TH040.B0.2"),
                           IPF_1 = c("IPF_1.TD011.B1.1", "IPF_1.TD012.B2.1", "IPF_1.TD013.B3.1"),
                           IPF_2 = c("IPF_2.TD021A.B1.1", "IPF_2.TD021B.B2.1", "IPF_2.TD022.B3.1"),
                           IPF_3 = c("IPF_3.TD031.B1.2", "IPF_3.TD032A.B2.2", "IPF_3.TD032B.B3.2"),
                           IPF_4 = c("IPF_4.TD041A.B1.1", "IPF_4.TD041B.B2.1", "IPF_4.TD042.B3.1")
                           )


pdf(file = file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors,"_HE_fib_factors", ".pdf")), width = 18, height = 9, useDingbats = F)
for(subject in c("HC",subject_ipf_names)){
  message(subject)
  se_he_ipf <- SubsetSTData(se.subset, expression = sample_name %in% sample_name_select[[subject]])
  se_he_ipf <- LoadImages(se_he_ipf, xdim = 1e3)
  p1 <- FeatureOverlay(se_he_ipf, 
                       sampleids = 1:3,
                       features = factors_plot,
                       label.by = "sample_name",
                       ncols = 3, 
                       blend = T,
                       add.alpha = T,
                       pt.size = 2, 
                       show.sb = F)
  p1 <- ggrastr::rasterize(p1, layers = "Point", dpi = 300)
  p_list <- lapply(factors_plot, function(fac){
    d_plot <- subset(factor_gene_loadings[[as.numeric(gsub("factor_", "", fac))]], rank <= 40)
    ggplot(d_plot, aes(x=reorder(gene, rank), y=gene_loading_scaled)) +
      geom_point() +
      labs(x="", y="", title=fac) +
      theme_classic() +
      theme(axis.text.x = element_text(angle=45, hjust=1), plot.title = element_text(hjust=.5),
            axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
            panel.grid = element_blank(), panel.grid.major.x = element_line(size = .3, color="black"))
  })
  p2 <- patchwork::wrap_plots(p_list, nrow = 1)
  p <- p1/p2+patchwork::plot_layout(heights = c(5,1))
  message("Saving plot...")
  print(p)
  }
dev.off()
rm(se_he_ipf)



factor_plot <- "factor_14"
se_he <- SubsetSTData(se.subset, expression = sample_name %in% unlist(sample_name_select))
se_he <- LoadImages(se_he, xdim = 2e3)
factor_max <- round(as.numeric(summary(se_he@reductions$NMF@cell.embeddings[,factor_plot])[6])-.5)
factor_min <- round(as.numeric(summary(se_he@reductions$NMF@cell.embeddings[,factor_plot])[2])+.5)
pdf(file = file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors,"_HE_", factor_plot, ".pdf")), width = 10, height = 10, useDingbats = F)
for(i in 1:length(unlist(sample_name_select))){  #:length(unlist(sample_name_select))
  message(i)
  p <- FeatureOverlay(se_he, 
                      sampleids = i,
                      features = factor_plot,
                      max.cutoff = factor_max, 
                      min.cutoff = factor_min,
                      label.by = "sample_name",
                       add.alpha = T,
                       pt.size = 2, 
                       show.sb = F);p
  p <- ggrastr::rasterize(p, layers = "Point", dpi = 300)
  print(p)
}
dev.off()



##### Fibrotic factor analysis ####
#' Further analysis of factors enriched in ECM and fibrosis-related genes
#' 
#' 
####### Heatmap of COLs and KRTs #####
fib_genes <- sort(c(grep("^COL[0-9]", rownames(se.subset@reductions$NMF@feature.loadings), value = T), 
                    "FN1", "FBN1", "THY1", "TNC", "VIM", "ACTA2", "TAGLN", "LUM", "DCN", "PDGRFA",
                    grep("^KRT[0-9]", rownames(se.subset@reductions$NMF@feature.loadings), value = T)))
nmf_gene_loads <- se.subset@reductions$NMF@feature.loadings
nmf_gene_loads_scaled <- apply(nmf_gene_loads, 2, Scale01)
nmf_gene_loads_scaled <- nmf_gene_loads_scaled[rownames(nmf_gene_loads_scaled)%in%fib_genes,]
nmf_gene_loads_scaled <- as.data.frame(t(nmf_gene_loads_scaled))
nmf_gene_loads_scaled$sum <- as.numeric(rowSums(nmf_gene_loads_scaled))
nmf_gene_loads_scaled_filtered <- nmf_gene_loads_scaled[nmf_gene_loads_scaled$sum>2, fib_genes[fib_genes%in%colnames(nmf_gene_loads_scaled)]]
nmf_gene_loads_scaled_filtered <- nmf_gene_loads_scaled_filtered[, c(colnames(nmf_gene_loads_scaled_filtered)[colSums(nmf_gene_loads_scaled_filtered)>0.3],"KRT5")]

p <- pheatmap(nmf_gene_loads_scaled_filtered, display_numbers = F, color = RColorBrewer::brewer.pal(8,"BuPu"));p
png(file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors, "_selected_gene_heatmap.png")), width = 10*fig_res, height = 4*fig_res, res = fig_res);p;dev.off()

# pz <- pheatmap(ZScoreMartrix(nmf_gene_loads_scaled_filtered, by_row = F), display_numbers = F, breaks = c(-3,-2,-1,0,1,2,3,4), color = rev(RColorBrewer::brewer.pal(7,"RdYlBu")))
# png(file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors, "_selected_gene_heatmap_zscore.png")), width = 10*fig_res, height = 4*fig_res, res = fig_res);pz;dev.off()


####### Label FibFactor-high spot  #####
#' Label spots top 1% most factor enriched spots of the fbrotic factors
se_nmf_emb <- se.subset@reductions$NMF@cell.embeddings
factor_xy <- names(sort(rowSums(nmf_gene_loads_scaled), decreasing = T))[1:4]
# factor_xy <- c("factor_14", "factor_4", "factor_16", "factor_10")

factor_cutoff <- lapply(factor_xy, function(x){
  as.numeric(quantile(se_nmf_emb[,x], c(.99)))
})
names(factor_cutoff) <- factor_xy
f_metadata_add <- se.subset@meta.data[,1:2]
for(f in factor_xy){
  f_metadata_add[[f]] <- paste0(f, "_low")
  spots_f_cutoff <- intersect(colnames(se.subset), rownames(se_nmf_emb[se_nmf_emb[,f]>factor_cutoff[[f]],]))
  f_metadata_add[spots_f_cutoff,f] <- paste0(f, "_high")
}
se.subset <- AddMetaData(se.subset, metadata = f_metadata_add[,factor_xy], col.name = paste0(factor_xy, "_cutoff"))


pdf(file = file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors, "_selected_cutoff_spatial.pdf")), width = 9, height = 10, useDingbats = F)
for(ft in paste0(factor_xy, "_cutoff")[1]){
  message(ft)
  p <- ST.FeaturePlot(se.subset, features = ft, pt.size = .6, ncol = 5, label.by = "sample_name", cols = c("grey80", "grey10")) & 
    theme(legend.position = "bottom") & guides(fill = guide_legend(override.aes = list(size=3)))
  p <- ggrastr::rasterize(p, layers = "Point", dpi = 300)
  print(p)
}
dev.off()


#' Plot proportions per sample
dat_summary <- se.subset@meta.data %>%
  group_by(subject_alias, factor_14_cutoff) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
dat_merge <- se.subset@meta.data %>% 
  group_by(subject_alias, condition) %>%
  summarise(n_spots = n())
dat_summary <- merge(dat_summary, dat_merge, by = "subject_alias")

ggplot(subset(dat_summary, factor_14_cutoff=="factor_14_high"), aes(x = subject_alias, y = freq*100, fill = as.factor(factor_14_cutoff))) + 
  geom_bar(stat = 'identity', colour=NA, position = "stack") +
  # facet_wrap(~condition, scales = "free_x") +
  facet_grid(~condition, scales = "free_x", space='free') +
  scale_fill_manual(values = c("grey10", "grey80")) +
  labs(x="", y="% of spots", fill="") +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        panel.grid = element_blank(),
        legend.position = "none", 
        legend.text = element_text(size=6))


##### F14 analysis ####
#' Further analysis of F14-enriched spots
summary(as.factor(se.subset$factor_14_cutoff))

###### Subcluster F14-high spots #####
spots_f14 <- rownames(se.subset@meta.data[se.subset$factor_14_cutoff %in% "factor_14_high",])
se.f14high <- SubsetSTData(se.subset, spots = spots_f14)

se.f14high <- RunPCA(se.f14high)
# ElbowPlot(se.f14high)
# DimHeatmap(se.f14high, reduction = "pca", dims = 1:12)
# VlnPlot(se.f14high, features = paste0("PC_",1:10), group.by = "subject_alias", ncol=5)

#' PCA only:
dims_select <- 1:8 #1:10
se.f14high <- RunUMAP(se.f14high, reduction = "pca", dims = dims_select)
se.f14high <- FindNeighbors(se.f14high, reduction = "pca", dims = dims_select)
se.f14high <- FindClusters(se.f14high, resolution = 0.4)

#' Save clusters
se.f14high$f14_subclusters <- se.f14high$seurat_clusters
n_clusters <- length(unique(se.f14high$f14_subclusters))
cols_clusters <- setNames(object = ggsci::pal_locuszoom(palette = "default")(n_clusters), #c(brewer.pal(n_clusters, "Set1")), 
                          nm = levels(se.f14high$f14_subclusters))

# Plot umap
DimPlot(se.f14high, reduction = "umap", group.by = "subject_alias")
DimPlot(se.f14high, reduction = "umap", group.by = "f14_subclusters", split.by = "subject_alias")

p1 <- DimPlot(se.f14high, reduction = "umap", group.by = "seurat_clusters", label = T, label.box = T, repel = T, cols = cols_clusters) & NoLegend() & theme_umap
p2 <- FeaturePlot(se.f14high, features = c("KRT5", "BPIFB1", "MT1E", "KRT17", "PRSS2", "MMP7"), ncol=3, cols = col_scale_mako, pt.size=0.25) & theme_void() & NoLegend() & theme(plot.title = element_text(face = "bold", size = 10, colour = "black"))
p <- (p1 & p2) & plot_annotation(title = "Subclustering of F14-high spots");p

png(file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors, "_f14_subclusters_umap.png")), width = 7*fig_res, height = 4*fig_res, res = fig_res);p;dev.off()

# Marker genes
se.f14high <- SetIdent(se.f14high, value = "f14_subclusters")
makers_f14_subcluster <- FindAllMarkers(se.f14high)

#'export
fname <- paste0("hs_visium_A_nmf_1-", n_factors, "_f14_subclusters_markers.csv")
write.csv(x = makers_f14_subcluster, file = file.path(DIR_RES, "objects", "A_NMF30", fname), row.names = F)
makers_f14_subcluster <- read.csv(file.path(DIR_RES, "objects", "A_NMF30", fname))

saveRDS(object = se.f14high, file = file.path(DIR_RES, "objects", "hs_visium_preproc_A_se_obj_f14high_subset.rds"))
se.f14high <- readRDS(file.path(DIR_RES, "objects", "hs_visium_preproc_A_se_obj_f14high_subset.rds"))


# Plot marker heatmap 
genes_plot <- subset(makers_f14_subcluster, p_val_adj<0.05)
genes_plot$updown <- "up"
genes_plot[genes_plot$avg_log2FC<0,"updown"] <- "down"
genes_plot <- genes_plot %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group=TRUE) %>%
  top_n(20, (avg_log2FC))

mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(256)
p <- DoHeatmap(se.f14high, features = c("KRT17", "KRT5", genes_plot$gene), group.by = "f14_subclusters", group.colors =  cols_clusters, slot = "scale.data") +
  scale_fill_gradientn(colours = rev(mapal))
png(file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors, "_f14_subclusters_marker_heatmap.png")), width = 12*fig_res, height = 12*fig_res, res = fig_res);p;dev.off()

mapal <- colorRampPalette(RColorBrewer::brewer.pal(8,"OrRd"))(256)
p <- DoHeatmap(se.f14high, features = c("KRT17", "KRT5", genes_plot$gene), 
               group.by = "f14_subclusters", 
               group.colors =  cols_clusters, 
               slot = "data") +
  scale_fill_gradientn(colours = mapal)
png(file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors, "_f14_subclusters_marker_heatmap_data.png")), width = 12*fig_res, height = 16*fig_res, res = fig_res);p;dev.off()

# DoHeatmap(se.f14high, features = c("KRT17", "KRT5", genes_plot$gene), group.by = "f14_subclusters", group.colors =  cols_clusters, slot = "counts", disp.max = 6) +
#   scale_fill_gradientn(colours = colorRampPalette(col_scale_rocket)(256))

# Plot marker dotplot (Ext Data Fig 3a)
genes_plot_dp <- genes_plot %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group=TRUE) %>%
  top_n(20, (avg_log2FC)) %>% 
  pull(gene) %>% 
  unique()

p <- DotPlot(se.f14high, features = rev(genes_plot_dp), 
             group.by = "f14_subclusters", scale = T) & 
  theme_dotplot &
  scale_color_gradientn(colours = col_scale_div_expr, limits=c(-2,2)) &
  scale_y_discrete(position = "right") &
  coord_flip();p
  # theme(
  #   # axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
  #   axis.title = element_blank(), 
  #   legend.position = "right") + 
  # coord_flip()
  # RotatedAxis()

png(file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors, "_f14_subclusters_marker_dotplot.png")), width = 4*fig_res, height = 8*fig_res, res = fig_res);p;dev.off()

pdf(file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors, "_f14_subclusters_marker_dotplot.pdf")), width = 4, height = 12);p;dev.off()


# Plot spatial HE
he_samples_select <- c("IPF_1.TD013.B3.1", "IPF_3.TD031.B1.2", "IPF_3.TD032B.B3.2", "IPF_4.TD042.B3.1") # "IPF_1.TD012.B2.1", "IPF_4.TD041B.B2.1", "IPF_4.TD042.B3.1"
se.f14high_he <- SubsetSTData(se.f14high, expression = sample_name %in% he_samples_select)
se.f14high_he <- LoadImages(se.f14high_he, xdim = 2000)

p <- FeatureOverlay(se.f14high_he, features = "f14_subclusters",
                    sampleids = 1:4, cols = cols_clusters, pt.alpha = 0.8, 
                    pt.size = 1, ncol=2, label.by = "sample_name") & 
  theme(legend.position = "bottom")
p <- ggrastr::rasterize(p, layers = "Point", dpi = fig_res)
fname <- paste0("hs_visium_A_nmf_1-", n_factors, "_f14_subclusters_spatial_HE")
pdf(file = file.path(DIR_FIG_NMF, paste0(fname, ".pdf")), width = 10, height = 12, useDingbats = F);p;dev.off()


rm(se.f14high_he)


# Add metadata to se.subset object
se.f14high$f14_subclusters <- se.f14high$seurat_clusters
se.subset <- AddMetaData(se.subset, se.f14high$f14_subclusters, "f14_subclusters")

se.subset$f14_subclusters <- as.character(se.subset$f14_subclusters)
se.subset$f14_subclusters[is.na(se.subset$f14_subclusters)] <- "other"
se.subset$f14_subclusters <- as.factor(se.subset$f14_subclusters)
summary(as.factor(se.subset$f14_subclusters))

p <- DotPlot(se.subset, features = rev(genes_plot_dp$gene), group.by = "f14_subclusters", cols="RdYlBu", scale = T) + 
  theme_dotplot +
  scale_y_discrete(position = "right") +
  coord_flip();p
png(file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors, "_f14_subclusters_marker_dotplot_other.png")), width = 4*fig_res, height = 8*fig_res, res = fig_res);p;dev.off()


# F14-subcluster markers vs all rest
se.subset_ipf <- SubsetSTData(se.subset, condition == "IPF")
se.subset_ipf <- SetIdent(se.subset_ipf, value = "f14_subclusters")
# clusters_test <- levels(se.f14high$seurat_clusters)

subc0_marker_genes <- FindMarkers(se.subset_ipf,
                                 ident.1 = "0", 
                                 min.pct = 0.25, 
                                 min.diff.pct = 0.1)
subc0_marker_genes$gene <- rownames(subc0_marker_genes)
subc0_marker_genes$cluster <- "f14_high_C0"
subc0_marker_genes <- subc0_marker_genes %>% 
  filter(p_val_adj<0.05) %>% 
  arrange(desc(abs(avg_log2FC)))

#'export
fname <- paste0("hs_visium_A_nmf_1-", n_factors, "_f14_subcluster_C0_markers_vs_all_IPF_only.csv")  # !Important
write.csv(x = subc0_marker_genes, file = file.path(DIR_RES, "objects", "A_NMF30", fname), row.names = F)



# Plot spatial
cols_clusters_all <- setNames(object = c(cols_clusters, "grey90"), nm = levels(as.factor(se.subset$f14_subclusters)))
p <- ST.FeaturePlot(se.subset, features = "f14_subclusters", ncol=5, cols = cols_clusters_all, label.by = "sample_name", pt.size = 0.75, show.sb = F) + 
  guides(fill = guide_legend(override.aes = list(size=2)))
# p <- ggrastr::rasterize(p, layers = "Point", dpi = fig_res)
png(file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors, "_f14_subclusters_spatial.png")), width = 12*fig_res, height = 12*fig_res, res = fig_res);p;dev.off()


# Plot proportions
d <- se.subset@meta.data %>%
  group_by(sample_name, subject_alias, f14_subclusters) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

p1 <- ggplot(subset(d, f14_subclusters!="other"), aes(x=sample_name, y=freq, fill=f14_subclusters)) +
  geom_bar(stat = 'identity', colour=NA, position = "stack") +
  facet_grid(~subject_alias, scales = "free_x", space='free') +
  scale_fill_manual(values = cols_clusters) +
  labs(x="", y="Frequency", fill="F14-high\nsubclusters", title="F14 subcluster spot proportions") +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        panel.grid = element_blank(),
        legend.position = "right", 
        text = element_text(size=10),
        plot.title = element_text(hjust=0.5));p1
png(file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors, "_f14_subclusters_proportions.png")), width = 7*fig_res, height = 3.5*fig_res, res = fig_res);p1;dev.off()

pdf(file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors, "_f14_subclusters_proportions.pdf")), width = 7, height = 3.5);p1;dev.off()


# clean up
rm(se.f14high)


###### F14-high C0 nbs spots #####
# Neighbouring spots
se.subset <- LoadImages(se.subset, xdim=50)
se.subset$f14_subclusters <- as.character(se.subset$f14_subclusters)
se.subset <- SetIdent(se.subset, value ="f14_subclusters")
se.subset <- IdentifyNNeighbors(se.subset, feature.column.name = "f14_subclusters", center.feature = "0", n.neighbors = 3)
unique(se.subset$d_c0)
cols_dist <- c("1"=col_scale_spec[2], "2"=col_scale_spec[4])

se_nbs <- SubsetSTData(se.subset, d_c0 %in% c(1,2))

se_nbs <- RunPCA(se_nbs)
# ElbowPlot(se_nbs)
# DimHeatmap(se_nbs, reduction = "pca", dims = 1:16)
# VlnPlot(se_nbs, features = paste0("PC_", 1:15), group.by = "subject_alias", ncol=5, pt.size = 0)

dims_select <- 1:9
se_nbs <- RunUMAP(se_nbs, reduction = "pca", dims = dims_select)
se_nbs <- FindNeighbors(se_nbs, reduction = "pca", dims = dims_select)
se_nbs <- FindClusters(se_nbs, resolution = 0.2)
n_clusters_nbs <- length(unique(se_nbs$seurat_clusters))
cluster_cols_nbs <- setNames(object = RColorBrewer::brewer.pal(n_clusters_nbs+2, "Spectral")[-c(1,4)], 
                             nm = sort(unique(se_nbs$seurat_clusters)))

p1 <- DimPlot(se_nbs, reduction = "umap", group.by = "seurat_clusters", pt.size=0.5,  cols = cluster_cols_nbs, label = T, label.box = T, repel = T) + 
  labs(title="F14-C0 nbs clusters") +
  theme_umap
p2 <- DimPlot(se_nbs, reduction = "umap", group.by = "d_c0", pt.size=0.5, cols = cols_dist) + 
  labs(title="Distance to F14-C0") +
  theme_umap
(p1/p2)

# FeaturePlot(se_nbs, reduction = "umap", features = "PC_15", cols = col_scale_spec)

# Marker genes
se_nbs <- SetIdent(se_nbs, value = "seurat_clusters")
makers_f14_subcluster_nbs <- FindAllMarkers(se_nbs)

#'export
fname <- paste0("hs_visium_A_nmf_1-", n_factors, "_f14_subclusters_nbs_cluster_markers.csv")  # !Important
write.csv(x = makers_f14_subcluster_nbs, file = file.path(DIR_RES, "objects", "A_NMF30", fname), row.names = F)
makers_f14_subcluster_nbs <- read.csv(file = file.path(DIR_RES, "objects", "A_NMF30", fname))


# Plot marker dotplot 
genes_plot <- subset(makers_f14_subcluster_nbs, p_val_adj<0.05)
genes_plot$updown <- "up"
genes_plot[genes_plot$avg_log2FC<0,"updown"] <- "down"
genes_plot_dp <- genes_plot %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group=TRUE) %>%
  top_n(8, (avg_log2FC))

p3 <- DotPlot(se_nbs, features = rev(unique(genes_plot_dp$gene)), group.by = "seurat_clusters", cols="RdYlBu", scale = T) + 
  theme_dotplot +
  scale_y_discrete(position = "right") +
  coord_flip() +
  labs(title = "F14 nbs cluster marker genes")
  
p <- (p1/p2)|p3;p
png(file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors, "_f14_subclusters_nbs_clusters.png")), width = 9*fig_res, height = 8*fig_res, res = fig_res);p;dev.off()


# Plot marker heatmap 
genes_plot <- subset(makers_f14_subcluster_nbs, p_val_adj<0.05)
genes_plot$updown <- "up"
genes_plot[genes_plot$avg_log2FC<0,"updown"] <- "down"
genes_plot <- genes_plot %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group=TRUE) %>%
  top_n(20, (avg_log2FC))

mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(256)
p <- DoHeatmap(se_nbs, features = c("KRT17", "KRT5", genes_plot$gene), group.by = "seurat_clusters", group.colors =  cluster_cols_nbs, slot = "scale.data") +
  scale_fill_gradientn(colours = rev(mapal))
png(file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors, "_f14_subclusters_nbs_clusters_marker_heatmap.png")), width = 12*fig_res, height = 12*fig_res, res = fig_res);p;dev.off()


mapal <- colorRampPalette(RColorBrewer::brewer.pal(8,"OrRd"))(256)
p <- DoHeatmap(se_nbs, features = c("KRT17", "KRT5", genes_plot$gene), 
               group.by = "seurat_clusters", 
               group.colors =  cluster_cols_nbs, 
               slot = "data") +
  scale_fill_gradientn(colours = mapal)
png(file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors, "_f14_subclusters_nbs_clusters_marker_heatmap_data.png")), width = 12*fig_res, height = 16*fig_res, res = fig_res);p;dev.off()



# Plot marker umap
genes_plot_umap <- genes_plot %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group=TRUE) %>%
  top_n(3, (avg_log2FC))

p <- FeaturePlot(se_nbs, reduction = "umap", features = unique(genes_plot_umap$gene), ncol = 3, cols = col_scale_rocket, pt.size=0.25) & theme_void() & NoLegend() & theme(plot.title = element_text(hjust=0.5, size=10))
png(file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors, "_f14_subclusters_nbs_clusters_marker_umap.png")), width = 6*fig_res, height = 12*fig_res, res = fig_res);p;dev.off()


# Plot marker violin
VlnPlot(se_nbs, 
        features = c(unique(genes_plot_umap$gene), "APOE", "KRT17", "KRT5"), 
        group.by = "seurat_clusters", ncol = 3, pt.size = 0, split.by = "d_c0", split.plot = TRUE, cols = cols_dist) & 
  theme_custom &
  theme(axis.title.y = element_blank(), plot.title = element_text(size=10), axis.text.x = element_text(angle=0))


# Violin plots
p1 <- VlnPlot(se_nbs, 
        features = paste0("factor_", c(6,3,1,21,25,2,12,4,14)), 
        group.by = "seurat_clusters", ncol = 9, pt.size = 0, cols = cluster_cols_nbs) & 
  theme_custom & 
  theme(axis.title.y = element_blank(), plot.title = element_text(size=10), axis.text.x = element_text(angle=0, hjust=0.5))

p2 <- VlnPlot(se_nbs, 
        features = c("SFTPB", "KRT17", "MUC5B", "KRT5", "MMP7", "CXCL10", "APOE", "IGLC3", "HBB"), 
        group.by = "seurat_clusters", ncol = 9, pt.size = 0, cols = cluster_cols_nbs) & 
  theme_custom & 
  theme(axis.title.y = element_blank(), plot.title = element_text(size=10), axis.text.x = element_text(angle=0, hjust=0.5))

p <- p1/p2;p
png(file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors, "_f14_subclusters_nbs_clusters_violin.png")), width = 12*fig_res, height = 4*fig_res, res = fig_res);p;dev.off()


# proportions plots
d <- se_nbs@meta.data %>%
  group_by(sample_name, subject_alias, seurat_clusters) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

p <- ggplot(d, aes(x=sample_name, y=freq, fill=seurat_clusters)) +
  geom_bar(stat = 'identity', colour=NA, position = "stack") +
  facet_grid(~subject_alias, scales = "free_x", space='free') +
  scale_fill_manual(values = cluster_cols_nbs) +
  labs(x="", y="Frequency", fill="", title="F14 subcluster 0 nbs cluster spot proportions") +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        panel.grid = element_blank(),
        legend.position = "right", 
        text = element_text(size=10),
        plot.title = element_text(hjust=0.5));p

png(file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors, "_f14_subclusters_nbs_clusters_proportions.png")), width = 7*fig_res, height = 3.5*fig_res, res = fig_res);p;dev.off()


# Add nbs cluster metadata to se.subset object
se_nbs$f14_nbs_clusters <- se_nbs$seurat_clusters
se.subset <- AddMetaData(se.subset, se_nbs$f14_nbs_clusters, "f14_nbs_clusters")

se.subset$f14_nbs_clusters <- as.character(se.subset$f14_nbs_clusters)
se.subset$f14_nbs_clusters[is.na(se.subset$f14_nbs_clusters)] <- "other"
se.subset$f14_nbs_clusters[se.subset$f14_subclusters=="0"] <- "F14_C0"
se.subset$f14_nbs_clusters <- as.factor(se.subset$f14_nbs_clusters)
summary(as.factor(se.subset$f14_nbs_clusters))

cluster_cols_nbs_all <- setNames(object = c(cluster_cols_nbs, "grey10", "grey90"), nm = levels(as.factor(se.subset$f14_nbs_clusters)))


# Plot all spatial
p <- ST.FeaturePlot(se.subset,
                    features = "f14_nbs_clusters", 
                    ncol=5, cols = cluster_cols_nbs_all, 
                    label.by = "sample_name", 
                    pt.size = 0.75, show.sb = F) + 
  guides(fill = guide_legend(override.aes = list(size=2)))
# p <- ggrastr::rasterize(p, layers = "Point", dpi = fig_res)
png(file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors, "_f14_subclusters_nbs_clusters_spatial.png")), width = 12*fig_res, height = 12*fig_res, res = fig_res);p;dev.off()


# Plot n spots and colors
d <- se.subset@meta.data %>%
  group_by(sample_name, subject_alias, f14_nbs_clusters) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

p <- ggplot(d, aes(x=sample_name, y=freq, fill=f14_nbs_clusters)) +
  geom_bar(stat = 'identity', colour=NA) +
  facet_grid(~subject_alias, scales = "free_x", space='free') +
  scale_fill_manual(values = cluster_cols_nbs_all) +
  labs(x="", y="Frequency", fill="", title="F14 subcluster 0 nbs cluster spot proportions") +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        panel.grid = element_blank(),
        legend.position = "right", 
        text = element_text(size=10),
        plot.title = element_text(hjust=0.5));p
png(file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors, "_f14_subclusters_nbs_clusters_proportions_all.png")), width = 8.5*fig_res, height = 4*fig_res, res = fig_res);p;dev.off()


# Marker dotplot with F14-C0 and others
p <- DotPlot(se.subset, features = rev(unique(genes_plot_dp$gene)), group.by = "f14_nbs_clusters", cols="RdYlBu", scale = T) + 
  theme_dotplot +
  scale_y_discrete(position = "right") +
  coord_flip();p
png(file.path(DIR_FIG_NMF, paste0("hs_visium_A_nmf_1-", n_factors, "_f14_subclusters_nbs_marker_dotplot_other.png")), width = 6*fig_res, height = 8*fig_res, res = fig_res);p;dev.off()



# Plot HE
he_samples_select <- sort(c("IPF_1.TD012.B2.2", "IPF_3.TD031.B1.2", "IPF_2.TD021A.B1.2", "IPF_4.TD042.B3.1")) # "IPF_1.TD012.B2.1", "IPF_4.TD041B.B2.1", "IPF_4.TD042.B3.1"
se.subset_he <- SubsetSTData(se.subset, expression = sample_name %in% he_samples_select)
se.subset_he <- LoadImages(se.subset_he, xdim = 2000)
se.subset_he <- SubsetSTData(se.subset_he, expression = f14_nbs_clusters != "other")

p <- FeatureOverlay(se.subset_he, 
                    features = "f14_nbs_clusters",
                    sampleids = 1:4, 
                    cols = cluster_cols_nbs_all, 
                    pt.alpha = 0.8, 
                    pt.size = 1, 
                    ncol=2, 
                    label.by = "sample_name") & 
  theme(legend.position = "bottom")

fname <- paste0("hs_visium_A_nmf_1-", n_factors, "_f14_subclusters_nbs_clusters_spatial_HE")
pdf(file = file.path(DIR_FIG_NMF, paste0(fname, ".pdf")), width = 10, height = 12, useDingbats = F);p;dev.off()

rm(se.subset_he)



#### Wrap up ####
#' Save rds
fname <- paste0("hs_visium_preproc_A_se_obj_nmf.rds")
saveRDS(se.subset, file = file.path(DIR_RES, "objects", fname))


#' Save metadata as table
fname <- paste0("hs_visium_preproc_A_se_obj_nmf_metadata.csv")
write.csv(se.subset@meta.data, file = file.path(DIR_RES, "objects", fname), row.names = T)
# se_A_md <- read.csv(file.path(DIR_RES, "objects", fname), row.names = 1)


