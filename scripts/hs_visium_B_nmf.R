#' [hs_visium_B_nmf.R]
#'
#' Run, analyse, and plot NMF results for hs visium B dataset
#'
#'
#' Oct 2022, L. Franzén [lovisa.franzen@scilifelab.se]

#### Set up ####
##### Define params. ####
set.seed(1)
SPECIES <- "human"
DIR_ROOT <- getwd()
DIR_DATA <- file.path(DIR_ROOT, "data", SPECIES, "visium")
DIR_RES <- file.path(DIR_ROOT, "results", SPECIES)
DIR_FIG <- file.path(DIR_RES, "figures")
fig_res <- 300

##### Load libs ####
library(STutility)
library(harmony)
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

##### Read objects ####
fname <- paste0("hs_visium_preproc_B_se_obj_list.rds")
se_subset_split <- readRDS(file = file.path(DIR_RES, "objects", fname))
subject_names <- sort(names(se_subset_split))


##### NMF ##### 
n_factors <- 30
DIR_FIG_NMF <- file.path(DIR_RES, "figures", paste0("nmf_preproc_B_", n_factors, "_factors"))
# dir.create(path = DIR_FIG_NMF)

# DIR_OBJ_NMF <- file.path(DIR_RES, "objects", paste0("nmf_preproc_B"))  # previous folder name
# dir.create(path = DIR_OBJ_NMF)


###### Gene loadings list ######
factor_gene_loadings_list <- lapply(subject_names, function(subject){
  factor_gene_loadings <- lapply(1:n_factors, function(factor_x){
    message(paste("Factor", factor_x))
    feat_loads <- as.data.frame(se_subset_split[[subject]]@reductions$NMF@feature.loadings[,paste0("factor_",factor_x)])
    colnames(feat_loads) <- "gene_loading"
    feat_loads$gene_loading_scaled <- Scale01(feat_loads$gene_loading)
    feat_loads$gene <- rownames(feat_loads)
    feat_loads$factor <- factor_x
    top_n_genes <- 100
    feat_loads <- feat_loads %>%
      dplyr::slice_max(order_by = gene_loading, n = top_n_genes) %>%
      mutate(rank = dense_rank(desc(gene_loading))) %>%
      as.data.frame()
  })
  factor_gene_loadings_m <- bind_rows(factor_gene_loadings)
})
names(factor_gene_loadings_list) <- subject_names

write_xlsx(
  x = factor_gene_loadings_list,
  path = file.path(DIR_RES, "objects", paste0("hs_visium_B_nmf_1-", n_factors, "_top100_gene_loadings.xlsx")),
  col_names = TRUE,
  format_headers = TRUE,
  use_zip64 = FALSE
)


###### Gene loading plot ######
pdf(file = file.path(DIR_FIG_NMF, paste0("hs_visium_B_nmf_1-", n_factors, "_gene_loadings", ".pdf")), 
    width = 9, height = 8)
for(subject in subject_names){
  message(subject)
  factor_gene_loadings_list[[subject]]
  d_plot <- factor_gene_loadings_list[[subject]]
  # d_plot$gene_loading_scaled <- 0  # if scale with only top 100
  # for(f in 1:n_factors){
  #   d_plot[d_plot$factor==f, "gene_loading_scaled"] <- Scale01(d_plot[d_plot$factor==f, "gene_loading"])
  # }
  p <- ggplot() +
    geom_line(data = d_plot, mapping = aes(x=rank, y=gene_loading_scaled)) +
    geom_point(data = subset(d_plot, rank<11), mapping = aes(x=rank, y=gene_loading_scaled), color="black", size=.5) +
    facet_wrap(~factor) +
    labs(title = subject) +
    theme_linedraw() +
    theme(panel.grid = element_blank())
  
  # png(file = file.path(DIR_FIG_NMF, paste0("hs_visium_B_nmf_1-", n_factors, "_gene_loadings_", subject, ".png")), 
  #     width = 9*fig_res, height = 8*fig_res, res = fig_res)
  print(p)
  # dev.off()
}
dev.off()


###### Spatial nmf plot ######
for(subject in subject_names){
  n_samples <- length(unique(se_subset_split[[subject]]$sample_name))
  if(n_samples>3){
    fig_sample_height <- 3
  } else {
    fig_sample_height <- 2
  }
  if(n_samples>1){
    fig_sample_width <- 24
  } else {
    fig_sample_width <- 12
  }
  png(file = file.path(DIR_FIG_NMF, paste0("hs_visium_B_nmf_1-", n_factors, "_", subject, ".png")), 
      width = fig_sample_width*fig_res, height = fig_sample_height*5*fig_res, res = fig_res)
  p <- ST.DimPlot(se_subset_split[[subject]], 
                  dims = 1:30,
                  # min.cutoff = 0.2,
                  max.cutoff = 6,
                  label.by = "sample_name",
                  center.zero = F, 
                  reduction = "NMF", 
                  ncol = 3, 
                  grid.ncol = 6,
                  pt.size = 0.55, 
                  pt.border = F,
                  cols = rev(RColorBrewer::brewer.pal(10, "Spectral")), 
                  dark.theme = T,
                  show.sb = F)
  print(p)
  dev.off()
}


###### NMF spatial HE ######
sample_name_select <- list(IPF_1 = c("IPF_1.TD011.B1.1", "IPF_1.TD012.B2.1", "IPF_1.TD013.B3.1"),
                           IPF_2 = c("IPF_2.TD021A.B1.1", "IPF_2.TD021B.B2.1", "IPF_2.TD022.B3.1"),
                           IPF_3 = c("IPF_3.TD031.B1.2", "IPF_3.TD032A.B2.2", "IPF_3.TD032B.B3.2"),
                           IPF_4 = c("IPF_4.TD041A.B1.1", "IPF_4.TD041B.B2.1", "IPF_4.TD042.B3.1"))

###### Fibrotic factors plot ######
fib_genes <- "FN1"
for(subject in subject_names){
  fib_genes <- c(fib_genes, grep("^COL[0-9]", rownames(se_subset_split[[subject]]@reductions$NMF@feature.loadings), value = T), "FN1", "FBN1", "THY1", "TNC", "VIM")
}
fib_genes <- sort(unique(fib_genes))
factor_select_list <- list()
for(subject in subject_names){
  nmf_gene_loads <- se_subset_split[[subject]]@reductions$NMF@feature.loadings
  nmf_gene_loads <- nmf_gene_loads[rownames(nmf_gene_loads)%in%fib_genes,]
  fib_factor_sum <- sort(colSums(nmf_gene_loads), decreasing = T)
  factor_select_list[[subject]] <- names(fib_factor_sum[as.numeric(fib_factor_sum)>0.5])
}

for(subject in subject_names[5:6]){
  message(subject)
  se_he_ipf <- SubsetSTData(se_subset_split[[subject]], expression = sample_name %in% sample_name_select[[subject]])
  se_he_ipf <- LoadImages(se_he_ipf, xdim = 1e3)
  factors_plot <- head(factor_select_list[[subject]],3)
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
    print(fac)
    d_plot <- subset(factor_gene_loadings_list[[subject]], factor == as.numeric(gsub("factor_", "", fac)) & rank <= 40)
    ggplot(d_plot, aes(x=reorder(gene, rank), y=gene_loading_scaled)) +
      geom_point() +
      labs(x="", y="", title=fac) +
      theme_classic() +
      theme(axis.text.x = element_text(angle=45, hjust=1), 
            axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
            panel.grid = element_blank(), panel.grid.major.x = element_line(size = .3, color="black"))
  })
  p2 <- patchwork::wrap_plots(p_list, nrow = 1)
  p <- p1/p2+patchwork::plot_layout(heights = c(5,1))
  pdf(file = file.path(DIR_FIG_NMF, paste0("hs_visium_B_nmf_1-", n_factors,"_HE_fib_factors_", subject, ".pdf")), width = 18, height = 9, useDingbats = F);print(p);dev.off()
  rm(se_he_ipf)
  }



###### Collagens & Keratins ######
sum_factor_cutoff <- 2 # _all: 0
sum_gene_cutoff <- 0.3
gene_add <- "KRT5"
nmf_gene_loads_list <- lapply(subject_names, function(subject){
  fib_genes <- sort(c(grep("^COL[0-9]", rownames(se_subset_split[[subject]]@reductions$NMF@feature.loadings), value = T), 
                      "FN1", "FBN1", "THY1", "TNC", "VIM", "ACTA2", "TAGLN", "LUM", "DCN", "PDGRFA", "SFTPC", "SFTPB",
                      grep("^KRT[0-9]", rownames(se_subset_split[[subject]]@reductions$NMF@feature.loadings), value = T)))
  nmf_gene_loads_scaled <- apply(se_subset_split[[subject]]@reductions$NMF@feature.loadings, 2, Scale01)
  nmf_gene_loads_scaled <- as.data.frame(t(nmf_gene_loads_scaled[rownames(nmf_gene_loads_scaled)%in%fib_genes,]))
  nmf_gene_loads_scaled$sum <- as.numeric(rowSums(nmf_gene_loads_scaled))
  nmf_gene_loads_scaled_filtered <- nmf_gene_loads_scaled[nmf_gene_loads_scaled$sum>sum_factor_cutoff, fib_genes[fib_genes%in%colnames(nmf_gene_loads_scaled)]]
  # nmf_gene_loads_scaled_filtered <- nmf_gene_loads_scaled_filtered[, c(colnames(nmf_gene_loads_scaled_filtered)[colSums(nmf_gene_loads_scaled_filtered)>sum_gene_cutoff])]
  # genes_keep <- colnames(nmf_gene_loads_scaled_filtered)
  # nmf_gene_loads_scaled_filtered <- nmf_gene_loads_scaled_filtered[,-grep("\\.[0-9]", colnames(nmf_gene_loads_scaled_filtered))]
  })
names(nmf_gene_loads_list) <- subject_names

pdf(file = file.path(DIR_FIG_NMF, paste0("hs_visium_B_nmf_1-", n_factors,"_selected_gene_heatmap", ".pdf")), width = 10, height =4, useDingbats = F)
for(subject in subject_names){
  message(subject)
  if(nrow(nmf_gene_loads_list[[subject]])>2){
    message("plot")
    p <- pheatmap(nmf_gene_loads_list[[subject]], display_numbers = F, color = RColorBrewer::brewer.pal(8,"BuPu"), main = subject)
    print(p)
  } else if (nrow(nmf_gene_loads_list[[subject]])>0) {
    p <- pheatmap(nmf_gene_loads_list[[subject]], display_numbers = F, color = RColorBrewer::brewer.pal(8,"BuPu"), cluster_rows = F, main = subject)
    print(p)
  }
}
dev.off()


factors_plot_list <- lapply(subject_names, function(subject){
  f_sum <- data.frame(factor_sum=as.numeric(rowSums(nmf_gene_loads_list[[subject]])),
                      factor=rownames(nmf_gene_loads_list[[subject]])
                      )
  f_sum$subject <- rep(subject, nrow(f_sum))
  f_sum <- f_sum[order(f_sum$factor_sum,decreasing = T),]
})
names(factors_plot_list) <- subject_names


#' Plot on HE
sample_name_select <- list(IPF_1 = c("IPF_1.TD011.B1.1", "IPF_1.TD012.B2.1", "IPF_1.TD013.B3.1"),
                           IPF_2 = c("IPF_2.TD021A.B1.1", "IPF_2.TD021B.B2.1", "IPF_2.TD022.B3.1"),
                           IPF_3 = c("IPF_3.TD031.B1.2", "IPF_3.TD032A.B2.2", "IPF_3.TD032B.B3.2"),
                           IPF_4 = c("IPF_4.TD041A.B1.1", "IPF_4.TD041B.B2.1", "IPF_4.TD042.B3.1"))

for(subject in names(sample_name_select)){
  se_he <- se_subset_split[[subject]]
  se_he <- SubsetSTData(se_he, expression = sample_name %in% sample_name_select[[subject]])
  se_he <- LoadImages(se_he, xdim = 1e3)
  
  f_plot <- factors_plot_list[[subject]] %>% head(4) %>% `[[`(.,'factor')
  
  pdf(file = file.path(DIR_FIG_NMF, paste0("hs_visium_B_nmf_1-", n_factors,"_HE_fibkrt_factors_", subject, ".pdf")), width = 14, height = 6, useDingbats = F)
  for(f in f_plot){
    p <- FeatureOverlay(se_he, 
                        sampleids = 1:3,
                        features = f, 
                        ncols = 3, 
                        value.scale = "all",
                        label.by = "sample_name",
                        cols = rev(col_scale_mako),
                        add.alpha = T,
                        pt.size = 1.2, 
                        show.sb = F) & theme(legend.position = "bottom")
    p <- ggrastr::rasterize(p, layers = "Point", dpi = 300)
    print(p)
  }
  dev.off()
}


###### Selected factors per donor ######
factor_select <- list(HC_3 = paste0("factor_", c(13, 9, 19)),
                      IPF_1 = paste0("factor_", c(6, 7, 16)),
                      IPF_2 = paste0("factor_", c(3, 29, 22)),
                      IPF_3 = paste0("factor_", c(8, 12, 5)),
                      IPF_4 = paste0("factor_", c(6, 15, 5)))

for(subject in names(factor_select)[3]){
  se_he <- se_subset_split[[subject]]
  if(subject %in% names(sample_name_select)){
    se_he <- SubsetSTData(se_he, expression = sample_name %in% sample_name_select[[subject]])
  } else {
    se_he <- SubsetSTData(se_he, expression = subject_alias %in% subject)
  }
  se_he <- LoadImages(se_he, xdim = 1e3)

  pdf(file = file.path(DIR_FIG_NMF, paste0("hs_visium_B_nmf_1-", n_factors,"_HE_fibkrt_factors_blend_", subject, ".pdf")), width = 14, height = 5, useDingbats = F)
  p1 <- FeatureOverlay(se_he, 
                 sampleids = 1:3,
                 features = factor_select[[subject]], 
                 ncols = 3,
                 value.scale = "all",
                 label.by = "sample_name",
                 add.alpha = T,
                 blend = T,
                 pt.size = 1.2, 
                 show.sb = F)
  p1 <- ggrastr::rasterize(p1, layers = "Point", dpi = 300)
  print(p1)
  for(i in 1:3){
    p2 <- FeatureOverlay(se_he, 
                         sampleids = i,
                         features = factor_select[[subject]], 
                         ncols = 3,
                         label.by = "sample_name",
                         cols = rev(col_scale_mako),
                         max.cutoff = 6,
                         add.alpha = T, 
                         pt.size = 1.2, 
                         show.sb = F) & NoLegend()
    p2 <- ggrastr::rasterize(p2, layers = "Point", dpi = 300)
    print(p2)
  }
  dev.off()
}


###### Selected immune factors plot #####
subject_plot <- "IPF_2"
factor_plot <- "factor_15"

se_he <- se_subset_split[[subject_plot]]
se_he <- SubsetSTData(se_he, expression = subject_alias == subject_plot)
se_he <- LoadImages(se_he, xdim = 2e3)

sample_ids_selected <- as.numeric(unique(se_he@tools$Staffli@meta.data[colnames(se_he),]$sample))

factor_max <- round(as.numeric(summary(se_he@reductions$NMF@cell.embeddings[,factor_plot])[6])-3)
factor_min <- round(as.numeric(summary(se_he@reductions$NMF@cell.embeddings[,factor_plot])[2])+.5)

pdf(file = file.path(DIR_FIG_NMF, paste0("hs_visium_B_nmf_1-", n_factors,"_HE_", subject_plot, "_", factor_plot, ".pdf")), width = 10, height = 10, useDingbats = F)
for(i in sample_ids_selected){  #:length(unique(se_he$sample_name))
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


###### TLS & B cell factors  #####
genes_tls <- c("CXCL13", "IL7R", "CXCR4", "PTPRC")
genes_bcell <- c("IGHG3", "IGHG1", "IGLC2", "IGLC3", "JCHAIN")
  
tls_factors <- list()
for(subject in grep("IPF", subject_names, value = T)){
  nmf_gene_loads_scaled <- apply(se_subset_split[[subject]]@reductions$NMF@feature.loadings, 2, Scale01)
  tls_factors[[subject]] <- nmf_gene_loads_scaled[rownames(nmf_gene_loads_scaled) %in% genes_tls,] %>% colSums() %>% sort(decreasing = T) %>% head(1) %>% names()
}

bcell_factors <- list()
for(subject in grep("IPF", subject_names, value = T)){
  nmf_gene_loads_scaled <- apply(se_subset_split[[subject]]@reductions$NMF@feature.loadings, 2, Scale01)
  bcell_factors[[subject]] <- nmf_gene_loads_scaled[rownames(nmf_gene_loads_scaled) %in% genes_bcell,] %>% colSums() %>% sort(decreasing = T) %>% head(1) %>% names()
}

# tls_factors <- list(IPF_1 = "factor_10", IPF_2 = "factor_15", IPF_3 = "factor_17")
# bcell_factors <- list(IPF_1 = "factor_22", IPF_2 = "factor_10", IPF_3 = "factor_1", IPF_4 = "factor_14")


#' TLS-factor gene loadings -- moved to script hs_visium_manus_figs
# pdf(file = file.path(DIR_FIG_NMF, paste0("hs_visium_B_nmf_1-", n_factors,"_tls_factor_gene_hm", ".pdf")), width = 12, height = 3, useDingbats = F)

#' Plot on HE for selected images -- moved to script hs_visium_manus_figs


#' B-cell-factor gene loadings
bcell_f_genes <- lapply(names(bcell_factors), function(s){
  se_subset_split[[s]]@reductions$NMF@feature.loadings[,bcell_factors[[s]]]
})
gene_ovelap <- Reduce(intersect, lapply(1:length(bcell_f_genes), function(i){names(bcell_f_genes[[i]])}))

bcell_f_df <- data.frame(gene = gene_ovelap,
                         IPF_1.F22 = bcell_f_genes[[1]][gene_ovelap],
                         IPF_2.F10 = bcell_f_genes[[2]][gene_ovelap],
                         IPF_3.F1 = bcell_f_genes[[3]][gene_ovelap],
                         IPF_4.F14 = bcell_f_genes[[4]][gene_ovelap])

bcell_f_df[,2:ncol(bcell_f_df)] <- apply(bcell_f_df[,2:ncol(bcell_f_df)], 2, Scale01)
bcell_f_df_filt <- bcell_f_df
bcell_f_df_filt$gene <- NULL
rsums <- rowSums(bcell_f_df_filt)
cutoff <- as.numeric(quantile(rsums, .95)[1])
bcell_f_df_filt <- bcell_f_df_filt[rsums>cutoff,]
head(bcell_f_df_filt);dim(bcell_f_df_filt)

pdf(file = file.path(DIR_FIG_NMF, paste0("hs_visium_B_nmf_1-", n_factors,"_bcell_factor_gene_hm", ".pdf")), width = 12, height = 3, useDingbats = F)
print(pheatmap(t(bcell_f_df_filt), 
               display_numbers = F, 
               cluster_cols = F, 
               treeheight_row = 20, 
               color = col_scale_mako,
               cellwidth = 14, 
               cellheight = 14, 
               angle_col = 90))
dev.off()


###### All factors for selected sample  #####
#' IPF 3
subject <- "IPF_3"
sample_select <- "IPF_3.TD031.B1.2"
se_he <- se_subset_split[[subject]]
se_he <- SubsetSTData(se_he, expression = sample_name %in% sample_select)
se_he <- LoadImages(se_he, xdim = 600)

pdf(file = file.path(DIR_FIG_NMF, paste0("hs_visium_B_nmf_1-", n_factors,"_all_factors_", sample_select, ".pdf")), width = 16, height = 9, useDingbats = F)
for(facs in list(1:10, 11:20, 21:30)){
  p <- FeatureOverlay(se_he,
                      features = paste0("factor_", facs), 
                      ncols = 5,
                      label.by = "sample_name",
                      cols = rev(col_scale_mako),
                      add.alpha = T,
                      pt.size = 0.8, 
                      show.sb = F) & theme(legend.position = "bottom")
  p <- ggrastr::rasterize(p, layers = "Point", dpi = 300)
  print(p)
}
dev.off()


#' IPF 2
subject <- "IPF_2"

sample_select <- "IPF_2.TD021B.B2.1"
se_he <- se_subset_split[[subject]]
se_he <- SubsetSTData(se_he, expression = sample_name %in% sample_select)
se_he <- LoadImages(se_he, xdim = 600)

pdf(file = file.path(DIR_FIG_NMF, paste0("hs_visium_B_nmf_1-", n_factors,"_all_factors_", sample_select, ".pdf")), width = 16, height = 9, useDingbats = F)
for(facs in list(1:10, 11:20, 21:30)){
  p <- FeatureOverlay(se_he,
                      features = paste0("factor_", facs), 
                      ncols = 5,
                      label.by = "sample_name",
                      cols = rev(col_scale_mako),
                      add.alpha = T,
                      pt.size = 0.8, 
                      show.sb = F) & theme(legend.position = "bottom")
  p <- ggrastr::rasterize(p, layers = "Point", dpi = 300)
  print(p)
}
dev.off()

sample_select <- "IPF_2.TD022.B3.1"
se_he <- se_subset_split[[subject]]
se_he <- SubsetSTData(se_he, expression = sample_name %in% sample_select)
se_he <- LoadImages(se_he, xdim = 600)

pdf(file = file.path(DIR_FIG_NMF, paste0("hs_visium_B_nmf_1-", n_factors,"_all_factors_", sample_select, ".pdf")), width = 16, height = 9, useDingbats = F)
for(facs in list(1:10, 11:20, 21:30)){
  p <- FeatureOverlay(se_he,
                      features = paste0("factor_", facs), 
                      ncols = 5,
                      label.by = "sample_name",
                      cols = rev(col_scale_mako),
                      add.alpha = T,
                      pt.size = 0.8, 
                      show.sb = F) & theme(legend.position = "bottom")
  p <- ggrastr::rasterize(p, layers = "Point", dpi = 300)
  print(p)
}
dev.off()

