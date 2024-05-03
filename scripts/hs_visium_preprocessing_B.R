#' [hs_visium_preprocessing_B.R]
#'
#' Preprocessing of spaceranger output from human lung Visium experiments, Workflow B
#' 
#' Processing workflows:
#' - A: Joint HC and IPF integration
#' - B: Individual donor
#'
#' B detailed description
#' Run SCT on donor samples separately, regressing for slide and sample when applicable.
#' Merge and harmonize data based on donor. Generate cluster labels.
#'
#' Aug-Nov 2022, L. Franzén [lovisa.franzen@scilifelab.se]

#### Set up ####
##### Define params. ####
set.seed(1)
SPECIES <- "human"
DIR_ROOT <- getwd()
DIR_DATA <- file.path(DIR_ROOT, "data", SPECIES, "visium")
DIR_RES <- file.path(DIR_ROOT, "results", SPECIES)
fig_res <- 300

##### Load libs #### 
library(STutility)
library(harmony)
library(RColorBrewer)
library(patchwork)
library(dplyr)
library(magrittr)
library(writexl)

##### Other ####
source(file.path(DIR_ROOT, "scripts", "custom_functions.R"))
source(file.path(DIR_ROOT, "scripts", "custom_colors.R"))
theme_custom <- theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))


analysis_txt_fname <- file.path(DIR_RES, paste0("log.hs_visium_preprocessing_B.txt"))
writeLines("Analysis of human lung visium data generated 2021 by script hs_visium_processing_B.R.", analysis_txt_fname)


##### Read tables ####
metadata <- read.table(file.path(DIR_DATA, "hs_visium_metadata_deposit.tsv"), sep = "\t", header = T)
metadata$fibrotic_extent_score_by_pathologist_0.3 <- as.character(metadata$fibrotic_extent_score_by_pathologist_0.3)
rownames(metadata) <- metadata$sample_id


#### Seurat/STutility processing ####
##### Read object ####
#' Use object created in A to process further 
fname <- paste0("hs_visium_preproc_A_se_obj.rds")
se.subset <- readRDS(file = file.path(DIR_RES, "objects", fname))

se_subset_staffli <- GetStaffli(se.subset)
se_subset_split <- SplitObject(se.subset, split.by = "subject_alias")

subject_names <- sort(names(se_subset_split))

##### SCT on obj ####
vars_reg <- c("sample_id", "slide_id")

for(subject in subject_names){
  message(subject)
  set <- se_subset_split[[subject]]
  if (length(unique(set$slide_id))>1) {
    vars_reg <- c("sample_id", "slide_id")
    message(paste0("SCT with vars to reg: ", paste(vars_reg, collapse = ", ")))
    set <- SCTransform(set, vars.to.regress = vars_reg, conserve.memory = T, return.only.var.genes = FALSE)
  } else {
    message("SCT without vars to reg")
    set <- SCTransform(set, conserve.memory = T, return.only.var.genes = FALSE)
  }
  se_subset_split[[subject]] <- set
  rm(set)
}

#' plot QC stats
for(subject in sort(names(se_subset_split))){
  p <- VlnPlot(se_subset_split[[subject]], 
          features=c("nFeature_SCT", "nCount_SCT", "GAPDH", "ACTB"), 
          ncol=2,
          pt.size=0,
          split.by = "tissue_alias", cols = RColorBrewer::brewer.pal(3, "Spectral"),
          group.by = "sample_name") & theme_custom & labs(subtitle = subject);print(p)
  rm(p)
}


##### NMF and PCA-Harmony-clustering of each set ####
for(subject in subject_names){
  message(subject)
  set <- se_subset_split[[subject]]
  n_dims <- 30
  cluster_res <- 0.4
  set <- RunNMF(set, nfactors = n_dims)
  set <- RunUMAP(set, reduction = "NMF", dims = 1:n_dims, reduction.name = "umap.nmf")
  set <- RunPCA(set, npcs = 50) %>%
    RunUMAP(reduction="pca", dims=1:n_dims)
  if (length(unique(set$slide_id))>1) {
    vars_reg <- c("slide_id", "tissue_alias")
    message(paste0("Harmony vars to reg: ", paste(vars_reg, collapse = ", ")))
    set <- RunHarmony(set,
                      group.by.vars = vars_reg,
                      reduction = "pca",
                      assay.use = "SCT",
                      epsilon.cluster=-Inf,  # prevent early stopping
                      epsilon.harmony=-Inf,  # prevent early stopping
                      max.iter.harmony = 21,
                      plot_convergence = T,
                      verbose = T)
    set <- RunUMAP(set, reduction = "harmony", dims = 1:n_dims, reduction.name = "umap.harmony")
    message("Clustering with Harmony-PCA")
    set <- FindNeighbors(set, reduction = "harmony", dims = 1:n_dims) %>%
      FindClusters(resolution=cluster_res)
  } else if (length(unique(set$tissue_alias))>1){
    vars_reg <- c("tissue_alias")
    message(paste0("Harmony vars to reg: ", paste(vars_reg, collapse = ", ")))
    set <- RunHarmony(set,
                      group.by.vars = vars_reg,
                      reduction = "pca",
                      assay.use = "SCT",
                      epsilon.cluster=-Inf,  # prevent early stopping
                      epsilon.harmony=-Inf,  # prevent early stopping
                      max.iter.harmony = 21,
                      plot_convergence = F,
                      verbose = T)
    set <- RunUMAP(set, reduction = "harmony", dims = 1:n_dims, reduction.name = "umap.harmony")
    message("Clustering with Harmony-PCA")
    set <- FindNeighbors(set, reduction = "harmony", dims = 1:n_dims) %>%
      FindClusters(resolution=cluster_res)
  }
  else {
    message("No Harmony")
    message("Clustering with PCA")
    set <- FindNeighbors(set, reduction = "pca", dims = 1:n_dims) %>%
      FindClusters(resolution=cluster_res)
  }
  se_subset_split[[subject]] <- set
  rm(set)
}

de_markers_list <- list()
for(subject in subject_names){
  message(subject)
  de_markers <- FindAllMarkers(se_subset_split[[subject]], only.pos = TRUE, logfc.threshold = 0.3, verbose = T)
  de_markers_list[[subject]] <- de_markers
  rm(de_markers)
}
names(de_markers_list) <- subject_names

write_xlsx(
  x = de_markers_list,
  path = file.path(DIR_RES, "objects", paste0("hs_visium_preproc_B_cluster_markers.xlsx")),
  col_names = TRUE,
  format_headers = TRUE,
  use_zip64 = FALSE
)


#### Plot and object export ####
#####  Plot UMAP and Harmony ##### 
p_list <- list()
for(subject in subject_names){
  size_pt <- 0.6
  if (length(unique(se_subset_split[[subject]]$tissue_alias))>1) {
    p1 <- FeaturePlot(se_subset_split[[subject]], reduction = "umap.harmony", features = "nFeature_SCT", 
                      pt.size = size_pt, cols = col_scale_rocket, raster = T)
    p2 <- DimPlot(se_subset_split[[subject]], reduction = "umap.harmony", group.by = "seurat_clusters", 
                  pt.size = size_pt, label = T, cols = cols_d3_20, raster = T)
    p3 <- DimPlot(se_subset_split[[subject]], reduction = "umap.harmony", group.by = "sample_name", 
                  pt.size = size_pt, cols = cols_d3_20, raster = T) & 
      theme(legend.position = "bottom", legend.text = element_text(size=5))
  } else {
    p1 <- FeaturePlot(se_subset_split[[subject]], reduction = "umap", features = "nFeature_SCT", 
                      pt.size = size_pt, cols = col_scale_rocket, raster = T)
    p2 <- DimPlot(se_subset_split[[subject]], reduction = "umap", group.by = "seurat_clusters", 
                  pt.size = size_pt, label = T, cols = cols_d3_20, raster = T)
    p3 <- DimPlot(se_subset_split[[subject]], reduction = "umap", group.by = "sample_name", 
                  pt.size = size_pt, cols = cols_d3_20, raster = T) & 
      theme(legend.position = "bottom", legend.text = element_text(size=5))
  }
  p_list[[subject]] <- (p1+p2+p3+patchwork::plot_annotation(title=subject))
}
pdf(file = file.path(DIR_RES, "figures", "hs_visium_preproc_B.umap_Nov22.pdf"), 
    width = 16, height = 12)
cowplot::plot_grid(plotlist = p_list[1:4], ncol=4)
cowplot::plot_grid(plotlist = p_list[5:8], ncol=4)
dev.off()


p_list <- list()
for(subject in subject_names){
  size_pt <- 0.6
  p1 <- FeaturePlot(se_subset_split[[subject]], reduction = "umap.nmf", features = "nFeature_SCT",
                    pt.size = size_pt, cols = col_scale_rocket, raster = T)
  p2 <- DimPlot(se_subset_split[[subject]], reduction = "umap.nmf", group.by = "seurat_clusters", 
                pt.size = size_pt, label = T, cols = cols_d3_20, raster = T)
  p3 <- DimPlot(se_subset_split[[subject]], reduction = "umap.nmf", group.by = "sample_name", 
                pt.size = size_pt, cols = cols_d3_20, raster = T) & 
    theme(legend.position = "bottom", legend.text = element_text(size=5))
  p_list[[subject]] <- (p1+p2+p3+patchwork::plot_annotation(title=subject))
}
pdf(file = file.path(DIR_RES, "figures", "hs_visium_preproc_B.umap_nmf_Nov22.pdf"), 
    width = 16, height = 12)
cowplot::plot_grid(plotlist = p_list[1:4], ncol=4)
cowplot::plot_grid(plotlist = p_list[5:8], ncol=4)
dev.off()


#####  Plot Cluster markers ##### 
# DotPlot(se.f14high, features = rev(genes_plot_dp$gene), group.by = "f14_subclusters", cols="RdYlBu", scale = T) + 
#   theme_dotplot +
#   scale_y_discrete(position = "right") +
#   coord_flip();
p_list <- list()
for(subject in names(de_markers_list)){
  top3_makers <- de_markers_list[[subject]]%>%
    dplyr::filter(p_val_adj < 0.01) %>%
    dplyr::group_by(cluster) %>%
    dplyr::top_n(wt = avg_log2FC, n = 3)
  p <- DotPlot(se_subset_split[[subject]], features = unique(top3_makers$gene), group.by = "seurat_clusters", cols="RdYlBu") +
    theme_dotplot +
    coord_flip() +
    labs(title=paste0(subject, " top 3 cluster markers")) +
    theme(axis.title = element_blank(), 
          legend.position = "none", 
          plot.title = element_text(face = "bold", size = 11, colour = "black"))
  p_list[[subject]] <- p # (p+patchwork::plot_annotation(title=subject))
  rm(p)
}
pdf(file = file.path(DIR_RES, "figures", "hs_visium_preproc_B.cluster_markers_Nov22.pdf"), 
    width = 17, height = 12)
patchwork::wrap_plots(p_list[1:4], nrow = 1)/patchwork::wrap_plots(p_list[5:8], nrow = 1)+patchwork::plot_layout(heights = c(1,3))
dev.off()


#####  Plot Clusters and Histopath Annotation ##### 
p_list <- list()
for(subject in subject_names){
  message(subject)
  set_all <- se_subset_split[[subject]]
  set <- SubsetSTData(set_all, expression = annotation!= "NA")
  
  p1 <- ST.FeaturePlot(set, features = "seurat_clusters", pt.size = 1.3, pt.border = F, 
                       cols = cols_d3_20, label.by = "sample_name", ncol=3, dark.theme = F) + 
    labs(color="Clusters") +
    guides(fill=guide_legend(ncol=2, override.aes = list(size=2))) +
    theme(legend.position = "right", plot.title = element_blank(), aspect.ratio = 1)
  p2 <- ST.FeaturePlot(set, features = "annotation", pt.size = 1.3, pt.border = F, 
                       cols = cols_d3_20, label.by = "sample_name", ncol=3, dark.theme = F) + 
    labs(color="Annotation") +
    guides(fill=guide_legend(override.aes = list(size=2))) +
    theme(legend.position = "right", plot.title = element_blank(), aspect.ratio = 1)
  p_spat <- p1 / p2 & patchwork::plot_annotation(title = paste0(subject, " spatial distribution"))
  p_spat <- ggrastr::rasterize(p_spat, layers = "Point", dpi = 300)
  
  anno_df1 <- set_all[[]] %>%
    group_by(seurat_clusters, sample_name) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n),
           pct = round(100 * n/sum(n), 0),
           rel.freq = paste0(round(100 * n/sum(n), 0), "%")) %>%
    as.data.frame() %>%
    arrange(sample_name, desc(freq))
  
  anno_df2 <- set[[]] %>%
    group_by(seurat_clusters, annotation) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n),
           pct = round(100 * n/sum(n), 0),
           rel.freq = paste0(round(100 * n/sum(n), 0), "%")) %>%
    as.data.frame() %>%
    arrange(seurat_clusters, desc(freq))
  
  anno_df3 <- set[[]] %>%
    group_by(annotation, seurat_clusters) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n),
           pct = round(100 * n/sum(n), 0),
           rel.freq = paste0(round(100 * n/sum(n), 0), "%")) %>%
    as.data.frame() %>%
    arrange(annotation, desc(freq))
  
  p1 <- ggplot(anno_df1, aes(x=seurat_clusters, y=freq, fill=sample_name)) +
    geom_bar(stat = 'identity', colour=NA, position = "stack") +
    scale_fill_manual(values = cols_d3_20) +
    labs(x="", y="Frequency", fill="Sample", title="Cluster proportions\n in samples") +
    theme_linedraw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "left", 
      legend.text = element_text(size=8),
      text = element_text(size=8),
      plot.title = element_text(hjust=0.5)) +
    coord_flip()
  
  p2 <- ggplot(anno_df2, aes(x=seurat_clusters, y=freq, fill=annotation)) +
    geom_bar(stat = 'identity', colour=NA, position = "stack") +
    scale_fill_manual(values = cols_d3_20) +
    labs(x="", y="Frequency", fill="Annotation", title="Annotation proportions\n in clusters") +
    theme_linedraw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom", 
      legend.text = element_text(size=8),
      text = element_text(size=8),
      plot.title = element_text(hjust=0.5)) +
    coord_flip();
  
  p3 <- ggplot(anno_df3, aes(x=annotation, y=freq, fill=seurat_clusters)) +
    geom_bar(stat = 'identity', colour=NA, position = "stack") +
    scale_fill_manual(values = cols_d3_20) +
    labs(x="", y="Frequency", fill="Clusters", title="Cluster proportions\n in annotated regions") +
    theme_linedraw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom", 
      legend.text = element_text(size=8),
      text = element_text(size=8),
      plot.title = element_text(hjust=0.5)) +
    coord_flip();
  
  p_freq <- (p1 | p2 | p3) & patchwork::plot_annotation(title = paste0(subject, " cluster proportions"))
  
  p_list[[subject]] <- cowplot::plot_grid(plotlist = list(p_spat, p_freq), ncol=1, rel_heights = c(5,4))
  rm(set_all);rm(set)
}

pdf(file = file.path(DIR_RES, "figures", "hs_visium_preproc_B.cluster_dist.pdf"), 
    width = 11, height = 11)
for( subject in subject_names ) { print(p_list[[subject]]) }
dev.off()

# pdf(file = file.path(DIR_RES, "figures", "hs_visium_preproc_B.cluster_dist_test.pdf"),
#     width = 11, height = 11)
# cowplot::plot_grid(plotlist = list(p_spat, p_freq), ncol=1, rel_heights = c(5,4))
# dev.off()

##### NMF ##### 
#' Can also be run using the script 'hs_visium_preprocessing_B.R'
n_factors <- 30
DIR_FIG_NMF <- file.path(DIR_RES, "figures", paste0("nmf_preproc_B_", n_factors, "_factors"))
# dir.create(path = DIR_FIG_NMF)

# DIR_OBJ_NMF <- file.path(DIR_RES, "objects", paste0("nmf_preproc_B"))
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

for(subject in subject_names[7:8]){
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
                       pt.size = 1.5, 
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
}
rm(se_he_ipf)



###### Immune factors plot #####
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
                      pt.size = 1.5, 
                      show.sb = F);p
  p <- ggrastr::rasterize(p, layers = "Point", dpi = 300)
  print(p)
}
dev.off()



#### Wrap up ####
#' Save rds
fname <- paste0("hs_visium_preproc_B_se_obj_list.rds")
saveRDS(se_subset_split, file = file.path(DIR_RES, "objects", fname))
# se_subset_split <- readRDS(file = file.path(DIR_RES, "objects", fname))

fname <- paste0("hs_visium_preproc_B_cluster_markers_list.rds")
saveRDS(de_markers_list, file = file.path(DIR_RES, "objects", fname))
# de_markers_list <- readRDS(file = file.path(DIR_RES, "objects", fname))

#' Session info
cat(capture.output(sessionInfo()), file = analysis_txt_fname, sep="\n", append=TRUE)


# renv::snapshot()
