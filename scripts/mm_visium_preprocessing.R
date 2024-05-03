#' [mm_visium_preprocessing.R]
#'
#' Preprocessing of spaceranger output from mouse lung Visium experiments
#'
#'
#' Aug 2022, L. Franzén [lovisa.franzen@scilifelab.se]

#### Set up ####
##### Define params. ####
set.seed(1)
SPECIES <- "mouse"
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

##### Other ####
source(file.path(DIR_ROOT, "scripts", "custom_functions.R"))
source(file.path(DIR_ROOT, "scripts", "custom_colors.R"))
theme_custom <- theme(axis.title.x = element_blank())

analysis_txt_fname <- file.path(DIR_RES, paste0("log.mm_visium_preprocessing.txt"))
writeLines("Analysis of mouse lung visium data generated 2021 by script mm_visium_preprocessing.R.", analysis_txt_fname)


##### Read tables ####
gene_anno <- read.table(file = file.path(DIR_ROOT, "data/misc", "mm_gene_ensembl.2022-05-23.tsv"), 
                        sep = "\t", 
                        header = T, 
                        stringsAsFactors = F)

metadata <- read.table(file.path(DIR_DATA, "mm_visium_metadata.tsv"), sep = "\t", header = T)
rownames(metadata) <- metadata$sample_id
metadata$n <- 1:nrow(metadata)

histopath_anno <- read.csv(file.path(DIR_DATA, "mm_visium_merged_histo_annotations.csv"), row.names = 1)
qc_spot_rm <- read.csv(file.path(DIR_DATA, "mm_visium_qc_spot_remove.csv"), row.names = 1)


#### Create Seurat/STutility object ####
infoTable <- createInfoTable(data.dir = DIR_DATA, 
                             sample.select.ids = metadata$sample_id, 
                             metadata.df = metadata, 
                             subfolder.name = NA,
                             use.raw = F)

se <- InputFromTable(infoTable, 
                     minUMICountsPerGene = 100,
                     minSpotsPerGene = 5, 
                     minUMICountsPerSpot = 200, 
                     platform = "Visium")

##### Histopath annotations #### 
#' Add annotations by histopathologist to metadata
histopath_anno_add <- histopath_anno$annotation
names(histopath_anno_add) <- rownames(histopath_anno)
se <- AddMetaData(se, metadata = histopath_anno_add, col.name = "annotation")

# Add slightly refined annotation column
se$annotation2 <- ifelse(se$annotation == "Inflammation", paste0(se$annotation, "_", se$day), se$annotation)

##### Add pct feature set ####
se <- PercentageFeatureSet(se, pattern = "^mt-", col.name = "percent.mt")
se <- PercentageFeatureSet(se, pattern = "^Hb.*-", col.name = "percent.hb")
se <- PercentageFeatureSet(se, pattern = "^Rps|^Rpl", col.name = "percent.rp")


#### QC ####
##### Plot QC 1 raw ####
p_out <- plotVisiumQC(seurat.object = se, extra.feat.plot = c("percent.mt", "percent.hb"))
pdf(file = file.path(DIR_RES, "figures", "mm_visium_preproc.qc1.pdf"), 
    width = 10, height = 14);p_out;dev.off()

# VlnPlot(se, 
#         features = c("percent.rp"), 
#         pt.size=0, 
#         group.by = "sample_name", 
#         split.by = "group", 
#         cols = cols_group) & theme_custom


##### Filter data ####
#' spot filter
msg1.1 <- paste0("removing ", nrow(subset(se@meta.data, percent.mt >= 30)), " spots due to mitochondial content")
msg1.2 <- paste0("removing ", nrow(subset(se@meta.data, percent.hb >= 30)), " spots due to mitochondial content")
msg1.3 <- paste0("removing ", nrow(subset(se@meta.data, nCount_RNA <= 300)), " spots due to too few unique transcripts") #???!
msg1.4 <- paste0("removing ", nrow(qc_spot_rm), " spots due to manual spot QC selection")

spots_keep <- rownames(subset(se@meta.data, 
                              percent.mt < 30 & percent.hb < 30 & nCount_RNA >300))
spots_keep <- spots_keep[!spots_keep %in% rownames(qc_spot_rm)]

#' gene filter
rp_genes <- grep(pattern = "^Rpl|^Rps", x = rownames(se), value = TRUE)
mt_genes <- grep(pattern = "^mt-", x = rownames(se), value = TRUE)
gene_bt_keep <- c("protein_coding", grep(pattern = "_gene", x = unique(gene_anno$gene_biotype), value = T))
nc_genes <- unique(gene_anno[!gene_anno$gene_biotype %in% gene_bt_keep, "mgi_symbol"])
genes_keep <- rownames(se)[!rownames(se) %in% c(rp_genes, mt_genes, nc_genes)]

msg2.1 <- paste0("removing ", length(rownames(se))-length(genes_keep), " genes")

cat(paste0("\n", msg1.1, "\n", msg1.2,"\n", msg1.3,"\n", msg1.4, "\n", msg2.1), file = analysis_txt_fname, sep="\n", append=TRUE)

#' subset data
se.subset <- SubsetSTData(se, spots = spots_keep, features = genes_keep)
cat(paste0("new dims: ", dim(se.subset)[1], " genes, ", dim(se.subset)[2], " spots"), file = analysis_txt_fname, sep="\n", append=TRUE)
rm(se)

##### Plot QC 1 filtered ####
p_out <- plotVisiumQC(seurat.object = se.subset, extra.feat.plot = c("percent.mt", "percent.hb"))
pdf(file = file.path(DIR_RES, "figures", "mm_visium_preproc.qc1_filtered.pdf"), 
    width = 10, height = 14);p_out;dev.off()

##### Plot QC 2 filtered ####
p1 <- ST.FeaturePlot(se.subset, 
               features = c("nCount_RNA", "nFeature_RNA"), 
               palette = "Spectral", 
               label.by = "sample_name", 
               ncol = 6, pt.size = 0.6)
p2 <- ST.FeaturePlot(se.subset, 
                     features = c("percent.mt", "percent.hb"), 
                     palette = "Spectral", 
                     label.by = "sample_name", 
                     ncol = 6, pt.size = 0.6)
png(file = file.path(DIR_RES, "figures", "mm_visium_preproc.qc2_filtered.png"), 
    width = 16*fig_res, height = 12*fig_res, res = fig_res);p1/p2;dev.off()

##### Plot QC 3 filtered ####
p_dat <- se.subset@meta.data

p_dat1 <- p_dat %>%
  group_by(sample_name, group) %>%
  summarise(n_spots = n())

p_dat2 <- p_dat %>%
  group_by(sample_name, annotation, group) %>%
  summarise(n_spots = n())

p1 <- ggplot(p_dat1, aes(x=sample_name, y=n_spots, fill=group)) +
  geom_col() +
  scale_fill_manual(values = cols_group) +
  coord_flip() +
  labs(y="n spots", x="", title="Spots inclued per sample") +
  theme_bw() +
  theme(legend.position = "top", 
        axis.text = element_text(colour = "black"), 
        plot.title = element_text(hjust=0.5, face = "bold"))

p2 <- ggplot(p_dat2, aes(x=annotation, y=n_spots, fill=group)) +
  geom_boxplot() +
  scale_fill_manual(values = cols_group) +
  labs(y="n spots", x="", title="Spots associated to annotated regions", fill="") +
  scale_y_log10()+
  theme_bw() +
  theme(legend.position = "top", 
        axis.text = element_text(colour = "black"), 
        axis.text.x = element_text(angle=45, hjust = 1),
        plot.title = element_text(hjust=0.5, face = "bold"))

p_out <- p1/p2
pdf(file = file.path(DIR_RES, "figures", "mm_visium_preproc.qc3_filtered.pdf"), 
    width = 6, height = 12);p_out;dev.off()


#### Process data ####
##### SCTransform ####
vars_reg <- c("animal")
cat(paste0("\nSCTranform: Default function, regressing out ", paste(vars_reg, collapse = " ")), file = analysis_txt_fname, sep="\n", append=TRUE)
se.subset <- SCTransform(se.subset, vars.to.regress = vars_reg, conserve.memory = T)

##### Dimensionality reduction ####
###### NMF ####
n_factors <- 20
se.subset <- RunNMF(se.subset, nfactors = n_factors)
cat(paste0("\nNMF: Default function, run with ", n_factors, " factors"), 
    file = analysis_txt_fname, sep="\n", append=TRUE)

###### PCA & UMAP ####
se.subset <- RunPCA(se.subset, npcs = 50) %>% 
  RunUMAP(reduction="pca", dims=1:30)

###### Harmony & UMAP ####
# PCA
red_use <- "pca"
dims_use <- 1:20
integrate_vars <- c("condition", "animal")

se.subset <- RunHarmony(se.subset,
                        group.by.vars = integrate_vars,
                        reduction = red_use,
                        assay.use = "SCT",
                        epsilon.cluster=-Inf,  # prevent early stopping
                        epsilon.harmony=-Inf,  # prevent early stopping
                        max.iter.harmony = 21, #39,  # increase from default 10.  But not much happening after ~iter 17, though alternating
                        plot_convergence = T, 
                        verbose = T)
se.subset <- RunUMAP(se.subset, reduction = "harmony", dims = dims_use, reduction.name = "umap.harmony")

cat(paste0("\nHarmony: Run with ", red_use, " as input and integrating data based on ", paste(integrate_vars, collapse = " "),
           "\nUMAP: Run with Harmony vectors ", paste(dims_use, collapse = " ")), 
    file = analysis_txt_fname, sep="\n", append=TRUE)


# NMF -- Harmany did not do much to integrate samples - was already quite integrated
# red_use <- "NMF"
# dims_use <- 1:20
# integrate_vars <- c("condition", "animal")
# 
# se.subset <- RunHarmony(se.subset,
#                         group.by.vars = integrate_vars,
#                         reduction = red_use,
#                         reduction.save = "harmony.nmf",
#                         assay.use = "SCT",
#                         epsilon.cluster=-Inf,  # prevent early stopping
#                         epsilon.harmony=-Inf,  # prevent early stopping
#                         max.iter.harmony = 21, #39,  # increase from default 10.  But not much happening after ~iter 17, though alternating
#                         plot_convergence = T, 
#                         verbose = T)
# se.subset <- RunUMAP(se.subset, reduction = "harmony.nmf", dims = dims_use, reduction.name = "umap.harmony.nmf")
# se.subset <- RunUMAP(se.subset, reduction = "NMF", dims = dims_use, reduction.name = "umap.nmf")

# DimPlot(se.subset, reduction = "umap.harmony.nmf", group.by = "animal") | DimPlot(se.subset, reduction = "umap.nmf", group.by = "animal")
# DimPlot(se.subset, reduction = "umap.harmony.nmf", group.by = "annotation2") | DimPlot(se.subset, reduction = "umap.nmf", group.by = "annotation2")
# DimPlot(se.subset, reduction = "umap.nmf", group.by = "annotation2", split.by = "annotation2", cols = cols_annotation)


###### Clustering ####
res <- 0.4
se.subset <- FindNeighbors(se.subset, reduction = "harmony", dims = dims_use)
se.subset <- FindClusters(se.subset, 
                          resolution = res, 
                          group.singletons=T)

se.subset$seurat_clusters <- as.numeric(se.subset@meta.data[, paste0("SCT_snn_res.", res)])
se.subset$seurat_clusters <- factor(se.subset$seurat_clusters, levels = sort(unique(se.subset$seurat_clusters)))
n_clusters <- length(levels(se.subset$seurat_clusters))
levels(se.subset$seurat_clusters); n_clusters

cat(paste0("\nClustering: ",
           "\nRun FindNeighbors with Harmony vectors ", paste(dims_use, collapse = " "),
           "\nRun FindClusters with resolution set to ", res,
           "\nGenerated ", n_clusters, " clusters"), file = analysis_txt_fname, sep="\n", append=TRUE)

p1 <- DimPlot(se.subset, reduction = "umap", group.by = "group", cols = cols_group) + labs(title="PCA-based")
p2 <- DimPlot(se.subset, reduction = "umap.harmony", group.by = "group", cols = cols_group) + labs(title="PCA-Harmony-based")
p3 <- DimPlot(subset(se.subset, annotation != "NA"), reduction = "umap.harmony", group.by = "annotation")
p4 <- DimPlot(se.subset, reduction = "umap.harmony", group.by = "seurat_clusters");p4
# p4 <- FeaturePlot(se.subset, reduction = "umap.harmony", features = c("nFeature_RNA"), cols = c("grey90", "black"))

p_out <- (p1+p2)/(p3+p4)
png(file = file.path(DIR_RES, "figures", "mm_visium_preproc.umap_filtered.png"), 
    width = 12*fig_res, height = 10*fig_res, res = fig_res);p_out;dev.off()


#### Other ####
##### Summary table #### 
t_dat <- se.subset@meta.data
t_dat <- t_dat %>%
  group_by(sample_name, sample_id, day, condition, group) %>%
  summarise(n_spots = n(),
            avg_genes_spot = mean(nFeature_RNA),
            avg_umis_spot = mean(nCount_RNA),
            avg_genes_sct_spot = mean(nFeature_RNA),
            avg_umis_sct_spot = mean(nCount_SCT),
            avg_pct_mt_spot = mean(percent.mt),
            avg_pct_hb_spot = mean(percent.hb)
  )
t_dat <- as.data.frame(t_dat)

write.csv(t_dat, file.path(DIR_RES, "objects", "mm_visium_preproc_sample_summary.csv"))


#### Wrap up ####
#' Save rds
fname <- paste0("mm_visium_preproc_se_obj.rds")
saveRDS(se.subset, file = file.path(DIR_RES, "objects", fname))

#' Session info
cat("\n", file = analysis_txt_fname, append=TRUE)
cat(capture.output(sessionInfo()), file = analysis_txt_fname, sep="\n", append=TRUE)

# renv::snapshot()
