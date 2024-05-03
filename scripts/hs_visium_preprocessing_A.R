#' [hs_visium_preprocessing_A.R]
#'
#' Preprocessing of spaceranger output from human lung Visium experiments, Workflow A
#' 
#' Processing workflows:
#' - A: Joint HC and IPF integration
#' - B: Individual donor
#'
#' Aug-Sept 2022, L. Franzén [lovisa.franzen@scilifelab.se]

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

##### Other ####
source(file.path(DIR_ROOT, "scripts", "custom_functions.R"))
source(file.path(DIR_ROOT, "scripts", "custom_colors.R"))
theme_custom <- theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))


#' Create a log file to store parameter settings and intermediate information
analysis_txt_fname <- file.path(DIR_RES, paste0("log.hs_visium_preprocessing_A.txt"))
writeLines("Analysis of human lung visium data generated 2021 by script hs_visium_preprocessing_A.R.", analysis_txt_fname)


##### Read tables ####
gene_anno <- read.table(file = file.path(DIR_ROOT, "data/misc", "hs_gene_biomart.2022-02-23.tsv"), 
                        sep = "\t", 
                        header = T, 
                        stringsAsFactors = F)

metadata <- read.table(file.path(DIR_DATA, "hs_visium_metadata.tsv"), sep = "\t", header = T)
metadata$fibrotic_extent_score_by_pathologist_0.3 <- as.character(metadata$fibrotic_extent_score_by_pathologist_0.3)
rownames(metadata) <- metadata$sample_id

histopath_anno <- read.csv(file.path(DIR_DATA, "hs_visium_merged_histo_annotations.csv"))
rownames(histopath_anno) <- histopath_anno$barcode_n


#### Create Seurat/STutility object ####
infoTable <- createInfoTable(data.dir = DIR_DATA, 
                             sample.select.ids = metadata$sample_id, 
                             metadata.df = metadata, 
                             subfolder.name = NA,
                             use.raw = F)

se <- InputFromTable(infotable = infoTable, 
                     minUMICountsPerSpot = 350,
                     minUMICountsPerGene = 100,
                     minGenesPerSpot = 10,
                     minSpotsPerGene = 5,
                     platform =  "Visium",
                     disable.subset = TRUE)

cat(paste0("Original dims: ", dim(se)[1], " genes, ", dim(se)[2], " spots"), file = analysis_txt_fname, sep="\n", append=TRUE)

#' QC plot
VlnPlot(se, features=c("nCount_RNA", "nFeature_RNA"), group.by="sample_name", pt.size=0, split.by = "subject_alias", cols = cols_donor) & geom_hline(yintercept = 1e3)
ST.FeaturePlot(se, features="nFeature_RNA", label.by="sample_name", pt.size=1, ncol=5)


##### Histopath annotations #### 
#' Add annotations by histopathologist to metadata
histopath_anno_add <- histopath_anno['annotation']
se <- AddMetaData(se, metadata = histopath_anno_add, col.name = "annotation")


##### Add pct feature set ####
se <- PercentageFeatureSet(se, pattern = "^MT-|^MTRNR", col.name = "percent.mt")
se <- PercentageFeatureSet(se, pattern = "^HB", col.name = "percent.hb")
se <- PercentageFeatureSet(se, pattern = "^RPS|^RPL", col.name = "percent.rp")

chry_genes <- intersect(subset(gene_anno, chromosome_name == "Y")$hgnc_symbol, rownames(se))
se <- PercentageFeatureSet(se, features = chry_genes, col.name = "percent.chrY")


#### QC ####
##### Plot QC 1 raw ####
p_out <- plotVisiumQC(seurat.object = se, extra.feat.plot = c("percent.mt", "percent.hb"), fill.by = "subject_alias", color.group = cols_donor)
p_out

##### Filter data ####
#' spot filter
msg1.1 <- paste0("removing ", nrow(subset(se@meta.data, percent.mt >= 30)), " spots due to mitochondial content")
msg1.2 <- paste0("removing ", nrow(subset(se@meta.data, percent.hb >= 30)), " spots due to hemoglobin content")
msg1.3 <- paste0("removing ", nrow(subset(se@meta.data, nCount_RNA <= 350)), " spots due to too few unique transcripts")

spots_keep <- rownames(subset(se@meta.data, percent.mt < 30 & percent.hb < 30 & nCount_RNA >350))  # nCount_RNA filter also done with InputFromTable

#' gene filter
rp_genes <- grep(pattern = "^RPS|^RPL", x = rownames(se), value = TRUE)
mt_genes <- grep(pattern = "^MT-|^MTRNR", x = rownames(se), value = TRUE)
xychr_genes <- unique(subset(gene_anno, chromosome_name %in% c("X", "Y"))$hgnc_symbol)
gene_bt_keep <- c("protein_coding", grep(pattern = "_gene", x = unique(gene_anno$gene_biotype), value = T))
nc_genes <- unique(c(gene_anno[!gene_anno$gene_biotype %in% gene_bt_keep, "hgnc_symbol"],
                   grep("^LINC[0-9][0-9][0-9]", rownames(se), value = T),
                   grep("^A[A-Z][0-9][0-9][0-9]", rownames(se), value = T))
                   )

genes_keep <- rownames(se)[!rownames(se) %in% c(rp_genes, mt_genes, xychr_genes, nc_genes)]

msg2.1 <- paste0("removing ", length(rownames(se))-length(genes_keep), " genes")

cat(paste0("\n", msg1.1, "\n", msg1.2,"\n", msg1.3, "\n", msg2.1), file = analysis_txt_fname, sep="\n", append=TRUE)

#' subset data
se.subset <- SubsetSTData(se, spots = spots_keep, features = genes_keep)
cat(paste0("new dims: ", dim(se.subset)[1], " genes, ", dim(se.subset)[2], " spots"), file = analysis_txt_fname, sep="\n", append=TRUE)


##### Plot QC 1 filtered ####
p_out <- plotVisiumQC(seurat.object = se.subset, extra.feat.plot = c("percent.mt", "percent.hb"), fill.by = "subject_alias", color.group = cols_donor)
pdf(file = file.path(DIR_RES, "figures", "hs_visium_preproc_A.qc1_filtered.pdf"), 
    width = 10, height = 14);p_out;dev.off()


##### Plot QC 2 filtered ####
p1 <- ST.FeaturePlot(se.subset, 
               features = c("nCount_RNA", "nFeature_RNA"), 
               cols = col_scale_rocket,
               label.by = "sample_name", 
               ncol = 5, pt.size = 0.6)
p2 <- ST.FeaturePlot(se.subset, 
                     features = c("percent.mt", "percent.hb"), 
                     cols = col_scale_rocket,
                     label.by = "sample_name", 
                     ncol = 5, pt.size = 0.6)
png(file = file.path(DIR_RES, "figures", "hs_visium_preproc_A.qc2_filtered.png"), 
    width = 16*fig_res, height = 16*fig_res, res = fig_res);p1/p2;dev.off()


##### Plot QC 3 filtered ####
p_dat <- se.subset@meta.data

p_dat1 <- p_dat %>%
  group_by(sample_name, condition) %>%
  summarise(n_spots = n()) %>%
  as.data.frame()

p1 <- ggplot(p_dat1, aes(x=sample_name, y=n_spots, fill=condition)) +
  geom_col() +
  scale_fill_manual(values = cols_cond) +
  coord_flip() +
  labs(y="n spots", x="", title="Spots included per sample") +
  theme_bw() +
  theme(legend.position = "top",
        axis.text = element_text(colour = "black"),
        plot.title = element_text(hjust=0.5, face = "bold"))

p2 <- ggplot(p_dat, aes(x=nCount_RNA, y=nFeature_RNA, color=paste0(B_tissue_selection, "-R", replicate))) +
  geom_point(size=0.2) +
  facet_wrap(~subject_alias, ncol = 4) +
  labs(color="Sample") +
  scale_color_manual(values = RColorBrewer::brewer.pal(8,"Set2")) +
  theme_linedraw() +
  theme(legend.position = "bottom", panel.grid = element_blank(), axis.text.x = element_text(angle = 90, vjust = .5)) +
  guides(colour = guide_legend(override.aes = list(size=2)))


p_out <- p1/p2 +patchwork::plot_layout(heights = c(3,2))
pdf(file = file.path(DIR_RES, "figures", "hs_visium_preproc_A.qc3_filtered.pdf"),
    width = 6, height = 12);p_out;dev.off()


#### Process data ####
##### SCTransform ####
vars_reg <- c("sample_id", "subject_alias")
cat(paste0("\nSCTranform: Default function, regressing out ", paste(vars_reg, collapse = " ")), file = analysis_txt_fname, sep="\n", append=TRUE)
se.subset <- SCTransform(se.subset, vars.to.regress = vars_reg, conserve.memory = T)

p_dat <- se.subset@meta.data
p1 <- ggplot(p_dat, aes_string(x="nCount_SCT")) +
  geom_histogram(aes(y=..density..), colour=NA, fill="black", bins = 75)+
  geom_density(data = p_dat, aes_string(x="nCount_SCT", y="..density..", fill = "subject_alias"), alpha=.5) +
  geom_vline(aes(xintercept=mean(nCount_SCT)), color="orange", linetype="dashed", size=1) +
  scale_fill_manual(values = cols_donor) +
  ggtitle("Normalized and corrected UMIs per spot") +
  theme_classic() + 
  theme(legend.position = "bottom", 
        axis.text = element_text(colour = "black"), 
        plot.title = element_text(hjust=0.5, face = "bold"))
p2 <- VlnPlot(se.subset, 
              features = c("nFeature_SCT", "nCount_SCT"), 
              pt.size=0,
              group.by = "sample_name", 
              split.by = "subject_alias", 
              cols = cols_donor) & theme_custom
p3 <- VlnPlot(se.subset, 
        assay = "SCT", 
        features = c("GAPDH", "ACTB"),
        pt.size=0,
        group.by = "sample_name", 
        split.by = "subject_alias", 
        cols = cols_donor) & theme_custom

p_out <- p1/p2/p3 #+patchwork::plot_layout(heights = c(1,2))
pdf(file = file.path(DIR_RES, "figures", "hs_visium_preproc_A.qc1_sct.pdf"), 
    width = 10, height = 14);p_out;dev.off()


##### Dimensionality reduction ####
###### NMF ####
n_factors <- 30
se.subset <- RunNMF(se.subset, nfactors = n_factors)
cat(paste0("\nNMF: Default function, run with ", n_factors, " factors"), 
    file = analysis_txt_fname, sep="\n", append=TRUE)


#### Other ####
##### Summary table #### 
t_dat <- se.subset@meta.data
t_dat <- t_dat %>%
  group_by(sample_name, sample_id, condition, subject_alias, fibrotic_extent_score_by_pathologist_0.3) %>%
  summarise(n_spots = n(),
            avg_genes_spot = mean(nFeature_RNA),
            avg_umis_spot = mean(nCount_RNA),
            avg_genes_sct_spot = mean(nFeature_RNA),
            avg_umis_sct_spot = mean(nCount_SCT),
            avg_pct_mt_spot = mean(percent.mt),
            avg_pct_hb_spot = mean(percent.hb)
  ) %>%
  as.data.frame()

write.csv(t_dat, file.path(DIR_RES, "objects", "hs_visium_preproc_A_sample_summary.csv"), row.names = F)


#### Wrap up ####
#' Save rds
fname <- paste0("hs_visium_preproc_A_se_obj_nmf.rds")
saveRDS(se.subset, file = file.path(DIR_RES, "objects", fname))

#' Session info
cat("\n", file = analysis_txt_fname, append=TRUE)
cat(capture.output(sessionInfo()), file = analysis_txt_fname, sep="\n", append=TRUE)

# renv::snapshot()
