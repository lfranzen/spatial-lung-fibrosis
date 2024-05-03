#' [hs_mm_visium_comparison_DEA.R]
#'
#' Comparison of pseudo-bulk regions by DEA in IPF vs BLM Visium data 
#'
#' L. Franzén [lovisa.franzen@scilifelab.se]

#### Set up ####
##### Define params. ####
SPECIES <- "translational"
DIR_ROOT <- getwd()
DIR_DATA <- file.path(DIR_ROOT, "data")
DIR_RES <- file.path(DIR_ROOT, "results", SPECIES)
DIR_FIG <- file.path(DIR_RES, "figures")
DIR_FIG_OUT <- file.path(DIR_FIG, "DEA")
DIR_OBJ_OUT <- file.path(DIR_RES, "objects", "DEA")
dir.create(DIR_FIG_OUT); dir.create(DIR_OBJ_OUT)

##### Load libs ####
# library(tidyverse)
library(dplyr)
library(STutility)
library(readxl)
library(writexl)
library(DESeq2)
library(patchwork)
library(gprofiler2)
library(ggVennDiagram)

##### Other ####
source(file.path(DIR_ROOT, "scripts", "custom_functions.R"))
source(file.path(DIR_ROOT, "scripts", "custom_colors.R"))
theme_custom <- theme(axis.title.x = element_blank())
fig_res <- 300

#### Gene conversion data ####
gene_conv_df <- read.csv(file.path(DIR_DATA, "misc", "orthogene_conv_combined_hs_mm.csv"), row.names = 1)
rownames(gene_conv_df) <- gene_conv_df$symbol_hs_mm

gene_conv <- setNames(gene_conv_df$symbol_mm, nm = gene_conv_df$symbol_hs)
gene_conv_mm <- setNames(gene_conv_df$symbol_hs, nm = gene_conv_df$symbol_mm)


#### Read pseudobulk data ####
SPECIES <- "mouse"
list.files(file.path(DIR_RES, "..", SPECIES, "objects", "DEA"))

SPECIES <- "human"
list.files(file.path(DIR_RES, "..", SPECIES, "objects", "DEA"))

##### Metadata ##### 
# Read files
SPECIES <- "mouse"
mdat_mm <- read.csv(file.path(DIR_RES, "..", SPECIES, "objects", "DEA", 
                              "mm_visium_pseudobulk_metadata_per_animal.csv"), 
                    row.names = 1)
mdat_mm <- mdat_mm[, -2]
mdat_mm$condition <- factor(mdat_mm$condition, levels = c("control","bleomycin"))

SPECIES <- "human"
mdat_hs <- read.csv(file.path(DIR_RES, "..", SPECIES, "objects", "DEA", 
                             "hs_visium_pseudobulk_metadata_per_donor.csv"), 
                   row.names = 1)
mdat_hs$condition <- factor(mdat_hs$condition, levels = c("control","IPF"))

# Format files to have same colnames
colnames(mdat_mm) <- colnames(mdat_hs) <- c("donor", "condition", "rep")


##### All ##### 
# mouse
SPECIES <- "mouse"
dat_mm_all <- read.csv(file.path(DIR_RES, "..", SPECIES, "objects", "DEA", 
                                 "mm_visium_pseudobulk_data_per_animal_all.csv"), 
                       row.names = 1)

# human
SPECIES <- "human"
dat_hs_all <- read.csv(file.path(DIR_RES, "..", SPECIES, "objects", "DEA", 
                                 "hs_visium_pseudobulk_data_grouped_by_donor_all.csv"), 
                       row.names = 1)

##### Alveolar ##### 
# mouse
SPECIES <- "mouse"
dat_mm_alv <- read.csv(file.path(DIR_RES, "..", SPECIES, "objects", "DEA", 
                                 "mm_visium_pseudobulk_data_per_animal_NormalAlveolarandOther.csv"), 
                       row.names = 1)
# dat_mm_alv <- dat_mm_alv[, colnames(dat_mm_alv) %>% grep(pattern = "d21_b", value = T)]

# human
SPECIES <- "human"
dat_hs_alv <- read.csv(file.path(DIR_RES, "..", SPECIES, "objects", "DEA", 
                                  "hs_visium_pseudobulk_data_per_donor_NormalAlveolarandOther.csv"), 
                        row.names = 1)
dat_hs_alv <- dat_hs_alv[, colnames(dat_hs_alv) %>% grep(pattern = "IPF", value = T)]

dat_hs_alv_all <- read.csv(file.path(DIR_RES, "..", SPECIES, "objects", "DEA", 
                                     "hs_visium_pseudobulk_data_per_donor_NormalAlveolarandOther.csv"), 
                           row.names = 1)

##### Fibrosis ##### 
# mouse
SPECIES <- "mouse"
dat_mm_fib <- read.csv(file.path(DIR_RES, "..", SPECIES, "objects", "DEA", 
                                 "mm_visium_pseudobulk_data_per_animal_SuspectFibrosisFibroplasia.csv"), 
                       row.names = 1)
# dat_mm_fib <- dat_mm_fib[, colnames(dat_mm_fib) %>% grep(pattern = "d21_b", value = T)]

dat_mm_fib_d7 <- read.csv(file.path(DIR_RES, "..", SPECIES, "objects", "DEA", 
                                    "mm_visium_pseudobulk_data_per_animal_Inflammation.csv"), 
                        row.names = 1) %>% select(contains("d7_"))

dat_mm_fib <- bind_cols(dat_mm_fib, dat_mm_fib_d7)

# human
SPECIES <- "human"
dat_hs_fib1 <- read.csv(file.path(DIR_RES, "..", SPECIES, "objects", "DEA", 
                                 "hs_visium_pseudobulk_data_per_donor_SuspectFibrosis.csv"), 
                       row.names = 1)
dat_hs_fib2 <- read.csv(file.path(DIR_RES, "..", SPECIES, "objects", "DEA", 
                                  "hs_visium_pseudobulk_data_per_donor_Diseased.csv"), 
                        row.names = 1)

shared_columns <- intersect(colnames(dat_hs_fib1), colnames(dat_hs_fib2))
shared_rows <- intersect(rownames(dat_hs_fib1), rownames(dat_hs_fib2))

dat_hs_fib_list <- lapply(shared_columns, function(c){
  summed_counts <- rowSums(cbind(dat_hs_fib1[shared_rows, c], 
                                 dat_hs_fib2[shared_rows, c]))
  names(summed_counts) <- shared_rows
  return(summed_counts)
}) %>% setNames(nm = shared_columns)

dat_hs_fib <- bind_cols(dat_hs_fib_list) %>% as.data.frame()
rownames(dat_hs_fib) <- shared_rows


##### Inflammation ##### 
# mouse
SPECIES <- "mouse"
dat_mm_infl <- read.csv(file.path(DIR_RES, "..", SPECIES, "objects", "DEA", 
                                 "mm_visium_pseudobulk_data_per_animal_Inflammation.csv"), 
                       row.names = 1)
dat_mm_infl <- dat_mm_infl[, colnames(dat_mm_infl) %>% grep(pattern = "d21_b", value = T)]


# human
SPECIES <- "human"
dat_hs_infl <- read.csv(file.path(DIR_RES, "..", SPECIES, "objects", "DEA", 
                                 "hs_visium_pseudobulk_data_per_donor_Inflammation.csv"), 
                       row.names = 1)
dat_hs_infl <- dat_hs_infl[, colnames(dat_hs_infl) %>% grep(pattern = "IPF", value = T)]



##### AbBa NMF factor hi subcluster ##### 
#' Pseudo-bulk of NMF30-F14-C0 for hs/mm
SPECIES <- "mouse"
dat_mm_f14c0 <- read.csv(file.path(DIR_RES, "..", SPECIES, "objects", "DEA", 
                                   "mm_visium_pseudobulk_data_per_animal_d21-NMF30-F14hi-C0.csv"), 
                         row.names = 1)

# human
SPECIES <- "human"
dat_hs_f14c0 <- read.csv(file.path(DIR_RES, "..", SPECIES, "objects", "DEA", 
                                   "hs_visium_pseudobulk_data_per_donor_NMF30-F14hi-C0.csv"), 
                        row.names = 1)


#### Run DEA ####
# Function for performing and extracting DEA results per region
RunRegionSpeciesDEA <- function(
    bulk_data_mm, 
    bulk_data_hs, 
    mdata_mm,
    mdata_hs,
    gene_conv_df){
  
  # Merge data sets
  # identify by overlapping genes
  genes_mm <- rownames(bulk_data_mm)[rownames(bulk_data_mm) %in% gene_conv_df$symbol_mm]
  genes_hs <- rownames(bulk_data_hs)[rownames(bulk_data_hs) %in% gene_conv_df$symbol_hs]
  genes_shared <- intersect(gene_conv_df[gene_conv_df$symbol_mm %in% genes_mm,]$symbol_hs_mm,
                            gene_conv_df[gene_conv_df$symbol_hs %in% genes_hs,]$symbol_hs_mm)
  # subset data to shared genes
  bulk_data_mm_subset <- bulk_data_mm[gene_conv_df[gene_conv_df$symbol_hs_mm %in% genes_shared,]$symbol_mm, ]
  bulk_data_hs_subset <- bulk_data_hs[gene_conv_df[gene_conv_df$symbol_hs_mm %in% genes_shared,]$symbol_hs, ]
  rownames(bulk_data_mm_subset) <- rownames(bulk_data_hs_subset) <- genes_shared
  
  # join data
  bulk_data <- cbind(bulk_data_mm_subset, bulk_data_hs_subset)
  message(paste0("Combined dataset with dim: rows ", paste(dim(bulk_data), collapse = ", cols ")))
  
  # Merge metadata
  mdata <- rbind(mdat_mm[colnames(bulk_data_mm_subset),],
                 mdat_hs[colnames(bulk_data_hs_subset),])
  
  # Run DEA
  message("Running DEA based on 'condition'")
  dds <- DESeqDataSetFromMatrix(countData = bulk_data,
                                colData = mdata,
                                design = ~ condition)
  dds <- DESeq(dds)
  
  res_name <- resultsNames(dds)[2]
  message(paste0("Fetching results for comparison '", res_name, "'"))
  res_dds <- results(dds, name=res_name)
  
  message("Returning results table")
  res_df <- res_dds %>% 
    as.data.frame()
  
  # add gene name columns
  res_df$gene <- rownames(res_df)
  res_df$gene_hs <- gene_conv_df[gene_conv_df$symbol_hs_mm %in% rownames(res_df), "symbol_hs"]
  res_df$gene_mm <- gene_conv_df[gene_conv_df$symbol_hs_mm %in% rownames(res_df), "symbol_mm"]
  
  res_df$design <- res_name
  
  return(res_df)
}

# res_df_test <- RunRegionSpeciesDEA(bulk_data_mm = dat_mm_fib,
#                                    bulk_data_hs = dat_hs_fib,
#                                    mdata_mm = mdat_mm,
#                                    mdata_hs = mdat_hs, 
#                                    gene_conv_df = gene_conv_df)


##### All ##### 
###### DEA mouse ######
#' d7
dat_mm_all_d7 <- dat_mm_all[,grep("d7", colnames(dat_mm_all), value = T)]
dds_all_mm_d7 <- DESeqDataSetFromMatrix(countData = dat_mm_all_d7,
                                     colData = mdat_mm[colnames(dat_mm_all_d7),],
                                     design = ~ condition)
dds_all_mm_d7 <- DESeq(dds_all_mm_d7)

res_name <- resultsNames(dds_all_mm_d7)[2];res_name
res_all_mm_d7 <- results(dds_all_mm_d7, name=res_name)
res_all_mm_d7$design <- paste0(res_name, ".d7")
res_all_mm_d7$species <- "mouse"

#' d21
dat_mm_all_d21 <- dat_mm_all[,grep("d21", colnames(dat_mm_all), value = T)]
dds_all_mm <- DESeqDataSetFromMatrix(countData = dat_mm_all_d21,
                                     colData = mdat_mm[colnames(dat_mm_all_d21),],
                                     design = ~ condition)
dds_all_mm <- DESeq(dds_all_mm)

res_name <- resultsNames(dds_all_mm)[2];res_name
res_all_mm <- results(dds_all_mm, name=res_name)
res_all_mm$design <- paste0(res_name, ".d21")
res_all_mm$species <- "mouse"

###### DEA human ######
dds_all_hs <- DESeqDataSetFromMatrix(countData = dat_hs_all,
                                     colData = mdat_hs[colnames(dat_hs_all),],
                                     design = ~ condition)
dds_all_hs <- DESeq(dds_all_hs)

res_name <- resultsNames(dds_all_hs)[2];res_name
res_all_hs <- results(dds_all_hs, name=res_name)
res_all_hs$design <- res_name
res_all_hs$species <- "human"

###### Join data ######
# identify by overlapping genes
genes_mm <- rownames(res_all_mm)[rownames(res_all_mm) %in% gene_conv_df$symbol_mm]
genes_mm <- c(genes_mm, rownames(res_all_mm_d7)[rownames(res_all_mm_d7) %in% gene_conv_df$symbol_mm]) %>% unique()
genes_hs <- rownames(res_all_hs)[rownames(res_all_hs) %in% gene_conv_df$symbol_hs]
genes_shared <- intersect(gene_conv_df[gene_conv_df$symbol_mm %in% genes_mm,]$symbol_hs_mm,
                          gene_conv_df[gene_conv_df$symbol_hs %in% genes_hs,]$symbol_hs_mm)
length(genes_shared)

# subset data to shared genes
res_all_mm_subset <- res_all_mm[gene_conv_df[gene_conv_df$symbol_hs_mm %in% genes_shared,]$symbol_mm, ]
res_all_mm_subset_d7 <- res_all_mm_d7[gene_conv_df[gene_conv_df$symbol_hs_mm %in% genes_shared,]$symbol_mm, ]
res_all_hs_subset <- res_all_hs[gene_conv_df[gene_conv_df$symbol_hs_mm %in% genes_shared,]$symbol_hs, ]
rownames(res_all_mm_subset) <- rownames(res_all_hs_subset) <- rownames(res_all_mm_subset_d7) <- genes_shared

res_all_mm_subset$gene <- res_all_mm_subset_d7$gene <- res_all_hs_subset$gene <- genes_shared
res_all_mm_subset$gene_hs <- res_all_mm_subset_d7$gene_hs <- res_all_hs_subset$gene_hs <- gene_conv_df[gene_conv_df$symbol_hs_mm %in% genes_shared, "symbol_hs"]
res_all_mm_subset$gene_mm <- res_all_mm_subset_d7$gene_mm <- res_all_hs_subset$gene_mm <- gene_conv_df[gene_conv_df$symbol_hs_mm %in% genes_shared, "symbol_mm"]

# Join data
res_all <- rbind(res_all_mm_subset, res_all_mm_subset_d7)
res_all <- rbind(res_all, res_all_hs_subset)

# save
write.csv(res_all, file.path(DIR_OBJ_OUT, "hs_ms_visium_pseudobulk_dea_res_joined_Control_IPF_or_BLMd21_or_BLMd7_region_all.csv"), row.names = T)
res_all <- read.csv(file.path(DIR_OBJ_OUT, "hs_ms_visium_pseudobulk_dea_res_joined_Control_IPF_or_BLMd21_or_BLMd7_region_all.csv"), row.names = 1)

# format
res_all_df <- res_all %>% 
  as.data.frame() %>% 
  dplyr::arrange(padj, abs(log2FoldChange))


###### Compare results ###### 
# Look at sign data
res_all_sign <- res_all_df %>% 
  dplyr::arrange(padj, abs(log2FoldChange)) %>%
  dplyr::filter(padj<0.01 & abs(log2FoldChange) > 1) # FC>2
dim(res_all_sign)


# IPF, d7, d21
unique_mm_d7 <- setdiff(subset(res_all_sign, species == "mouse" & design == "condition_bleomycin_vs_control.d7")$gene,
                        c(subset(res_all_sign, species == "mouse" & design == "condition_bleomycin_vs_control.d21")$gene,
                          subset(res_all_sign, species == "human")$gene))
unique_mm_d21 <- setdiff(subset(res_all_sign, species == "mouse" & design == "condition_bleomycin_vs_control.d21")$gene,
                        c(subset(res_all_sign, species == "mouse" & design == "condition_bleomycin_vs_control.d7")$gene,
                          subset(res_all_sign, species == "human")$gene))
unique_hs <- setdiff(subset(res_all_sign, species == "human")$gene,
                     c(subset(res_all_sign, species == "mouse" & design == "condition_bleomycin_vs_control.d21")$gene,
                       subset(res_all_sign, species == "mouse" & design == "condition_bleomycin_vs_control.d7")$gene))
shared_mm_hs <- intersect(
  intersect(
    subset(res_all_sign, species == "mouse" & design == "condition_bleomycin_vs_control.d7")$gene,
    subset(res_all_sign, species == "mouse" & design == "condition_bleomycin_vs_control.d21")$gene),
  subset(res_all_sign, species == "human")$gene)
length(unique_mm_d7);length(unique_mm_d21);length(unique_hs);length(shared_mm_hs)


# Data frame summarizing DEGs (IPF vs BLM d7 / IPF vs BLM d21)
dea_comp_df_all_groups <- data.frame(area = c("all", "fibrosis", "alveolar"),
                                     DEGs_hs = rep(0,3),
                                     DEGs_mmd7 = rep(0,3),
                                     DEGs_mmd21 = rep(0,3),
                                     DEGs_shared_all = rep(0,3),
                                     DEGs_shared_hs_mmd7 = rep(0,3),
                                     DEGs_shared_hs_mmd21 = rep(0,3),
                                     DEGs_shared_mmd7_mmd21 = rep(0,3),
                                     shared_up = rep(0,3),
                                     shared_down = rep(0,3))

res_all_shared_df <- res_all_sign %>% 
  filter(gene %in% shared_mm_hs) %>% 
  arrange(gene) %>% 
  select(gene, design, log2FoldChange) %>% 
  pivot_wider(names_from = design, values_from = log2FoldChange)

colnames(res_all_shared_df) <- c("gene", "mm_d7", "mm_d21", "hs")

res_all_shared_df <- res_all_shared_df %>% 
  mutate(up_mmd7_hs = ifelse(mm_d7>0 & hs>0, 1, 0),
         down_mmd7_hs = ifelse(mm_d7<0 & hs<0, 1, 0),
         diff_mmd7_hs = ifelse((mm_d7>0 & hs<0) | (mm_d7<0 & hs>0), 1, 0),
         up_mmd21_hs = ifelse(mm_d21>0 & hs>0, 1, 0),
         down_mmd21_hs = ifelse(mm_d21<0 & hs<0, 1, 0),
         diff_mmd21_hs = ifelse((mm_d21>0 & hs<0) | (mm_d21<0 & hs>0), 1, 0),
         up_mmd7_mmd21 = ifelse(mm_d7>0 & mm_d21>0, 1, 0),
         down_mmd7_mmd21 = ifelse(mm_d7<0 & mm_d21<0, 1, 0),
         diff_mmd7_mmd21 = ifelse((mm_d7>0 & mm_d21<0) | (mm_d7<0 & mm_d21>0), 1, 0),
         up_all = ifelse(mm_d7>0 & mm_d21>0 & hs>0, 1, 0),
         down_all = ifelse(mm_d7<0 & mm_d21<0 & hs<0, 1, 0)) 
  # %>% 
  # replace(is.na(.), 0)

dea_comp_df_all_groups[dea_comp_df_all_groups$area == "all", -1] <- c(
  unique_hs %>% length(),          # DEGs_hs
  unique_mm_d7 %>% length(),       # DEGs_mmd7
  unique_mm_d21 %>% length(),      # DEGs_mmd21
  shared_mm_hs %>% length(),       # DEGs_shared_all
  sum(res_all_shared_df$up_mmd7_hs)+sum(res_all_shared_df$down_mmd7_hs),   # DEGs_shared_hs_mmd7
  sum(res_all_shared_df$up_mmd21_hs)+sum(res_all_shared_df$down_mmd21_hs), # DEGs_shared_hs_mmd21
  sum(res_all_shared_df$up_mmd7_mmd21)+sum(res_all_shared_df$down_mmd7_mmd21), # DEGs_shared_mmd7_mmd21
  sum(res_all_shared_df$up_all),   # shared_up
  sum(res_all_shared_df$down_all)  # shared_down
)


# IPF, d21
unique_mm_d21 <- setdiff(subset(res_all_sign, species == "mouse" & design == "condition_bleomycin_vs_control.d21")$gene,
                         subset(res_all_sign, species == "human")$gene)
unique_hs <- setdiff(subset(res_all_sign, species == "human")$gene,
                     subset(res_all_sign, species == "mouse" & design == "condition_bleomycin_vs_control.d21")$gene)

shared_mm_hs <- intersect(subset(res_all_sign, species == "mouse" & design == "condition_bleomycin_vs_control.d21")$gene,
                          subset(res_all_sign, species == "human")$gene)

length(unique_mm_d21);length(unique_hs);length(shared_mm_hs)

share_df <- data.frame(shared_gene = shared_mm_hs,
                       shared_gene_hs = gene_conv_df[gene_conv_df$symbol_hs_mm %in% shared_mm_hs,][shared_mm_hs,]$symbol_hs,
                       shared_gene_mm = gene_conv_df[gene_conv_df$symbol_hs_mm %in% shared_mm_hs,][shared_mm_hs,]$symbol_mm)
share_df <- share_df %>% arrange(shared_gene)
View(share_df)

# Data frame summarizing DEGs (IPF vs BLM d21)
dea_comp_df <- data.frame(area = c("all", "fibrosis", "alveolar"),
                          DEGs_hs = rep(0,3),
                          DEGs_mm = rep(0,3),
                          DEGs_shared = rep(0,3),
                          shared_up = rep(0,3),
                          shared_down = rep(0,3),
                          shared_diff = rep(0,3))

res_all_shared_df <- subset(res_all_sign, design != "condition_bleomycin_vs_control.d7") %>% 
  filter(gene %in% share_df$shared_gene) %>% 
  arrange(gene) %>% 
  select(gene, species, log2FoldChange) %>% 
  pivot_wider(names_from = species, values_from = log2FoldChange) %>% 
  mutate(up = ifelse(mouse>0 & human>0, 1, 0),
         down = ifelse(mouse<0 & human<0, 1, 0),
         diff = ifelse((mouse>0 & human<0) | (mouse<0 & human>0), 1, 0))

dea_comp_df[dea_comp_df$area == "all", -1] <- c(
  unique_hs %>% length(),
  unique_mm_d21 %>% length(),
  shared_mm_hs %>% length(),
  sum(res_all_shared_df$up),
  sum(res_all_shared_df$down),
  sum(res_all_shared_df$diff)
)

# library(ggVennDiagram)
# Venn: all three groups
d_venn <- list(
  BLM.d7 = subset(res_all_sign, species == "mouse" & design == "condition_bleomycin_vs_control.d7")$gene,
  IPF = subset(res_all_sign, species == "human")$gene,
  BLM.d21 = subset(res_all_sign, species == "mouse" & design == "condition_bleomycin_vs_control.d21")$gene
  )

# Save data
venn_df_save <- data.frame(group = c(rep("IPF", length(d_venn$IPF)), rep("BLM.d7", length(d_venn$BLM.d7)), rep("BLM.d21", length(d_venn$BLM.d21))),
                           gene = c(d_venn$IPF, d_venn$BLM.d7, d_venn$BLM.d21))
venn_df_save$overlap <- ifelse(duplicated(venn_df_save$gene)|duplicated(venn_df_save$gene, fromLast=T),
                               "two_groups", 
                               venn_df_save$group)

write.csv(venn_df_save,
          file.path(DIR_OBJ_OUT, "hs_ms_visium_pseudobulk_dea_res_joined_Control_IPF_or_BLMd21_oe_BLMd7_region_all_DEGs_venn.csv"), row.names = T)

# png(filename = file.path(DIR_FIG_OUT, "hs_ms_visium_pseudobulk_dea_res_joined_Control_IPF_or_BLMd21_or_BLMd7_region_all_DEGs_venn.png"),
#     width = 3*fig_res, height = 3*fig_res, res = fig_res)
pdf(file.path(DIR_FIG_OUT, "hs_ms_visium_pseudobulk_dea_res_joined_Control_IPF_or_BLMd21_or_BLMd7_region_all_DEGs_venn.pdf"), 
    width = 5, height = 3)
ggVennDiagram(d_venn, edge_color = "black", edge_lty = "dashed") +
  scale_color_manual(values = c("#215288", "#B26694", "#215288")) +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  theme(legend.position = "none")
dev.off()

# Venn: IPF and BLM d21
d_venn2 <- list(
  IPF = subset(res_all_sign, species == "human")$gene,
  BLM = subset(res_all_sign, species == "mouse" & design == "condition_bleomycin_vs_control.d21")$gene)

# Save data
venn_df_save <- data.frame(group = c(rep("IPF", length(d_venn2$IPF)), rep("BLM", length(d_venn2$BLM))),
                           gene = c(d_venn2$IPF, d_venn2$BLM))
venn_df_save$overlap <- ifelse(duplicated(venn_df_save$gene)|duplicated(venn_df_save$gene, fromLast=T),
                               "both", 
                               venn_df_save$group)
write.csv(venn_df_save,
          file.path(DIR_OBJ_OUT, "hs_ms_visium_pseudobulk_dea_res_joined_Control_IPF_or_BLMd21_region_all_DEGs_venn.csv"), row.names = T)

# Plot
pdf(file.path(DIR_FIG_OUT, "hs_ms_visium_pseudobulk_dea_res_joined_Control_IPF_or_BLMd21_region_all_DEGs_venn.pdf"), 
    width = 1.5, height = 1.5)
ggVennDiagram(d_venn2, label_alpha = 0,
              # edge_color = "black", 
              edge_lty = "solid", 
              label = "count") +
  scale_color_manual(values = c("#B26694", "#215288")) +
  scale_fill_distiller(palette = "Greys", direction = -1) +
  theme(legend.position = "none")
dev.off()


##### Alveolar ##### 
rownames(dat_mm_alv)[rowSums(dat_mm_alv)>10] %>% length()
rownames(dat_hs_alv)[rowSums(dat_hs_alv)>10] %>% length()

###### DEA ms vs hs ######
dat_mm_alv_blm <- dat_mm_alv[, colnames(dat_mm_alv) %>% grep(pattern = "d21_b", value = T)]
res_alv_df <- RunRegionSpeciesDEA(bulk_data_mm = dat_mm_alv_blm,
                                  bulk_data_hs = dat_hs_alv,
                                  mdata_mm = mdat_mm,
                                  mdata_hs = mdat_hs, 
                                  gene_conv_df = gene_conv_df)

write.csv(res_alv_df, file.path(DIR_OBJ_OUT, "hs_ms_visium_pseudobulk_dea_res_IPF_BLMd21_region_alveolar.csv"), row.names = T)
res_alv_df <- read.csv(file.path(DIR_OBJ_OUT, "hs_ms_visium_pseudobulk_dea_res_IPF_BLMd21_region_alveolar.csv"), row.names = 1)

res_alv_sign <- res_alv_df %>% 
  dplyr::arrange(padj, abs(log2FoldChange)) %>%
  dplyr::filter(padj<0.01)

head(res_alv_sign)


###### DEA mouse ######
#' d7
dat_mm_alv_d7_blm <- dat_mm_alv[,grep("d7_b", colnames(dat_mm_alv), value = T)]
dat_mm_all_d7_ctrl <- dat_mm_all[,grep("d7_c", colnames(dat_mm_all), value = T)]
dat_mm_alv_vs_ctrl_d7 <- bind_cols(dat_mm_all_d7_ctrl, dat_mm_alv_d7_blm)

dds_alv_mm <- DESeqDataSetFromMatrix(countData = dat_mm_alv_vs_ctrl_d7,
                                     colData = mdat_mm[colnames(dat_mm_alv_vs_ctrl_d7),],
                                     design = ~ condition)
dds_alv_mm <- DESeq(dds_alv_mm)

res_name <- resultsNames(dds_alv_mm)[2];res_name
res_alv_mm <- results(dds_alv_mm, name=res_name)
res_alv_mm$design <- paste0(res_name, ".alv.d7")
res_alv_mm$species <- "mouse"
res_alv_mm_d7_df <- as.data.frame(res_alv_mm)
res_alv_mm_d7_df$gene_mm <- rownames(res_alv_mm_d7_df)
res_alv_mm_d7_df$gene_hs <- "NA"

write.csv(res_alv_mm_d7_df, file.path(DIR_OBJ_OUT, "mm_visium_pseudobulk_dea_res_BLMd7_region_alveolar_vs_all_CTRLd7.csv"), row.names = T)
res_alv_mm_d7_df <- read.csv(file.path(DIR_OBJ_OUT, "mm_visium_pseudobulk_dea_res_BLMd7_region_alveolar_vs_all_CTRLd7.csv"), row.names = 1)


#' d21
dat_mm_alv_d21_blm <- dat_mm_alv[,grep("d21_b", colnames(dat_mm_alv), value = T)]
dat_mm_all_d21_ctrl <- dat_mm_all[,grep("d21_c", colnames(dat_mm_all), value = T)]
dat_mm_alv_vs_ctrl_d21 <- bind_cols(dat_mm_all_d21_ctrl, dat_mm_alv_d21_blm)

dds_alv_mm <- DESeqDataSetFromMatrix(countData = dat_mm_alv_vs_ctrl_d21,
                                     colData = mdat_mm[colnames(dat_mm_alv_vs_ctrl_d21),],
                                     design = ~ condition)
dds_alv_mm <- DESeq(dds_alv_mm)

res_name <- resultsNames(dds_alv_mm)[2];res_name
res_alv_mm <- results(dds_alv_mm, name=res_name)
res_alv_mm$design <- paste0(res_name, ".alv.d21")
res_alv_mm$species <- "mouse"
res_alv_mm_df <- as.data.frame(res_alv_mm)
res_alv_mm_df$gene_mm <- rownames(res_alv_mm_df)
res_alv_mm_df$gene_hs <- "NA"

write.csv(res_alv_mm_df, file.path(DIR_OBJ_OUT, "mm_visium_pseudobulk_dea_res_BLMd21_region_alveolar_vs_all_CTRLd21.csv"), row.names = T)
res_alv_mm_df <- read.csv(file.path(DIR_OBJ_OUT, "mm_visium_pseudobulk_dea_res_BLMd21_region_alveolar_vs_all_CTRLd21.csv"), row.names = 1)


###### DEA human - vs all HC ######
dat_hs_alv_vs_ctrl <- bind_cols(dat_hs_all[,grep("HC_", colnames(dat_hs_all), value = T)], dat_hs_alv)

dds_alv_hs <- DESeqDataSetFromMatrix(countData = dat_hs_alv_vs_ctrl,
                                     colData = mdat_hs[colnames(dat_hs_alv_vs_ctrl),],
                                     design = ~ condition)
dds_alv_hs <- DESeq(dds_alv_hs)

res_name <- resultsNames(dds_alv_hs)[2];res_name
res_alv_hs <- results(dds_alv_hs, name=res_name)
res_alv_hs$design <- res_name
res_alv_hs$species <- "human"
res_alv_hs_df <- as.data.frame(res_alv_hs)
res_alv_hs_df$gene_mm <- "NA"
res_alv_hs_df$gene_hs <- rownames(res_alv_hs_df)

write.csv(res_alv_hs_df, file.path(DIR_OBJ_OUT, "hs_visium_pseudobulk_dea_res_IPF_region_alveolar_vs_all_CTRL.csv"), row.names = T)
res_alv_hs_df <- read.csv(file.path(DIR_OBJ_OUT, "hs_visium_pseudobulk_dea_res_IPF_region_alveolar_vs_all_CTRL.csv"), row.names = 1)


###### DEA human - vs Alv HC ######
dds_alv_only_hs <- DESeqDataSetFromMatrix(countData = dat_hs_alv_all,
                                          colData = mdat_hs[colnames(dat_hs_alv_all),],
                                          design = ~ condition)

dds_alv_only_hs <- DESeq(dds_alv_only_hs)

res_name <- resultsNames(dds_alv_only_hs)[2];res_name
res_alv_only_hs <- results(dds_alv_only_hs, name=res_name)
res_alv_only_hs$design <- res_name
res_alv_only_hs$species <- "human"
res_alv_only_hs_df <- as.data.frame(res_alv_only_hs)
res_alv_only_hs_df$gene_mm <- "NA"
res_alv_only_hs_df$gene_hs <- rownames(res_alv_hs_df)

write.csv(res_alv_only_hs_df, file.path(DIR_OBJ_OUT, "hs_visium_pseudobulk_dea_res_IPF_region_alv_vs_alv_CTRL.csv"), row.names = T)
res_alv_only_hs_df <- read.csv(file.path(DIR_OBJ_OUT, "hs_visium_pseudobulk_dea_res_IPF_region_alv_vs_alv_CTRL.csv"), row.names = 1)


###### Join data ######
# Compare d7 and d21
res_alv_mm_d7d21_df <- full_join(x = res_alv_mm_d7_df %>% mutate(logFC_d7 = log2FoldChange, padj_d7 = padj) %>% select(gene_mm, logFC_d7, padj_d7),
                                 y = res_alv_mm_df %>% mutate(logFC_d21 = log2FoldChange, padj_d21 = padj) %>% select(gene_mm, logFC_d21, padj_d21),
                                 by = "gene_mm") %>%
  as_tibble()

res_alv_mm_d7d21_df <- res_alv_mm_d7d21_df %>% 
  mutate(mm_diff = ifelse((logFC_d7>0 & logFC_d21<0)|(logFC_d7<0 & logFC_d21>0), 1, 0),
         sign_005 = ifelse(padj_d7<0.05 & padj_d21<0.05, TRUE, FALSE),
         sign_001 = ifelse(padj_d7<0.01 & padj_d21<0.01, TRUE, FALSE)) 


# identify by overlapping genes
genes_mm <- rownames(res_alv_mm_df)[rownames(res_alv_mm_df) %in% gene_conv_df$symbol_mm]
# genes_mm <- c(genes_mm, rownames(res_alv_mm_d7)[rownames(res_alv_mm_d7) %in% gene_conv_df$symbol_mm]) %>% unique()
genes_hs <- rownames(res_alv_hs_df)[rownames(res_alv_hs_df) %in% gene_conv_df$symbol_hs]
genes_shared <- intersect(gene_conv_df[gene_conv_df$symbol_mm %in% genes_mm,]$symbol_hs_mm,
                          gene_conv_df[gene_conv_df$symbol_hs %in% genes_hs,]$symbol_hs_mm)
length(genes_shared)

# subset data to shared genes
res_alv_mm_subset <- res_alv_mm_df[gene_conv_df[gene_conv_df$symbol_hs_mm %in% genes_shared,]$symbol_mm, ]
# res_alv_mm_subset_d7 <- res_alv_mm_d7[gene_conv_df[gene_conv_df$symbol_hs_mm %in% genes_shared,]$symbol_mm, ]
res_alv_hs_subset <- res_alv_hs_df[gene_conv_df[gene_conv_df$symbol_hs_mm %in% genes_shared,]$symbol_hs, ]
rownames(res_alv_mm_subset) <- rownames(res_alv_hs_subset) <- genes_shared

res_alv_mm_subset$gene <- res_alv_hs_subset$gene <- genes_shared
res_alv_mm_subset$gene_hs <- res_alv_hs_subset$gene_hs <- gene_conv_df[gene_conv_df$symbol_hs_mm %in% genes_shared, "symbol_hs"]
res_alv_mm_subset$gene_mm <- res_alv_hs_subset$gene_mm <- gene_conv_df[gene_conv_df$symbol_hs_mm %in% genes_shared, "symbol_mm"]

# Join data
# res_alv <- rbind(res_alv_mm_subset, res_alv_mm_subset_d7)
res_alv <- rbind(res_alv_mm_subset, res_alv_hs_subset)

# save
write.csv(res_alv, file.path(DIR_OBJ_OUT, "hs_ms_visium_pseudobulk_dea_res_joined_Control_IPF_or_BLMd21_region_alv.csv"), row.names = T)
res_alv <- read.csv(file.path(DIR_OBJ_OUT, "hs_ms_visium_pseudobulk_dea_res_joined_Control_IPF_or_BLMd21_region_alv.csv"), row.names = 1)

# format
res_alv_df <- res_alv %>% 
  as.data.frame() %>% 
  dplyr::arrange(padj, abs(log2FoldChange))

###### Compare results ###### 
# Look at sign data
res_alv_sign <- res_alv_df %>% 
  dplyr::arrange(padj, abs(log2FoldChange)) %>%
  dplyr::filter(padj<0.01 & abs(log2FoldChange) > 1) # FC>2
dim(res_alv_sign)

unique_mm_d21 <- setdiff(subset(res_alv_sign, species == "mouse" & design == "condition_bleomycin_vs_control.alv.d21")$gene, 
                         subset(res_alv_sign, species == "human")$gene)
unique_hs <- setdiff(subset(res_alv_sign, species == "human")$gene, 
                     subset(res_alv_sign, species == "mouse" & design == "condition_bleomycin_vs_control.alv.d21")$gene)

shared_mm_hs <- intersect(
  subset(res_alv_sign, species == "mouse" & design == "condition_bleomycin_vs_control.alv.d21")$gene,
  subset(res_alv_sign, species == "human")$gene)

length(unique_mm_d21);length(unique_hs);length(shared_mm_hs)

share_df <- data.frame(shared_gene = shared_mm_hs,
                       shared_gene_hs = gene_conv_df[gene_conv_df$symbol_hs_mm %in% shared_mm_hs,][shared_mm_hs,]$symbol_hs,
                       shared_gene_mm = gene_conv_df[gene_conv_df$symbol_hs_mm %in% shared_mm_hs,][shared_mm_hs,]$symbol_mm)
share_df <- share_df %>% arrange(shared_gene)
# View(share_df)


# Data frame summarizing DEGs (IPF vs BLM d21)
res_alv_shared_df <- subset(res_alv_sign, design != "condition_bleomycin_vs_control.d7") %>% 
  filter(gene %in% share_df$shared_gene) %>% 
  arrange(gene) %>% 
  select(gene, species, log2FoldChange) %>% 
  pivot_wider(names_from = species, values_from = log2FoldChange) %>% 
  mutate(up = ifelse(mouse>0 & human>0, 1, 0),
         down = ifelse(mouse<0 & human<0, 1, 0),
         diff = ifelse((mouse>0 & human<0) | (mouse<0 & human>0), 1, 0))
subset(res_alv_shared_df, diff==1)

dea_comp_df[dea_comp_df$area == "alveolar", -1] <- c(
  unique_hs %>% length(),
  unique_mm_d21 %>% length(),
  shared_mm_hs %>% length(),
  sum(res_alv_shared_df$up),
  sum(res_alv_shared_df$down),
  sum(res_alv_shared_df$diff)
)


# Venn: IPF and BLM d21
d_venn2 <- list(
  IPF = subset(res_alv_sign, species == "human")$gene,
  BLM = subset(res_alv_sign, species == "mouse" & design == "condition_bleomycin_vs_control.alv.d21")$gene)

# Save data
venn_df_save <- data.frame(group = c(rep("IPF", length(d_venn2$IPF)), rep("BLM", length(d_venn2$BLM))),
                           gene = c(d_venn2$IPF, d_venn2$BLM))
venn_df_save$overlap <- ifelse(duplicated(venn_df_save$gene)|duplicated(venn_df_save$gene, fromLast=T),
                               "both", 
                               venn_df_save$group)
subset(venn_df_save, overlap=="both" & group=="IPF")

write.csv(venn_df_save,
          file.path(DIR_OBJ_OUT, "hs_ms_visium_pseudobulk_dea_res_joined_Control_IPF_or_BLMd21_region_alv_DEGs_venn.csv"), row.names = T)

# plot
pdf(file.path(DIR_FIG_OUT, "hs_ms_visium_pseudobulk_dea_res_joined_Control_IPF_or_BLMd21_region_alv_DEGs_venn.pdf"), 
    width = 1.5, height = 1.5)
ggVennDiagram(d_venn2, label_alpha = 0,
              # edge_color = "black", 
              edge_lty = "solid", 
              label = "count") +
  scale_color_manual(values = c("#B26694", "#215288")) +
  scale_fill_distiller(palette = "Greys", direction = -1) +
  theme(legend.position = "none")
dev.off()


####### Box plot top DEGs #######
g_plot <- c(
  "BDKRB2",
  "FSTL3",
  "CXCL9",
  "TAP1",
  "PHLDA1",
  "UNC5B",
  "PALD1",
  "ITPKC",
  "SAP30",
  "SYDE1",
  "FN1",
  "ZNF189",
  "HMOX2",
  "OAS3",
  "TNFAIP1",
  "APOL1"
)

# Box plot of top res_alv_only_hs_df
# Box plot of top res_alv_hs_df
g_plot <- bind_rows(
  res_alv_hs_df %>% 
    dplyr::filter(padj<0.01, log2FoldChange > 2) %>% 
    dplyr::arrange(desc(log2FoldChange), padj) %>%
    head(10),
  res_alv_hs_df %>% 
    dplyr::filter(padj<0.01, log2FoldChange < -2) %>% 
    dplyr::arrange((log2FoldChange), padj) %>%
    head(10),
) %>% 
  rownames(g_plot)

res_alv_hs_df[g_plot,]
d_plot <- t(dat_hs_alv_all[g_plot, ])
d_plot <- cbind(d_plot, mdat_hs)

plot_list <- lapply(g_plot, function(g){
  p <- ggplot(d_plot, aes_string(x = "condition", y = paste0("log2(", g,"+1)"), fill = "condition")) +
    geom_boxplot(color="black", size=0.2) +
    geom_point(size=1) +
    labs(y="log2(count+1)", title=g) +
    scale_fill_manual(values = cols_cond) +
    theme_classic() +
    theme(legend.position = "none", 
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size=8),
          axis.text = element_text(size=8),
          plot.title = element_text(hjust=0.5, size=8, face = "bold"))
})
p_box_top_all <- wrap_plots(plot_list, ncol = 5) + 
  patchwork::plot_annotation(title = "DEGs: IPF alv. vs HC alv.", theme = theme(plot.title = element_text(hjust=0.5)))
p_box_top_all


##### Fibrosis ##### 
rownames(dat_mm_fib)[rowSums(dat_mm_fib)>10] %>% length()
rownames(dat_hs_fib)[rowSums(dat_hs_fib)>10] %>% length()

###### DEA ms vs hs ######
dat_mm_fib_blm <- dat_mm_fib[, colnames(dat_mm_fib) %>% grep(pattern = "d21_b", value = T)]
res_fib_df <- RunRegionSpeciesDEA(bulk_data_mm = dat_mm_fib_blm,
                                  bulk_data_hs = dat_hs_fib,
                                  mdata_mm = mdat_mm,
                                  mdata_hs = mdat_hs, 
                                  gene_conv_df = gene_conv_df)

write.csv(res_fib_df, file.path(DIR_OBJ_OUT, "hs_ms_visium_pseudobulk_dea_res_IPF_BLMd21_region_fibrosis.csv"), row.names = T)
res_fib_df <- read.csv(file.path(DIR_OBJ_OUT, "hs_ms_visium_pseudobulk_dea_res_IPF_BLMd21_region_fibrosis.csv"), row.names = 1)

res_fib_sign <- res_fib_df %>% 
  dplyr::arrange(padj, abs(log2FoldChange)) %>%
  dplyr::filter(padj<0.01)

head(res_fib_sign)


####### Box plot top DEGs #######
# Box plot of top res_alveoli_d21_sign
g_plot <- bind_rows(
  res_fib_df %>% 
    dplyr::filter(padj<0.01, log2FoldChange > 2) %>% 
    dplyr::arrange(desc(log2FoldChange), padj) %>%
    head(10),
  res_fib_df %>% 
    dplyr::filter(padj<0.01, log2FoldChange < -2) %>% 
    dplyr::arrange((log2FoldChange), padj) %>%
    head(10),
) %>% 
  rownames(g_plot)

d_plot <- t(dat_fib[g_plot, ])
d_plot <- cbind(d_plot, mdat_fib)

plot_list <- lapply(g_plot, function(g){
  p <- ggplot(d_plot, aes_string(x = "condition", y = paste0("log2(", g,"+1)"), fill = "condition")) +
    geom_boxplot(color="black", size=0.2) +
    geom_point(size=1) +
    labs(y="log2(count+1)", title=g) +
    scale_fill_manual(values = cols_cond) +
    theme_classic() +
    theme(legend.position = "none", 
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size=8),
          axis.text = element_text(size=8),
          plot.title = element_text(hjust=0.5, size=8, face = "bold"))
})
p_box_top_all <- wrap_plots(plot_list, ncol = 5) + 
  patchwork::plot_annotation(title = "Top DEGs: BLM vs IPF", theme = theme(plot.title = element_text(hjust=0.5)))
p_box_top_all


###### DEA mouse ######
#' d7
dat_mm_fib_d7_blm <- dat_mm_fib[,grep("d7_b", colnames(dat_mm_fib), value = T)]
dat_mm_all_d7_ctrl <- dat_mm_all[,grep("d7_c", colnames(dat_mm_all), value = T)]
dat_mm_fib_vs_ctrl_d7 <- bind_cols(dat_mm_all_d7_ctrl, dat_mm_fib_d7_blm)

dds_fib_mm_d7 <- DESeqDataSetFromMatrix(countData = dat_mm_fib_vs_ctrl_d7,
                                        colData = mdat_mm[colnames(dat_mm_fib_vs_ctrl_d7),],
                                        design = ~ condition)
dds_fib_mm_d7 <- DESeq(dds_fib_mm_d7)

res_name <- resultsNames(dds_fib_mm_d7)[2];res_name
res_fib_mm_d7 <- results(dds_fib_mm_d7, name=res_name)
res_fib_mm_d7$design <- paste0(res_name, ".fib.d7")
res_fib_mm_d7$species <- "mouse"
res_fib_mm_d7_df <- as.data.frame(res_fib_mm_d7)
res_fib_mm_d7_df$gene_mm <- rownames(res_fib_mm_d7_df)
res_fib_mm_d7_df$gene_hs <- "NA"

write.csv(res_fib_mm_d7_df, file.path(DIR_OBJ_OUT, "mm_visium_pseudobulk_dea_res_BLMd7_region_fibrosis_vs_all_CTRLd7.csv"), row.names = T)
res_fib_mm_d7_df <- read.csv(file.path(DIR_OBJ_OUT, "mm_visium_pseudobulk_dea_res_BLMd7_region_fibrosis_vs_all_CTRLd7.csv"), row.names = 1)


#' d21
dat_mm_fib_d21_blm <- dat_mm_fib[,grep("d21_b", colnames(dat_mm_fib), value = T)]
dat_mm_all_d21_ctrl <- dat_mm_all[,grep("d21_c", colnames(dat_mm_all), value = T)]
dat_mm_fib_vs_ctrl_d21 <- bind_cols(dat_mm_all_d21_ctrl, dat_mm_fib_d21_blm)

dds_fib_mm <- DESeqDataSetFromMatrix(countData = dat_mm_fib_vs_ctrl_d21,
                                     colData = mdat_mm[colnames(dat_mm_fib_vs_ctrl_d21),],
                                     design = ~ condition)
dds_fib_mm <- DESeq(dds_fib_mm)

res_name <- resultsNames(dds_fib_mm)[2];res_name
res_fib_mm <- results(dds_fib_mm, name=res_name)
res_fib_mm$design <- paste0(res_name, ".fib.d21")
res_fib_mm$species <- "mouse"
res_fib_mm_df <- as.data.frame(res_fib_mm)
res_fib_mm_df$gene_mm <- rownames(res_fib_mm_df)
res_fib_mm_df$gene_hs <- "NA"

write.csv(res_fib_mm_df, file.path(DIR_OBJ_OUT, "mm_visium_pseudobulk_dea_res_BLMd21_region_fibrosis_vs_all_CTRLd21.csv"), row.names = T)
res_fib_mm_df <- read.csv(file.path(DIR_OBJ_OUT, "mm_visium_pseudobulk_dea_res_BLMd21_region_fibrosis_vs_all_CTRLd21.csv"), row.names = 1)
res_fib_mm_d21_df <- res_fib_mm_df

###### DEA human - vs all HC ######
dat_hs_fib_vs_ctrl <- bind_cols(dat_hs_all[,grep("HC_", colnames(dat_hs_all), value = T)], dat_hs_fib)

dds_fib_hs <- DESeqDataSetFromMatrix(countData = dat_hs_fib_vs_ctrl,
                                     colData = mdat_hs[colnames(dat_hs_fib_vs_ctrl),],
                                     design = ~ condition)
dds_fib_hs <- DESeq(dds_fib_hs)

res_name <- resultsNames(dds_fib_hs)[2];res_name
res_fib_hs <- results(dds_fib_hs, name=res_name)
res_fib_hs$design <- res_name
res_fib_hs$species <- "human"
res_fib_hs_df <- as.data.frame(res_fib_hs)
res_fib_hs_df$gene_mm <- "NA"
res_fib_hs_df$gene_hs <- rownames(res_fib_hs_df)

write.csv(res_fib_hs_df, file.path(DIR_OBJ_OUT, "hs_visium_pseudobulk_dea_res_IPF_region_fibrosis_vs_all_CTRL.csv"), row.names = T)
res_fib_hs_df <- read.csv(file.path(DIR_OBJ_OUT, "hs_visium_pseudobulk_dea_res_IPF_region_fibrosis_vs_all_CTRL.csv"), row.names = 1)


###### DEA human - vs Alv HC ######
dat_hs_fib_vs_alv_ctrl <- bind_cols(dat_hs_alv_all[,grep("HC_", colnames(dat_hs_alv_all), value = T)], dat_hs_fib)

dds_fib_alv_hs <- DESeqDataSetFromMatrix(countData = dat_hs_fib_vs_alv_ctrl,
                                         colData = mdat_hs[colnames(dat_hs_fib_vs_alv_ctrl),],
                                         design = ~ condition)
dds_fib_alv_hs <- DESeq(dds_fib_alv_hs)

res_name <- resultsNames(dds_fib_alv_hs)[2];res_name
res_fib_alv_hs <- results(dds_fib_alv_hs, name=res_name)
res_fib_alv_hs$design <- res_name
res_fib_alv_hs$species <- "human"
res_fib_alv_hs_df <- as.data.frame(res_fib_alv_hs)
res_fib_alv_hs_df$gene_mm <- "NA"
res_fib_alv_hs_df$gene_hs <- rownames(res_fib_alv_hs_df)

write.csv(res_fib_alv_hs_df, file.path(DIR_OBJ_OUT, "hs_visium_pseudobulk_dea_res_IPF_region_fibrosis_vs_alveolar_CTRL.csv"), row.names = T)
res_fib_alv_hs_df <- read.csv(file.path(DIR_OBJ_OUT, "hs_visium_pseudobulk_dea_res_IPF_region_fibrosis_vs_alveolar_CTRL.csv"), row.names = 1)


###### Join data ######
# IPF - BLM d21
res_fib_mm_df <- rbind(res_fib_mm_d7_df, res_fib_mm_d21_df)

# identify by overlapping genes
genes_mm <- res_fib_mm_df$gene_mm %>% unique() %>% intersect(gene_conv_df$symbol_mm)  # rownames(res_fib_mm_df)[rownames(res_fib_mm_df) %in% gene_conv_df$symbol_mm]
genes_hs <- res_fib_hs_df$gene_hs %>% unique() %>% intersect(gene_conv_df$symbol_hs)  # rownames(res_fib_hs_df)[rownames(res_fib_hs_df) %in% gene_conv_df$symbol_hs]
genes_shared <- intersect(gene_conv_df[gene_conv_df$symbol_mm %in% genes_mm,]$symbol_hs_mm,
                          gene_conv_df[gene_conv_df$symbol_hs %in% genes_hs,]$symbol_hs_mm)
length(genes_shared)

# subset data to shared genes
res_fib_mm_subset <- subset(res_fib_mm_df, gene_mm %in% gene_conv_df[gene_conv_df$symbol_hs_mm %in% genes_shared,]$symbol_mm)
res_fib_hs_subset <- subset(res_fib_hs_df, gene_hs %in% gene_conv_df[gene_conv_df$symbol_hs_mm %in% genes_shared,]$symbol_hs)

res_fib_mm_subset$gene <- setNames(gene_conv_df$symbol_hs_mm, nm = gene_conv_df$symbol_mm)[res_fib_mm_subset$gene_mm]
res_fib_hs_subset$gene <- setNames(gene_conv_df$symbol_hs_mm, nm = gene_conv_df$symbol_hs)[res_fib_hs_subset$gene_hs]

res_fib_mm_subset <- res_fib_mm_subset %>% arrange(design, gene)
res_fib_hs_subset <- res_fib_hs_subset %>% arrange(design, gene)

# res_fib_mm_subset <- res_fib_mm_df[gene_conv_df[gene_conv_df$symbol_hs_mm %in% genes_shared,]$symbol_mm, ]
# res_fib_hs_subset <- res_fib_hs_df[gene_conv_df[gene_conv_df$symbol_hs_mm %in% genes_shared,]$symbol_hs, ]
# rownames(res_fib_mm_subset) <- rownames(res_fib_hs_subset) <- genes_shared
# res_fib_mm_subset$gene <- res_fib_hs_subset$gene <- genes_shared
# res_fib_mm_subset$gene_hs <- res_fib_hs_subset$gene_hs <- gene_conv_df[gene_conv_df$symbol_hs_mm %in% genes_shared, "symbol_hs"]
# res_fib_mm_subset$gene_mm <- res_fib_hs_subset$gene_mm <- gene_conv_df[gene_conv_df$symbol_hs_mm %in% genes_shared, "symbol_mm"]

res_fib_mm_subset$gene_hs <- setNames(gene_conv_df$symbol_hs, nm = gene_conv_df$symbol_hs_mm)[res_fib_mm_subset$gene]
res_fib_hs_subset$gene_mm <- setNames(gene_conv_df$symbol_mm, nm = gene_conv_df$symbol_hs_mm)[res_fib_hs_subset$gene]

# Join data
res_fib_joined <- rbind(res_fib_mm_subset, res_fib_hs_subset)
head(res_fib_joined);dim(res_fib_joined)

# save
write.csv(res_fib_joined, file.path(DIR_OBJ_OUT, "hs_ms_visium_pseudobulk_dea_res_joined_Control_IPF_or_BLMd21_or_BLMd7_region_fib.csv"), row.names = T)
# res_fib_joined <- read.csv(file.path(DIR_OBJ_OUT, "hs_ms_visium_pseudobulk_dea_res_joined_Control_IPF_or_BLMd21_region_fib.csv"), row.names = 1)
res_fib_joined <- read.csv(file.path(DIR_OBJ_OUT, "hs_ms_visium_pseudobulk_dea_res_joined_Control_IPF_or_BLMd21_or_BLMd7_region_fib.csv"), row.names = 1)


# format
res_fib_joined_df <- res_fib_joined %>% 
  as.data.frame() %>% 
  dplyr::arrange(padj, abs(log2FoldChange))


###### Compare results ###### 
# Look at sign data
res_fib_joined_sign <- res_fib_joined_df %>% 
  dplyr::filter(padj<0.01 & abs(log2FoldChange) > 1) %>% # FC>2
  dplyr::arrange(desc(abs(log2FoldChange)), padj)
head(res_fib_joined_sign);dim(res_fib_joined_sign)

# IPF vs d7
unique_mm_d7_fib <- setdiff(
  subset(res_fib_joined_sign, species == "mouse" & design == "condition_bleomycin_vs_control.fib.d7")$gene,
  subset(res_fib_joined_sign, species == "human")$gene)

unique_hs_d7_fib <- setdiff(
  subset(res_fib_joined_sign, species == "human")$gene,
  subset(res_fib_joined_sign, species == "mouse" & design == "condition_bleomycin_vs_control.fib.d7")$gene)

shared_mm_hs_d7_fib <- intersect(
  subset(res_fib_joined_sign, species == "mouse" & design == "condition_bleomycin_vs_control.fib.d7")$gene,
  subset(res_fib_joined_sign, species == "human")$gene)

length(unique_mm_d7_fib);length(unique_hs_d7_fib);length(shared_mm_hs_d7_fib)


share_df_d7_fib <- data.frame(
                      shared_gene = shared_mm_hs_d7_fib,
                      shared_gene_hs = gene_conv_df[gene_conv_df$symbol_hs_mm %in% shared_mm_hs_d7_fib,][shared_mm_hs_d7_fib,]$symbol_hs,
                      shared_gene_mm = gene_conv_df[gene_conv_df$symbol_hs_mm %in% shared_mm_hs_d7_fib,][shared_mm_hs_d7_fib,]$symbol_mm)
share_df_d7_fib <- share_df_d7_fib %>% arrange(shared_gene)


# Data frame summarizing DEGs (IPF vs BLM d21)
res_fib_shared_d7_df <- subset(res_fib_joined_sign, design == "condition_bleomycin_vs_control.fib.d7" | design == "condition_IPF_vs_control") %>% 
  filter(gene %in% share_df_d7_fib$shared_gene) %>% 
  arrange(gene) %>% 
  select(gene, species, log2FoldChange) %>% 
  pivot_wider(names_from = species, values_from = log2FoldChange) %>% 
  mutate(up = ifelse(mouse>0 & human>0, 1, 0),
         down = ifelse(mouse<0 & human<0, 1, 0),
         diff = ifelse((mouse>0 & human<0) | (mouse<0 & human>0), 1, 0))
subset(res_fib_shared_d7_df, diff==1) %>% print(n=40)

write.csv(res_fib_shared_d7_df, file = file.path(DIR_OBJ_OUT, "hs_ms_visium_pseudobulk_dea_res_joined_Control_IPF_or_BLMd7_shared_DEGs.csv"), row.names = F)


# IPF vs d21

# unique_mm_d7 <- setdiff(subset(res_all_sign, species == "mouse" & design == "condition_bleomycin_vs_control.d7")$gene,
#                         subset(res_all_sign, species == "human")$gene)
# unique_hs <- setdiff(subset(res_all_sign, species == "human")$gene,
#                      subset(res_all_sign, species == "mouse" & design == "condition_bleomycin_vs_control.d21")$gene)
# 
# shared_mm_hs <- intersect(subset(res_all_sign, species == "mouse" & design == "condition_bleomycin_vs_control.d21")$gene,
#                           subset(res_all_sign, species == "human")$gene)
# 
# length(unique_mm_d21);length(unique_hs);length(shared_mm_hs)
# 
# share_df <- data.frame(shared_gene = shared_mm_hs,
#                        shared_gene_hs = gene_conv_df[gene_conv_df$symbol_hs_mm %in% shared_mm_hs,][shared_mm_hs,]$symbol_hs,
#                        shared_gene_mm = gene_conv_df[gene_conv_df$symbol_hs_mm %in% shared_mm_hs,][shared_mm_hs,]$symbol_mm)
# share_df <- share_df %>% arrange(shared_gene)
# Look at sign data
# res_fib_sign <- res_fib_df %>% 
#   dplyr::arrange(padj, abs(log2FoldChange)) %>%
#   dplyr::filter(padj<0.01 & abs(log2FoldChange) > 1) # FC>2
# dim(res_fib_sign)

unique_mm_d21_fib <- setdiff(subset(res_fib_joined_sign, species == "mouse" & design == "condition_bleomycin_vs_control.fib.d21")$gene, 
                             subset(res_fib_joined_sign, species == "human")$gene)
unique_hs_d21_fib <- setdiff(subset(res_fib_joined_sign, species == "human")$gene, 
                             subset(res_fib_joined_sign, species == "mouse" & design == "condition_bleomycin_vs_control.fib.d21")$gene)

shared_mm_hs_d21_fib <- intersect(
  subset(res_fib_joined_sign, species == "mouse" & design == "condition_bleomycin_vs_control.fib.d21")$gene,
  subset(res_fib_joined_sign, species == "human")$gene)

length(unique_mm_d21_fib);length(unique_hs_d21_fib);length(shared_mm_hs_d21_fib)


share_df_d21_fib <- data.frame(
                       shared_gene = shared_mm_hs_d21_fib,
                       shared_gene_hs = gene_conv_df[gene_conv_df$symbol_hs_mm %in% shared_mm_hs_d21_fib,][shared_mm_hs_d21_fib,]$symbol_hs,
                       shared_gene_mm = gene_conv_df[gene_conv_df$symbol_hs_mm %in% shared_mm_hs_d21_fib,][shared_mm_hs_d21_fib,]$symbol_mm)
share_df_d21_fib <- share_df_d21_fib %>% arrange(shared_gene)


# Data frame summarizing DEGs (IPF vs BLM d21)
res_fib_shared_d21_df <- subset(res_fib_joined_sign, design == "condition_bleomycin_vs_control.fib.d21" | design == "condition_IPF_vs_control") %>% 
  filter(gene %in% share_df_d21_fib$shared_gene) %>% 
  arrange(gene) %>% 
  select(gene, species, log2FoldChange) %>% 
  pivot_wider(names_from = species, values_from = log2FoldChange) %>% 
  mutate(up = ifelse(mouse>0 & human>0, 1, 0),
         down = ifelse(mouse<0 & human<0, 1, 0),
         diff = ifelse((mouse>0 & human<0) | (mouse<0 & human>0), 1, 0))
subset(res_fib_shared_d21_df, diff==1) %>% print(n=40)

write.csv(res_fib_shared_d21_df, file = file.path(DIR_OBJ_OUT, "hs_ms_visium_pseudobulk_dea_res_joined_Control_IPF_or_BLMd21_shared_DEGs.csv"), row.names = F)

# dea_comp_df[dea_comp_df$area == "fibrosis", -1] <- c(
#   unique_hs %>% length(),
#   unique_mm_d21 %>% length(),
#   shared_mm_hs %>% length(),
#   sum(res_fib_shared_df$up),
#   sum(res_fib_shared_df$down),
#   sum(res_fib_shared_df$diff)
# )
# 
# write.csv(dea_comp_df, file = file.path(DIR_OBJ_OUT, "hs_ms_visium_pseudobulk_dea_res_joined_Control_IPF_or_BLMd21_DEG_comparison.csv"), row.names = F)
# dea_comp_df <- read.csv(file.path(DIR_OBJ_OUT, "hs_ms_visium_pseudobulk_dea_res_joined_Control_IPF_or_BLMd21_DEG_comparison.csv"))


# Plot differing DEGs (IPF, BLMd21)
res_fib_shared_df <- res_fib_shared_d21_df
d_plot <- subset(res_fib_shared_df, diff==1)
d_plot <- d_plot %>% arrange((human))
d_plot$gene <- gsub("_", ":", d_plot$gene)
d_plot$gene <- factor(d_plot$gene, levels = d_plot$gene)

pdf(file.path(DIR_FIG_OUT, "hs_ms_visium_pseudobulk_dea_res_joined_Control_IPF_or_BLMd21_region_fib_DEGs_diff_bars.pdf"), 
  width = 2.25, height = 3)
ggplot(d_plot) +
  geom_col(aes(x=gene, y=human), fill="#B26694") +
  geom_col(aes(x=gene, y=mouse), fill="#215288") +
  geom_hline(yintercept = 0, color="black", linewidth=0.25) +
  scale_x_discrete(position = "top") +
  coord_flip() +
  # ylim(-5,5) +
  labs(y="log2FC", x="") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color="black", size=10),
        axis.text.y = element_text(face = "italic"),
        axis.ticks = element_line(color="black", linewidth=0.25),
        panel.background = element_rect(color="black", linewidth=0.25, fill=NA))
dev.off()


# Plot differing DEGs (opposing directions) - shared (x2) DEGs between IPF-BLMd21 and IPF-BLMd7
genes_plot <- intersect(subset(res_fib_shared_d7_df, diff==1)$gene, subset(res_fib_shared_d21_df, diff==1)$gene)

d_plot <- bind_cols(subset(dplyr::rename(res_fib_shared_d7_df, mouse_d7 = mouse), gene %in% genes_plot) %>% select(gene, human, mouse_d7),
                    subset(dplyr::rename(res_fib_shared_d21_df, mouse_d21 = mouse, human_d21 = human), gene %in% genes_plot) %>% select(mouse_d21, human_d21))
d_plot <- d_plot %>% arrange((human))
d_plot$gene <- gsub("_", ":", d_plot$gene)
d_plot$gene <- factor(d_plot$gene, levels = rev(d_plot$gene))
d_plot_long <- pivot_longer(d_plot, cols = c("mouse_d7", "mouse_d21", "human"), names_to = "group", values_to = "log2FC") %>% select(gene, group, log2FC)
d_plot_long$group <- factor(d_plot_long$group, levels = c("mouse_d7", "mouse_d21", "human"))

custom_cols <- c(hcl.colors(n = length(genes_plot), palette = "inferno")[1:length(genes_plot)-1], "#48C2B4")

pdf(file.path(DIR_FIG_OUT, "hs_ms_visium_pseudobulk_dea_res_joined_Control_IPF_or_BLMd21_or_BLMd7_region_fib_DEGs_diff_lines.pdf"), 
    width = 4, height = 4)
ggplot(d_plot_long, aes(x=group, y=log2FC, color=gene, group=gene)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color="black") +
  geom_line(linetype = "dotted", linewidth=1) +
  geom_point(size=3) +
  geom_text(aes(label=ifelse(group=="human", as.character(gene),'')), hjust=-0.1, vjust=0.5) +
  scale_color_manual(values = custom_cols) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.text = element_text(color="black", size=10),
        axis.ticks = element_line(color="black", linewidth=0.25),
        panel.background = element_rect(color="black", linewidth=0.25, fill=NA),
        legend.position = "none")
dev.off()

write.csv(d_plot_long,
          file.path(DIR_OBJ_OUT, "hs_ms_visium_pseudobulk_dea_res_joined_Control_IPF_or_BLMd21_or_BLMd7_region_fib_DEGs_diff_lines.csv"), row.names = T)


# Plot differing DEGs between IPF-BLMd21 and IPF-BLMd7
intersect(subset(res_fib_shared_d7_df, up==0)$gene, subset(res_fib_shared_d21_df, up==1)$gene)
intersect(subset(res_fib_shared_d7_df, up==1)$gene, subset(res_fib_shared_d21_df, up==0)$gene)
intersect(subset(res_fib_shared_d7_df, diff==1)$gene, subset(res_fib_shared_d21_df, diff==0)$gene)
intersect(subset(res_fib_shared_d7_df, diff==0)$gene, subset(res_fib_shared_d21_df, diff==1)$gene)

subset(res_fib_shared_d7_df, gene %in% "SHH_Shh")
subset(res_fib_shared_d21_df, gene %in% "SHH_Shh")

res_fib_shared_d7d21_df <- full_join(
                              dplyr::rename(res_fib_shared_d7_df, mouse_d7 = mouse) %>% select(gene, mouse_d7),
                              dplyr::rename(res_fib_shared_d21_df, mouse_d21 = mouse) %>% select(gene, mouse_d21, human), 
                              by = "gene")

res_fib_shared_d7d21_df <- res_fib_shared_d7d21_df %>% 
  mutate(
    duu = ifelse(mouse_d7<0 & mouse_d21>0 & human>0, 1, 0),
    udd = ifelse(mouse_d7>0 & mouse_d21<0 & human<0, 1, 0),
    udu = ifelse(mouse_d7>0 & mouse_d21<0 & human>0, 1, 0),
    dud = ifelse(mouse_d7<0 & mouse_d21>0 & human<0, 1, 0),
    mm_diff = ifelse((mouse_d7<0 & mouse_d21>0) | (mouse_d7>0 & mouse_d21<0), 1, 0),
    )
res_fib_shared_d7d21_df_filt <- res_fib_shared_d7d21_df %>% 
  na.omit() %>% 
  mutate(sums = rowSums(across(c(duu,udd,udu,dud))))
res_fib_shared_d7d21_df_filt %>% 
  filter(sums>0)



# Venn: IPF and BLM d21 / BLM d7
d_venn <- list(
  IPF = subset(res_fib_joined_sign, species == "human")$gene,
  BLM_d7 = subset(res_fib_joined_sign, species == "mouse" & design == "condition_bleomycin_vs_control.fib.d7")$gene,
  BLM_d21 = subset(res_fib_joined_sign, species == "mouse" & design == "condition_bleomycin_vs_control.fib.d21")$gene)

p_v1 <- ggVennDiagram(d_venn[c("IPF", "BLM_d7")], label_alpha = 0,
              # edge_color = "black", 
              edge_lty = "solid", 
              label = "count") +
  scale_color_manual(values = c("#B26694", "#215288")) +
  scale_fill_distiller(palette = "Greys", direction = -1) +
  theme(legend.position = "none")

p_v2 <- ggVennDiagram(d_venn[c("IPF", "BLM_d21")], label_alpha = 0,
                      # edge_color = "black", 
                      edge_lty = "solid", 
                      label = "count") +
  scale_color_manual(values = c("#B26694", "#215288")) +
  scale_fill_distiller(palette = "Greys", direction = -1) +
  theme(legend.position = "none")

pdf(file.path(DIR_FIG_OUT, "hs_ms_visium_pseudobulk_dea_res_joined_Control_IPF_or_BLMd21_or_BLMd7_region_fib_DEGs_venn.pdf"), 
    width = 1.5, height = 3)
p_v1 / p_v2
dev.off()


# Venn: IPF and BLM d21
d_venn2 <- list(
  IPF = subset(res_fib_sign, species == "human")$gene,
  BLM = subset(res_fib_sign, species == "mouse" & design == "condition_bleomycin_vs_control.fib.d21")$gene)

# Save data
venn_df_save <- data.frame(group = c(rep("IPF", length(d_venn2$IPF)), rep("BLM", length(d_venn2$BLM))),
                           gene = c(d_venn2$IPF, d_venn2$BLM))
venn_df_save$overlap <- ifelse(duplicated(venn_df_save$gene)|duplicated(venn_df_save$gene, fromLast=T),
                               "both", 
                               venn_df_save$group)
subset(venn_df_save, overlap=="both" & group=="IPF")

write.csv(venn_df_save,
          file.path(DIR_OBJ_OUT, "hs_ms_visium_pseudobulk_dea_res_joined_Control_IPF_or_BLMd21_region_fib_DEGs_venn.csv"), row.names = T)

# plot
pdf(file.path(DIR_FIG_OUT, "hs_ms_visium_pseudobulk_dea_res_joined_Control_IPF_or_BLMd21_region_fib_DEGs_venn.pdf"), 
    width = 1.5, height = 1.5)
ggVennDiagram(d_venn2, label_alpha = 0,
              # edge_color = "black", 
              edge_lty = "solid", 
              label = "count") +
  scale_color_manual(values = c("#B26694", "#215288")) +
  scale_fill_distiller(palette = "Greys", direction = -1) +
  theme(legend.position = "none")
dev.off()


###### Box plot top DEGs #######
g_plot <- c(
  "BDKRB2",
  "FSTL3",
  "CXCL9",
  "TAP1",
  "PHLDA1",
  "UNC5B",
  "PALD1",
  "ITPKC",
  "SAP30",
  "SYDE1",
  "FN1",
  "ZNF189",
  "HMOX2",
  "OAS3",
  "TNFAIP1",
  "APOL1"
)

# Box plot of top res_fib_alv_hs_df
# Box plot of top res_fib_hs_df
g_plot <- bind_rows(
  res_fib_alv_hs_df %>% 
    dplyr::filter(padj<0.01, log2FoldChange > 2) %>% 
    dplyr::arrange(desc(log2FoldChange), padj) %>%
    head(10),
  res_fib_alv_hs_df %>% 
    dplyr::filter(padj<0.01, log2FoldChange < -2) %>% 
    dplyr::arrange((log2FoldChange), padj) %>%
    head(10),
) %>% 
  rownames(g_plot)

d_plot <- t(dat_hs_fib_vs_alv_ctrl_dds[g_plot, ])
d_plot <- cbind(d_plot, mdat_hs)

plot_list <- lapply(g_plot, function(g){
  p <- ggplot(d_plot, aes_string(x = "condition", y = paste0("log2(", g,"+1)"), fill = "condition")) +
    geom_boxplot(color="black", size=0.2) +
    geom_point(size=1) +
    labs(y="log2(count+1)", title=g) +
    scale_fill_manual(values = cols_cond) +
    theme_classic() +
    theme(legend.position = "none", 
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size=8),
          axis.text = element_text(size=8),
          plot.title = element_text(hjust=0.5, size=8, face = "bold"))
})
p_box_top_all <- wrap_plots(plot_list, ncol = 5) + 
  patchwork::plot_annotation(title = "DEGs: IPF fib. vs HC alv.  (norm data)", theme = theme(plot.title = element_text(hjust=0.5)))
p_box_top_all


##### Inflammation ##### 
rownames(dat_mm_infl)[rowSums(dat_mm_infl)>10] %>% length()
rownames(dat_hs_infl)[rowSums(dat_hs_infl)>10] %>% length()

###### DEA ms vs hs ######
res_infl_df <- RunRegionSpeciesDEA(bulk_data_mm = dat_mm_infl,
                                   bulk_data_hs = dat_hs_infl,
                                   mdata_mm = mdat_mm,
                                   mdata_hs = mdat_hs, 
                                   gene_conv_df = gene_conv_df)

write.csv(res_infl_df, file.path(DIR_OBJ_OUT, "hs_ms_visium_pseudobulk_dea_res_IPF_BLMd21_region_inflammation.csv"), row.names = T)
res_infl_df <- read.csv(file.path(DIR_OBJ_OUT, "hs_ms_visium_pseudobulk_dea_res_IPF_BLMd21_region_inflammation.csv"), row.names = 1)

res_infl_sign <- res_infl_df %>% 
  dplyr::arrange(padj, abs(log2FoldChange)) %>%
  dplyr::filter(padj<0.01)

head(res_infl_sign); dim(res_infl_sign)


##### All/Region comparison ###### 
###### DEA ms vs hs ######
degs_all_regions <- list(alveolar = res_alv_sign,
                    fibrosis = res_fib_sign,
                    inflammation = res_infl_sign
                    )
degs_up_ipf <- lapply(names(degs_all_regions), function(r){
  d <- degs_all_regions[[r]]
  d <- d %>% filter(padj<0.01 & log2FoldChange > 2)
  d$region <- r
  return(d)
})
degs_up_ipf <- bind_rows(degs_up_ipf) %>% as_tibble()

degs_up_blm <- lapply(names(degs_all_regions), function(r){
  d <- degs_all_regions[[r]]
  d <- d %>% filter(padj<0.01 & log2FoldChange < -2)
  d$region <- r
  return(d)
})
degs_up_blm <- bind_rows(degs_up_blm) %>% as_tibble()


# plot venn
d_venn_up_ipf <- list(Alveoli = subset(degs_up_ipf, region == "alveolar")$gene, 
                      Fibrosis = subset(degs_up_ipf, region == "fibrosis")$gene,
                      Inflammation = subset(degs_up_ipf, region == "inflammation")$gene)
d_venn_up_blm <- list(Alveoli = subset(degs_up_blm, region == "alveolar")$gene, 
                      Fibrosis = subset(degs_up_blm, region == "fibrosis")$gene,
                      Inflammation = subset(degs_up_blm, region == "inflammation")$gene)

pv1 <- ggVennDiagram(d_venn_up_ipf, edge_color = "black", edge_lty = "dashed") +
  scale_color_manual(values = c("#1F77B4FF", "#B26694", "#E3A946")) +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  labs(title = "DEGs up in IPF") +
  theme(legend.position = "bottom", legend.text = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust=0.5));pv1
pv2 <- ggVennDiagram(d_venn_up_blm, edge_color = "black", edge_lty = "dashed") +
  scale_color_manual(values = c("#1F77B4FF", "#B26694", "#E3A946")) +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  labs(title = "DEGs up in BLM d21") +
  theme(legend.position = "bottom", legend.text = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust=0.5));pv2


png(filename = file.path(DIR_FIG_OUT, "hs_ms_visium_pseudobulk_dea_res_IPF_BLMd21_region_comparison_DEGs_venn.png"), 
    width = 8*fig_res, height = 4*fig_res, res = fig_res)
(pv1|plot_spacer()|pv2) + plot_layout(widths = c(1,0.1,1))
dev.off()


###### DEA region vs all control ######
#' 1. filter data based on shared genes
#' 2. join tables
#' 3. look at shared/differing DEGs between comparisons
genes_res_mm <- rownames(res_alv_mm_df) %>% intersect(gene_conv_df$symbol_mm)
genes_res_hs <- rownames(res_alv_hs_df) %>% intersect(gene_conv_df$symbol_hs)
genes_shared <- subset(gene_conv_df, symbol_mm %in% genes_res_mm & symbol_hs %in% genes_res_hs)

res_alv_vs_ctrl <- bind_rows(
  res_alv_mm_df[genes_shared$symbol_mm, ],
  res_alv_hs_df[genes_shared$symbol_hs, ])
res_alv_vs_ctrl$symbol_hs_mm <- ifelse(res_alv_vs_ctrl$gene_hs %in% genes_shared$symbol_hs | res_alv_vs_ctrl$gene_mm %in% genes_shared$symbol_mm, 
                                       genes_shared$symbol_hs_mm, NA)
res_alv_vs_ctrl$region <- "alveolar"

res_fib_vs_ctrl <- bind_rows(
  res_fib_mm_df[genes_shared$symbol_mm, ],
  res_fib_hs_df[genes_shared$symbol_hs, ])
res_fib_vs_ctrl$symbol_hs_mm <- ifelse(res_fib_vs_ctrl$gene_hs %in% genes_shared$symbol_hs | res_fib_vs_ctrl$gene_mm %in% genes_shared$symbol_mm, 
                                       genes_shared$symbol_hs_mm, NA)
res_fib_vs_ctrl$region <- "fibrosis"

res_region_vs_ctrl <- bind_rows(res_alv_vs_ctrl, res_fib_vs_ctrl)

res_region_vs_ctrl %>% 
  filter(padj < 0.001 & abs(log2FoldChange) > 1) %>% 
  arrange(desc(abs(log2FoldChange)), padj) %>% 
  head(n=15)


##### AbBa NMF factor hi subcluster ######
###### Direct species DEA ###### 
res_f14c0_df <- RunRegionSpeciesDEA(bulk_data_mm = dat_mm_f14c0,
                                    bulk_data_hs = dat_hs_f14c0,
                                    mdata_mm = mdat_mm[colnames(dat_mm_f14c0),],
                                    mdata_hs = mdat_hs[colnames(dat_hs_f14c0),], 
                                    gene_conv_df = gene_conv_df)

write.csv(res_f14c0_df, file.path(DIR_OBJ_OUT, "hs_ms_visium_pseudobulk_dea_res_IPF_BLMd21_NMF30-F14-C0.csv"), row.names = T)
res_f14c0_df <- read.csv(file.path(DIR_OBJ_OUT, "hs_ms_visium_pseudobulk_dea_res_IPF_BLMd21_NMF30-F14-C0.csv"), row.names = 1)

res_f14c0_sign <- res_f14c0_df %>% 
  dplyr::arrange(padj, abs(log2FoldChange)) %>%
  dplyr::filter(padj<0.01)
res_f14c0_sign$gene <- gsub("_", ":", res_f14c0_sign$gene)
res_f14c0_sign$updown <- ifelse(res_f14c0_sign$log2FoldChange>0, "up", "down")
head(res_f14c0_sign);dim(res_f14c0_sign)

###### Plot results ###### 
#' Least significant genes (shared genes)
res_f14c0_ns <- res_f14c0_df %>%
  dplyr::filter(padj>0.05) %>% 
  dplyr::arrange(desc(padj))


#' Box plot top DEGs
g_plot <- bind_rows(
  res_f14c0_df %>% 
    dplyr::filter(padj<0.01, log2FoldChange > 2) %>% 
    dplyr::arrange(desc(log2FoldChange), padj) %>%
    head(10),
  res_f14c0_df %>% 
    dplyr::filter(padj<0.01, log2FoldChange < -2) %>% 
    dplyr::arrange((log2FoldChange), padj) %>%
    head(10),
) %>% 
  rownames(g_plot)

d_plot_mm <- t(dat_mm_f14c0[gene_conv_df[gene_conv_df$symbol_hs_mm %in% g_plot, "symbol_mm"],]) %>% as.data.frame()
d_plot_mm$condition <- "mm_BLM_d21"
d_plot_hs <- t(dat_hs_f14c0[gene_conv_df[gene_conv_df$symbol_hs_mm %in% g_plot, "symbol_hs"],]) %>% as.data.frame()
d_plot_hs$condition <- "hs_IPF"
colnames(d_plot_mm) <- colnames(d_plot_hs) <- c(g_plot, "condition")
d_plot <- rbind(d_plot_mm, d_plot_hs)

plot_list <- lapply(g_plot, function(g){
  p <- ggplot(d_plot, aes_string(x = "condition", y = paste0("log2(", g,"+1)"), fill = "condition")) +
    geom_boxplot(color="black", size=0.2) +
    geom_point(size=1) +
    labs(y="log2(count+1)", title=g) +
    scale_fill_manual(values = c("#B26694", "#1F77B4FF")) +
    theme_classic() +
    theme(legend.position = "none", 
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size=8),
          axis.text = element_text(size=8),
          plot.title = element_text(hjust=0.5, size=8, face = "bold"))
})
p_box_top_all <- wrap_plots(plot_list, ncol = 5) + 
  patchwork::plot_annotation(title = "Top DEGs: BLM vs IPF", theme = theme(plot.title = element_text(hjust=0.5)))
p_box_top_all

#' Plot heatmap and volcano
res_f14c0_sign_plot <- res_f14c0_sign
res_f14c0_sign_plot <- res_f14c0_sign_plot %>% 
  group_by(updown) %>% 
  slice_max(n=60, order_by = abs(log2FoldChange))
res_f14c0_sign_plot$gene <- factor(res_f14c0_sign_plot$gene, levels = res_f14c0_sign_plot$gene)

lim_max <- round(max(abs(res_f14c0_sign_plot$log2FoldChange))+0.5);lim_max
p1 <- ggplot(res_f14c0_sign_plot, aes(x=reorder(gene, log2FoldChange), y="log2FC", fill=log2FoldChange)) +
  geom_tile(width=0.9) +
  scale_fill_gradientn(colors = col_scale_div_custom2, limits=c(-lim_max, lim_max)) +
  theme_minimal() +
  theme(legend.direction="horizontal", legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(), 
        axis.text.x = element_text(angle=90, hjust = 1))

p_vol <- ggplot(res_f14c0_sign, aes(x = log2FoldChange, y = -log10(padj), color = log2FoldChange)) +
  geom_hline(yintercept = 0, linetype="solid", color = "grey90") +
  geom_vline(xintercept = 0, linetype="solid", color = "grey90") +
  # geom_hline(yintercept = 2, linetype="dashed", color = "grey90") +
  # geom_vline(xintercept = c(-2, 2), linetype="dashed", color = "grey90") +
  geom_point(size = .5, alpha = 0.8) +
  scale_color_gradientn(colors = col_scale_div_custom2, limits=c(-lim_max, lim_max)) +
  # scale_color_manual(values = cols_sign_up_down) +
  # xlim(c(-minmax_xilm, minmax_xilm)) +
  labs(title = "DEGs BLMd21 F14-C0 vs IPF F14-C0") +
  theme_classic() +
  theme(aspect.ratio = 1,
        legend.position = "none",
        plot.title = element_text(size = rel(1), hjust = 0.5),
        axis.title = element_text(size = rel(1)),
        axis.line = element_line(colour = "black"));p_vol


#' Plot barplot
n_genes <- 30
d_pos <- res_f14c0_sign_plot %>% filter(log2FoldChange>0) %>% head(n_genes)
d_neg <- res_f14c0_sign_plot %>% filter(log2FoldChange<0) %>% head(n_genes)
lim_max <- round(max(d_pos$log2FoldChange)+0.5);lim_max
lim_min <- round(min(d_neg$log2FoldChange)-0.5);lim_min
axis_extend <- 5

p_pos <- ggplot(d_pos, aes(x=reorder(gene_hs, log2FoldChange), y=log2FoldChange, fill = log2FoldChange)) +
  geom_col(width = 0.8) +
  geom_text(aes(label = gene_hs), hjust = 0, colour = "black", size=2.8, nudge_y = 0.25) +
  scale_fill_gradientn(colors = col_scale_div_custom2, limits=c(-11, lim_max)) +
  scale_x_discrete(position = "top") +
  scale_y_continuous(position = "right", limits = c(0, lim_max+axis_extend)) +
  coord_flip() +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 8, color="black"),
        axis.text.y = element_blank(),
        axis.line.x = element_line(linewidth = 0.25, color="black"), 
        axis.ticks.x = element_line(linewidth = 0.25, color="black"),
        axis.title = element_blank(), 
        legend.position = "none")

p_neg <- ggplot(d_neg, aes(x=reorder(gene_mm, abs(log2FoldChange)), y=log2FoldChange, fill = log2FoldChange)) +
  geom_col(width = 0.8) +
  geom_text(aes(label = gene_mm), hjust = 1, colour = "black", size=2.8, nudge_y = -0.25) +
  scale_fill_gradientn(colors = col_scale_div_custom2, limits=c(-11, lim_max)) +
  scale_y_continuous(position = "right", limits = c(lim_min-axis_extend, 0)) +
  coord_flip() +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 8, color="black"),
        axis.text.y = element_blank(),
        axis.line.x = element_line(linewidth = 0.25), 
        axis.ticks.x = element_line(linewidth = 0.25),
        axis.title = element_blank(), 
        legend.position = "none")

(p_neg|plot_spacer()|p_pos) + plot_layout(widths = c(1,-0.2,1))

png(filename = file.path(DIR_FIG_OUT, "hs_ms_visium_pseudobulk_dea_res_IPF_BLMd21_NMF30-F14-C0_topDEG_bars.png"), 
    width = 3*fig_res, height = 3.8*fig_res, res = fig_res)
((p_neg|plot_spacer()|p_pos) + plot_layout(widths = c(1,-0.2,1))) + 
  plot_annotation(title = "Log2 fold change", theme = theme(plot.title = element_text(hjust=0.5, size=10)))
dev.off()



###### Region vs Control, species separately ###### 
###### DEA mouse ######
#' d21
dat_mm_f14c0_d21 <- bind_cols(dat_mm_all[,grep("d21_c", colnames(dat_mm_all), value = T)], 
                              dat_mm_f14c0)

dds_f14c0_mm <- DESeqDataSetFromMatrix(countData = dat_mm_f14c0_d21,
                                     colData = mdat_mm[colnames(dat_mm_f14c0_d21),],
                                     design = ~ condition)
dds_f14c0_mm <- DESeq(dds_f14c0_mm)

res_name <- resultsNames(dds_f14c0_mm)[2];res_name
res_f14c0_mm <- results(dds_f14c0_mm, name=res_name)
res_f14c0_mm$design <- paste0(res_name, ".d21")
res_f14c0_mm$species <- "mouse"

###### DEA human ######
dat_hs_f14c0 <- bind_cols(dat_hs_all[,grep("HC_", colnames(dat_hs_all), value = T)],
                          dat_hs_f14c0)
dds_f14c0_hs <- DESeqDataSetFromMatrix(countData = dat_hs_f14c0,
                                     colData = mdat_hs[colnames(dat_hs_f14c0),],
                                     design = ~ condition)
dds_f14c0_hs <- DESeq(dds_f14c0_hs)

res_name <- resultsNames(dds_f14c0_hs)[2];res_name
res_f14c0_hs <- results(dds_f14c0_hs, name=res_name)
res_f14c0_hs$design <- res_name
res_f14c0_hs$species <- "human"


###### Join data ######
# identify by overlapping genes
genes_mm <- rownames(res_f14c0_mm)[rownames(res_f14c0_mm) %in% gene_conv_df$symbol_mm]
genes_hs <- rownames(res_f14c0_hs)[rownames(res_f14c0_hs) %in% gene_conv_df$symbol_hs]
genes_shared <- intersect(gene_conv_df[gene_conv_df$symbol_mm %in% genes_mm,]$symbol_hs_mm,
                          gene_conv_df[gene_conv_df$symbol_hs %in% genes_hs,]$symbol_hs_mm)
length(genes_shared)

# subset data to shared genes
res_f14c0_mm_subset <- res_f14c0_mm[gene_conv_df[gene_conv_df$symbol_hs_mm %in% genes_shared,]$symbol_mm, ]
res_f14c0_hs_subset <- res_f14c0_hs[gene_conv_df[gene_conv_df$symbol_hs_mm %in% genes_shared,]$symbol_hs, ]
rownames(res_f14c0_mm_subset) <- rownames(res_f14c0_hs_subset) <- genes_shared

res_f14c0_mm_subset$gene <- res_f14c0_hs_subset$gene <- genes_shared
res_f14c0_mm_subset$gene_hs <- res_f14c0_hs_subset$gene_hs <- gene_conv_df[gene_conv_df$symbol_hs_mm %in% genes_shared, "symbol_hs"]
res_f14c0_mm_subset$gene_mm <- res_f14c0_hs_subset$gene_mm <- gene_conv_df[gene_conv_df$symbol_hs_mm %in% genes_shared, "symbol_mm"]

# Join data
res_f14c0 <- rbind(res_f14c0_mm_subset, res_f14c0_hs_subset)

# save
write.csv(res_f14c0, file.path(DIR_OBJ_OUT, "hs_ms_visium_pseudobulk_dea_res_joined_Control_IPF_or_BLMd21_region_f14c0.csv"), row.names = T)
res_f14c0 <- read.csv(file.path(DIR_OBJ_OUT, "hs_ms_visium_pseudobulk_dea_res_joined_Control_IPF_or_BLMd21_region_f14c0.csv"), row.names = 1)

# format
res_f14c0_df <- res_f14c0 %>% 
  as.data.frame() %>% 
  dplyr::arrange(padj, abs(log2FoldChange))


###### Compare results ###### 
# Look at sign data
res_f14c0_sign <- res_f14c0_df %>% 
  dplyr::arrange(padj, abs(log2FoldChange)) %>%
  dplyr::filter(padj<0.01 & abs(log2FoldChange) > 1) # FC>2
dim(res_f14c0_sign)

# IPF, d21
unique_mm_d21 <- setdiff(subset(res_f14c0_sign, species == "mouse" & design == "condition_bleomycin_vs_control.d21")$gene,
                         subset(res_f14c0_sign, species == "human")$gene)
unique_hs <- setdiff(subset(res_f14c0_sign, species == "human")$gene,
                     subset(res_f14c0_sign, species == "mouse" & design == "condition_bleomycin_vs_control.d21")$gene)

shared_mm_hs <- intersect(subset(res_f14c0_sign, species == "mouse" & design == "condition_bleomycin_vs_control.d21")$gene,
                          subset(res_f14c0_sign, species == "human")$gene)

length(unique_mm_d21);length(unique_hs);length(shared_mm_hs)

share_df <- data.frame(shared_gene = shared_mm_hs,
                       shared_gene_hs = gene_conv_df[gene_conv_df$symbol_hs_mm %in% shared_mm_hs,][shared_mm_hs,]$symbol_hs,
                       shared_gene_mm = gene_conv_df[gene_conv_df$symbol_hs_mm %in% shared_mm_hs,][shared_mm_hs,]$symbol_mm)
share_df <- share_df %>% arrange(shared_gene)
View(share_df)

# Data frame summarizing DEGs (IPF vs BLM d21)
dea_comp_df2 <- data.frame(area = c("f14C0"),
                          DEGs_hs = rep(0,1),
                          DEGs_mm = rep(0,1),
                          DEGs_shared = rep(0,1),
                          shared_up = rep(0,1),
                          shared_down = rep(0,1),
                          shared_diff = rep(0,1))


res_f14c0_shared_df <- subset(res_f14c0_sign, design != "condition_bleomycin_vs_control.d7") %>% 
  filter(gene %in% share_df$shared_gene) %>% 
  arrange(gene) %>% 
  select(gene, species, log2FoldChange) %>% 
  pivot_wider(names_from = species, values_from = log2FoldChange) %>% 
  mutate(up = ifelse(mouse>0 & human>0, 1, 0),
         down = ifelse(mouse<0 & human<0, 1, 0),
         diff = ifelse((mouse>0 & human<0) | (mouse<0 & human>0), 1, 0))

dea_comp_df2[dea_comp_df2$area == "f14C0", -1] <- c(
  unique_hs %>% length(),
  unique_mm_d21 %>% length(),
  shared_mm_hs %>% length(),
  sum(res_all_shared_df$up),
  sum(res_all_shared_df$down),
  sum(res_all_shared_df$diff)
)
dea_comp_df <- bind_rows(dea_comp_df, dea_comp_df2)

# look at closer at shared genes
res_f14c0_shared_df[res_f14c0_shared_df$diff==1,]

res_f14c0_shared_df <- res_f14c0_shared_df %>% arrange((human))
res_f14c0_shared_df$gene <- gsub("_", ":", res_f14c0_shared_df$gene)
res_f14c0_shared_df$gene <- factor(res_f14c0_shared_df$gene, levels = res_f14c0_shared_df$gene)

# Plot differing genes
# d_plot <- subset(res_f14c0_shared_df, diff==1)
# d_plot <- d_plot %>% arrange((human))
# d_plot$gene <- gsub("_", ":", d_plot$gene)
# d_plot$gene <- factor(d_plot$gene, levels = d_plot$gene)

theme_custom <- theme(panel.grid = element_blank(),
                      axis.text = element_text(color="black", size=10),
                      axis.text.y = element_text(face = "italic"),
                      axis.ticks = element_line(color="black", linewidth=0.25),
                      panel.background = element_rect(color="black", linewidth=0.25, fill=NA))

p_up <- ggplot(subset(res_f14c0_shared_df, up==1)) +
  geom_col(aes(x=gene, y=human), fill="#B26694", alpha = 1) +
  geom_col(aes(x=gene, y=mouse), fill="#215288", alpha = 0.7) +
  geom_hline(yintercept = 0, color="black", linewidth=0.25) +
  scale_x_discrete(position = "top") +
  coord_flip() +
  labs(y="log2FC", x="") +
  theme_bw() + theme_custom

p_diff <- ggplot(subset(res_f14c0_shared_df, diff==1)) +
  geom_col(aes(x=gene, y=human), fill="#B26694") +
  geom_col(aes(x=gene, y=mouse), fill="#215288") +
  geom_hline(yintercept = 0, color="black", linewidth=0.25) +
  scale_x_discrete(position = "top") +
  coord_flip() +
  # ylim(-5,5) +
  labs(y="log2FC", x="") +
  theme_bw() + theme_custom

p_down <- ggplot(subset(res_f14c0_shared_df, down==1)) +
  geom_col(aes(x=gene, y=human), fill="#B26694", alpha = 1) +
  geom_col(aes(x=gene, y=mouse), fill="#215288", alpha = 0.7) +
  geom_hline(yintercept = 0, color="black", linewidth=0.25) +
  scale_x_discrete(position = "top") +
  coord_flip() +
  labs(y="log2FC", x="") +
  theme_bw() + theme_custom

pdf(file.path(DIR_FIG_OUT, "hs_ms_visium_pseudobulk_dea_res_joined_Control_IPF_or_BLMd21_region_f14c0_DEGs_shared_bars.pdf"), 
    width = 2.25*4, height = 6)
((p_down|p_diff|p_up)+plot_layout(widths = c(1,1.5,1)))+plot_annotation(title = "Shared DEGs within f14c0 region DEA (vs control)", theme = theme(plot.title=element_text(hjust=0.5)))
dev.off()


# library(ggVennDiagram)
# Venn: IPF and BLM d21
d_venn2 <- list(
  IPF = subset(res_f14c0_sign, species == "human")$gene,
  BLM = subset(res_f14c0_sign, species == "mouse" & design == "condition_bleomycin_vs_control.d21")$gene)

# Save data
venn_df_save <- data.frame(group = c(rep("IPF", length(d_venn2$IPF)), rep("BLM", length(d_venn2$BLM))),
                           gene = c(d_venn2$IPF, d_venn2$BLM))
venn_df_save$overlap <- ifelse(duplicated(venn_df_save$gene)|duplicated(venn_df_save$gene, fromLast=T),
                               "both", 
                               venn_df_save$group)
write.csv(venn_df_save,
          file.path(DIR_OBJ_OUT, "hs_ms_visium_pseudobulk_dea_res_joined_Control_IPF_or_BLMd21_region_f14c0_DEGs_venn.csv"), row.names = T)

# Plot
pdf(file.path(DIR_FIG_OUT, "hs_ms_visium_pseudobulk_dea_res_joined_Control_IPF_or_BLMd21_region_f14c0_DEGs_venn.pdf"), 
    width = 1.5, height = 1.5)
ggVennDiagram(d_venn2, label_alpha = 0,
              # edge_color = "black", 
              edge_lty = "solid", 
              label = "count") +
  scale_color_manual(values = c("#B26694", "#215288")) +
  scale_fill_distiller(palette = "Greys", direction = -1) +
  theme(legend.position = "none")
dev.off()


#### Pathway enrichment of DEGS #### 
#' Using g:Profiler (and fGSVA?)
org <- "mmusculus"
org <- "hsapiens"

##### Alveolar ###### 
gene_query_list <- list(
  up_IPF = gene_conv_df[gene_conv_df$symbol_hs_mm %in% rownames(res_alv_sign[res_alv_sign$log2FoldChange > 2,]), "symbol_hs"],
  up_BLM = gene_conv_df[gene_conv_df$symbol_hs_mm %in% rownames(res_alv_sign[res_alv_sign$log2FoldChange < -2,]), "symbol_hs"]
)

gostres_list <- list()
for (f in names(gene_query_list)){
  message(f)
  gostres_list[[f]] <- gprofiler2::gost(query = gene_query_list[[f]],
                                        organism = org, 
                                        ordered_query = F)
}

gostplot(gostres_list$up_IPF)
gostplot(gostres_list$up_BLM)

gostres_export_list <- list()
for (f in names(gostres_list)){
  gostres_export_list[[f]] <- gostres_list[[f]]$result
  gostres_export_list[[f]]$query <- f
}

gostres_list_alv <- gostres_list
gostres_export_list_alv <- gostres_export_list

write_xlsx(
  x = gostres_export_list_alv,
  path = file.path(DIR_OBJ_OUT, paste0("hs_ms_visium_pseudobulk_dea_res_IPF_BLMd21_deg_gprofiler_Hs_region_alveolar.csv", ".xlsx")),
  col_names = TRUE,
  format_headers = TRUE,
  use_zip64 = FALSE
)

##### Fibrosis ###### 
# gene_query_list <- list(
#   up_IPF = gene_conv_df[gene_conv_df$symbol_hs_mm %in% rownames(res_fib_sign[res_fib_sign$log2FoldChange > 2,]), "symbol_mm"],
#   up_BLM = gene_conv_df[gene_conv_df$symbol_hs_mm %in% rownames(res_fib_sign[res_fib_sign$log2FoldChange < -2,]), "symbol_mm"]
# )
gene_query_list <- list(
  up_IPF = gene_conv_df[gene_conv_df$symbol_hs_mm %in% rownames(res_fib_sign[res_fib_sign$log2FoldChange > 2,]), "symbol_hs"],
  up_BLM = gene_conv_df[gene_conv_df$symbol_hs_mm %in% rownames(res_fib_sign[res_fib_sign$log2FoldChange < -2,]), "symbol_hs"]
)

gostres_list <- list()
for (f in names(gene_query_list)){
  message(f)
  gostres_list[[f]] <- gprofiler2::gost(query = gene_query_list[[f]],
                                        organism = org, 
                                        ordered_query = F)
}

gostplot(gostres_list$up_IPF)
gostplot(gostres_list$up_BLM)

gostres_export_list <- list()
for (f in names(gostres_list)){
  gostres_export_list[[f]] <- gostres_list[[f]]$result
  gostres_export_list[[f]]$query <- f
}

gostres_list_fib <- gostres_list
gostres_export_list_fib <- gostres_export_list

write_xlsx(
  x = gostres_export_list_fib,
  path = file.path(DIR_OBJ_OUT, paste0("hs_ms_visium_pseudobulk_dea_res_IPF_BLMd21_deg_gprofiler_Hs_region_fibrosis.csv", ".xlsx")),
  col_names = TRUE,
  format_headers = TRUE,
  use_zip64 = FALSE
)


##### Inflammation ###### 
gene_query_list <- list(
  up_IPF = gene_conv_df[gene_conv_df$symbol_hs_mm %in% rownames(res_infl_sign[res_infl_sign$log2FoldChange > 2,]), "symbol_hs"],
  up_BLM = gene_conv_df[gene_conv_df$symbol_hs_mm %in% rownames(res_infl_sign[res_infl_sign$log2FoldChange < -2,]), "symbol_hs"]
)

gostres_list <- list()
for (f in names(gene_query_list)){
  message(f)
  gostres_list[[f]] <- gprofiler2::gost(query = gene_query_list[[f]],
                                        organism = org, 
                                        ordered_query = F)
}

gostplot(gostres_list$up_IPF)
gostplot(gostres_list$up_BLM)

gostres_export_list <- list()
for (f in names(gostres_list)){
  gostres_export_list[[f]] <- gostres_list[[f]]$result
  gostres_export_list[[f]]$query <- f
}

gostres_list_infl <- gostres_list
gostres_export_list_infl <- gostres_export_list

write_xlsx(
  x = gostres_export_list_infl,
  path = file.path(DIR_OBJ_OUT, paste0("hs_ms_visium_pseudobulk_dea_res_IPF_BLMd21_deg_gprofiler_Hs_region_inflammtion.csv", ".xlsx")),
  col_names = TRUE,
  format_headers = TRUE,
  use_zip64 = FALSE
)

##### Comparison ###### 
# plot venn
d_venn_pw_alv <- list(up_BLM = gostres_list_alv$up_BLM$result$term_name,
                      up_IPF = gostres_list_alv$up_IPF$result$term_name)
d_venn_pw_fib <- list(up_BLM = gostres_list_fib$up_BLM$result$term_name,
                      up_IPF = gostres_list_fib$up_IPF$result$term_name)
d_venn_pw_infl <- list(up_BLM = gostres_list_infl$up_BLM$result$term_name,
                      up_IPF = gostres_list_infl$up_IPF$result$term_name)

pv1 <- ggVennDiagram(d_venn_pw_alv, edge_color = "black", edge_lty = "dashed") +
  scale_color_manual(values = c("#1F77B4FF", "#B26694", "#E3A946")) +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  labs(title = "Alveolar region", subtitle = "gProfiler enriched terms") +
  theme(legend.position = "bottom", 
        legend.text = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5));pv1

pv2 <- ggVennDiagram(d_venn_pw_fib, edge_color = "black", edge_lty = "dashed") +
  scale_color_manual(values = c("#1F77B4FF", "#B26694", "#E3A946")) +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  labs(title = "Fibrosis region", subtitle = "gProfiler enriched terms") +
  theme(legend.position = "bottom", 
        legend.text = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5));pv2

pv3 <- ggVennDiagram(d_venn_pw_infl, edge_color = "black", edge_lty = "dashed") +
  scale_color_manual(values = c("#1F77B4FF", "#B26694", "#E3A946")) +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  labs(title = "Inflammation region", subtitle = "gProfiler enriched terms") +
  theme(legend.position = "bottom", 
        legend.text = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5));pv3

png(filename = file.path(DIR_FIG_OUT, "hs_ms_visium_pseudobulk_dea_res_IPF_BLMd21_region_comparison_gprofiler_Hs_terms_venn.png"), 
    width = 12*fig_res, height = 4*fig_res, res = fig_res)
pv1|pv2|pv3
dev.off()






