#' [hs_visium_nichenet.R]
#'
#' NicheNet analysis in human IPF Visium data
#'
#' Based on vignette for NicheNet analysis on Seurat objects
#' https://github.com/saeyslab/nichenetr/blob/master/vignettes/seurat_steps.md
#' 
#' NicheNet analyses performed in four rounds (1-4)
#' Round 1) receiver: hsNMF-F14hi-C0 | senders: nb. clusters 0-5
#' Round 2) receiver: nb. cluster 0 ("F14C0_nbs_C0") (AT cells) | sender: F14hi-C0
#' Round 3) receiver: nb. cluster 1 ("F14C0_nbs_C1") (Fib cells) | sender: F14hi-C0
#' Round 4) receiver: nb. cluster 2 ("F14C0_nbs_C2") (Mac cells) | sender: F14hi-C0
#'
#' And lastly, all results are combined and compared.
#'
#' Oct 2022 - Nov 2023, L. Franzén [lovisa.franzen@scilifelab.se] 

#### Set up ####
##### Define params. ####
set.seed(1)
SPECIES <- "human"
DIR_ROOT <- getwd()
DIR_DATA <- file.path(DIR_ROOT, "data", SPECIES, "visium")
DIR_RES <- file.path(DIR_ROOT, "results", SPECIES)
DIR_FIG <- file.path(DIR_RES, "figures")
DIR_NN <- file.path(DIR_FIG, "NicheNet_A")
DIR_OBJ_NN <- file.path(DIR_RES, "objects", "NicheNet")
fig_res <- 300

##### Load libs ####
library(STutility)
library(patchwork)
library(dplyr)
library(tidyverse)
library(tidyr)
library(writexl)
library(pheatmap)
library(nichenetr)

##### Other ####
# source(file.path(DIR_ROOT, "scripts", "colors.R"))
source(file.path(DIR_ROOT, "scripts", "custom_functions.R"))
source(file.path(DIR_ROOT, "scripts", "custom_colors.R"))
theme_custom <- theme(axis.title.x = element_blank())
col_blue <- "#2171B5"
col_orange <- "#D94701"

##### Read objects ####
# fname <- paste0("hs_visium_preproc_A_se_obj.rds")
# se.subset <- readRDS(file = file.path(DIR_RES, "objects", fname))
fname <- paste0("hs_visium_preproc_A_se_obj_nmf.rds")
se.subset <- readRDS(file = file.path(DIR_RES, "objects", fname))


##### Install NicheNet #####
# Install dependencies:
# BiocManager::install("limma")
# BiocManager::install("ComplexHeatmap")
# urlPackage <- "https://cran.r-project.org/src/contrib/Archive/randomForest/randomForest_4.6-14.tar.gz"
# install.packages(urlPackage, repos=NULL, type="source")  # install.packages("randomForest")
# install.packages("lifecycle")
# install.packages("recipes")
# install.packages("tidymodels")

# Install NicheNet:
# devtools::install_github("saeyslab/nichenetr")  # install using devtools

# Sys.unsetenv("GITHUB_PAT")  # may need to unset GITHUB_PAT
# renv::install("saeyslab/nichenetr")  # install via renv
# library(nichenetr)


#### NicheNet Set Up  ####
#' Following this vignette: 
#' https://github.com/saeyslab/nichenetr/blob/master/vignettes/seurat_steps.md

##### Read R-L priors ##### 
# ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
# lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
# weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))


# dir.create(path = DIR_OBJ_NN)  # Saved 2023-03-01
# saveRDS(ligand_target_matrix, file = file.path(DIR_OBJ_NN, "ligand_target_matrix.rds"))
# saveRDS(lr_network, file = file.path(DIR_OBJ_NN, "lr_network.rds"))
# saveRDS(weighted_networks, file = file.path(DIR_OBJ_NN, "weighted_networks.rds"))

ligand_target_matrix <- readRDS(file = file.path(DIR_OBJ_NN, "ligand_target_matrix.rds"))
lr_network <- readRDS(file = file.path(DIR_OBJ_NN, "lr_network.rds"))
weighted_networks <- readRDS(file = file.path(DIR_OBJ_NN, "weighted_networks.rds"))

weighted_networks_lr <- weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))


##### Prep se obj ##### 
#' Define column containing hsNMF30-F14-hi-C0 and nb. cluster information
se.subset$f14_nbs_clusters2 <- as.character(se.subset$f14_nbs_clusters)
rename_spots_clusters <- grep("^[0-9]$", se.subset$f14_nbs_clusters, value=T)
se.subset$f14_nbs_clusters2[names(rename_spots_clusters)] <- paste0("F14C0_nbs_C", rename_spots_clusters)
se.subset$f14_nbs_clusters2 <- factor(se.subset$f14_nbs_clusters2, 
                                      levels = c(paste0("F14C0_nbs_C", 0:5), "F14_C0", "other"))
se.subset <- SetIdent(se.subset, value = "f14_nbs_clusters2")


#### NicheNet Round 1 ####
#' receiver: hsNMF-F14hi-C0 
#' senders: nb. clusters 0-5

##### Step 1. Define sender/receiver ##### 
## receiver
receiver = "F14_C0"
expressed_genes_receiver = get_expressed_genes(receiver, se.subset, pct = 0.10, assay_oi = "SCT")
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## sender
sender_celltypes = grep("F14C0_nbs_C", levels(se.subset$f14_nbs_clusters2), value = T) #sort(grep("nbs_cluster", unique(se.subset$f14_nbs_clusters2), value = T))
list_expressed_genes_sender = sender_celltypes %>% 
  unique() %>% 
  lapply(get_expressed_genes, se.subset, pct = 0.10, assay_oi = "SCT") # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()


##### Step 2. Define receiver gene set ##### 
#' 2. Define a gene set of interest
condition_oi = "F14_C0"
condition_reference = "other" 

DE_table_receiver = FindMarkers(object = se.subset, 
                                ident.1 = condition_oi, ident.2 = condition_reference, 
                                min.pct = 0.10) %>% 
  tibble::rownames_to_column("gene")

geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

#' Export DE table
write.csv(DE_table_receiver,
          file.path(DIR_RES, "objects", "A_NMF30", "hs_visium_A_nmf_1-30_f14high_C0_vs_other_NicheNet.DE_table_receiver.csv"), 
          row.names = F)
DE_table_receiver <- read.csv(file.path(DIR_RES, "objects", "A_NMF30", "hs_visium_A_nmf_1-30_f14high_C0_vs_other_NicheNet.DE_table_receiver.csv"))


##### Step 3. Define ligand gene set ##### 
#' 3. Define a set of potential ligands: these are ligands that are expressed 
#' by the “sender/niche” cell population and bind a (putative) receptor 
#' expressed by the “receiver/target” population
ligands = lr_network %>% 
  pull(from) %>% 
  unique()
receptors = lr_network %>% 
  pull(to) %>% 
  unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% 
  filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% 
  pull(from) %>% unique()


##### Step 4. Perform NicheNet ##### 
#' 4. Perform NicheNet ligand activity analysis: rank the potential ligands 
#' based on the presence of their target genes in the gene set of interest 
#' (compared to the background set of genes)
ligand_activities = predict_ligand_activities(geneset = geneset_oi, 
                                              background_expressed_genes = background_expressed_genes, 
                                              ligand_target_matrix = ligand_target_matrix, 
                                              potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% 
  arrange(-pearson) %>% 
  mutate(rank = rank(desc(pearson)))

#' Export results
write.csv(ligand_activities,
          file.path(DIR_RES, "objects", "A_NMF30", "hs_visium_A_nmf_1-30_f14high_C0_vs_nbs_NicheNet.ligand_activities.csv"), 
          row.names = F)
ligand_activities <- read.csv(file.path(DIR_RES, "objects", "A_NMF30", "hs_visium_A_nmf_1-30_f14high_C0_vs_nbs_NicheNet.ligand_activities.csv"))

#' View results
best_upstream_ligands = ligand_activities %>% top_n(40, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()

p <- DotPlot(se.subset, features = best_upstream_ligands, cols = "RdYlBu") + 
  RotatedAxis() + 
  theme_dotplot +
  scale_y_discrete(position = "right") +
  coord_flip() +
  theme(axis.text.x = element_text(angle=45, hjust=0));p


##### Step 5. Infer R-L activity ##### 
#' 5. Infer receptors and top-predicted target genes of ligands that are 
#' top-ranked in the ligand activity analysis

#' Active target gene inference.
#' Visualize regulatory potential in predicted upstream targets for our
#' best (pioritized) upstream ligands.
#' Adjust the cutoff to include more/less genes.
active_ligand_target_links_df = best_upstream_ligands %>% 
  lapply(get_weighted_ligand_target_links, geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% 
  bind_rows() %>% 
  tidyr::drop_na()
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, 
                                                                 ligand_target_matrix = ligand_target_matrix, 
                                                                 cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% 
  rev() %>% 
  make.names()
order_targets = active_ligand_target_links_df$target %>% 
  unique() %>% 
  intersect(rownames(active_ligand_target_links)) %>% 
  make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

#' Visualize results
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands", "Predicted target genes", 
                                                                    color = col_orange,
                                                                    legend_position = "bottom", 
                                                                    x_axis_position = "top",
                                                                    legend_title = "Regulatory potential") + 
  theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "grey80",  high = col_orange)#, breaks = c(0,0.0045,0.0090))
p_ligand_target_network

pdf(file = file.path(DIR_NN, "hs_visium_A_NN_F14-C0_vs_nbs_upsteam_ligand_target_hm.pdf"), width = 15, height = 5);p_ligand_target_network;dev.off()

#' Export results
write.csv(active_ligand_target_links_df,
          file.path(DIR_RES, "objects", "A_NMF30", "hs_visium_A_nmf_1-30_f14high_C0_vs_nbs_NicheNet.active_ligand_target_links.csv"), 
          row.names = F)


##### Step 6. Include LogFC ##### 
#' 6. Add log fold change information of ligands from sender cells
# DE analysis for each sender cell type
# this uses a new nichenetr function - reinstall nichenetr if necessary!
se.subset$f14_nbs_clusters2 <- factor(se.subset$f14_nbs_clusters2, levels = sort(unique(se.subset$f14_nbs_clusters2)))
se.subset <- SetIdent(se.subset, value = "f14_nbs_clusters2")

#' Find markers for each cluster (set against all the rest)
DE_table_all <- FindAllMarkers(se.subset,
                               min.pct = 0.25,
                               max.cells.per.ident = 2e3)
# colnames(DE_table_all) <- c("gene", levels(se.subset$f14_nbs_clusters2))

DE_table_all <- DE_table_all %>% tibble() %>% filter(p_val_adj < 0.01) %>% select(avg_log2FC, cluster, gene)
DE_table_all <- spread(DE_table_all, key = cluster, value = avg_log2FC)
DE_table_all[is.na(DE_table_all)] = 0

#' Export DE table
write.csv(DE_table_all,
          file.path(DIR_RES, "objects", "A_NMF30", "hs_visium_A_nmf_1-30_f14high_C0_nbs_clusters_NicheNet.DE_table_all.csv"),
          row.names = F)
DE_table_all <- read.csv(file.path(DIR_RES, "objects", "A_NMF30", "hs_visium_A_nmf_1-30_f14high_C0_nbs_clusters_NicheNet.DE_table_all.csv"))

#' Combine ligand activities with DE information
ligand_activities_de = ligand_activities %>% select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE_table_all %>% rename(ligand = gene))
ligand_activities_de[is.na(ligand_activities_de)] = 0

#' Make LFC heatmap (C)
lfc_matrix = ligand_activities_de %>% select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
vis_ligand_lfc = lfc_matrix[order_ligands,]

colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()

p_ligand_lfc = vis_ligand_lfc %>% 
  make_threecolor_heatmap_ggplot("Prioritized ligands",
                                 "Log2 Fold-Change", 
                                 low_color = col_blue,
                                 mid_color = "grey95", 
                                 mid = median(vis_ligand_lfc), 
                                 high_color = col_orange,
                                 legend_position = "right", 
                                 x_axis_position = "top", 
                                 legend_title = "Log2FC")
p_ligand_lfc <- p_ligand_lfc + theme(axis.text.y = element_text(face = "italic", size = 9, colour = "grey30"));p_ligand_lfc

pdf(file = file.path(DIR_NN, "hs_visium_A_NN_F14-C0_vs_nbs_ligand_logfc_hm.pdf"), width = 3.5, height = 7);p_ligand_lfc;dev.off()

##### Step 7. Visualize all ##### 
#' Summarize visualizations of results for
#' Prioritized ligands and their exåression in sender spots

#' A. Ligand activity heatmap
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% 
  make_heatmap_ggplot("Prioritized ligands", 
                      "Ligand activity", 
                      color = col_orange,
                      legend_position = "left", 
                      x_axis_position = "top", 
                      # x_axis = F,
                      legend_title = "Pearson correlation \ncoefficient (target gene\nprediction ability") + 
  theme(legend.text = element_text(size = 10),
        axis.title.y = element_text(size = 10, colour = "black"))
p_ligand_pearson


#' B. Ligand expression Seurat dotplot
order_ligands_adapted = order_ligands
rotated_dotplot = DotPlot(se.subset, features = order_ligands_adapted, cols = "RdYlBu") + 
  theme_dotplot +
  scale_y_discrete(position = "right") +
  coord_flip() +
  ylab("Expression") +
  theme(legend.position = "bottom",
        legend.box="vertical", 
        legend.key.size = unit(10, "pt"),
        legend.text = element_text(size = 9, colour = "grey30"),
        legend.title = element_text(size = 9, colour = "grey30"),
        # panel.background = element_rect(fill = NA, color = "grey30", size = .5),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=0, vjust=1, size = 9, colour = "grey30"),
        axis.title.x = element_text(size = 10, colour = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic", size = 9, colour = "grey30"))
rotated_dotplot

#' A + B + C (p_ligand_lfc)
p_ligand_lfc2 <- p_ligand_lfc + theme(axis.title.y = element_blank())
p_patch <- p_ligand_pearson + rotated_dotplot + p_ligand_lfc2 + patchwork::plot_layout(nrow = 1, widths = c(1, 4, 4))
p_patch <- p_patch & theme(axis.title.x = element_text(size = 10, colour = "black"), 
                           axis.ticks = element_blank(),
                           legend.text = element_text(size = 9, colour = "grey30"),
                           legend.title = element_text(size = 9, colour = "grey30"));p_patch

pdf(file = file.path(DIR_NN, "hs_visium_A_NN_F14-C0_vs_nbs_grid.pdf"), width = 8, height = 8);p_patch;dev.off()


#### NicheNet Round 2 ####
#' receiver: nb. cluster 0 ("F14C0_nbs_C0") (AT cells)
#' sender: hsNMF-F14hi-C0

fname_nn_comparison <- "senderF14-C0_receiverF14-C0-Nbs-C0"
se.subset <- SetIdent(se.subset, value = "f14_nbs_clusters2")

##### Step 1 ##### 
## receiver
receiver = "F14C0_nbs_C0"
expressed_genes_receiver = get_expressed_genes(receiver, se.subset, pct = 0.10, assay_oi = "SCT")
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## sender
sender_celltypes = "F14_C0"
list_expressed_genes_sender = sender_celltypes %>% 
  unique() %>% 
  lapply(get_expressed_genes, se.subset, pct = 0.10, assay_oi = "SCT") # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()


##### Step 2-3 ##### 
#' Define genes: receivers
condition_oi = sender_celltypes
condition_reference = "other" 
DE_table_receiver = FindMarkers(object = se.subset, 
                                ident.1 = condition_oi, 
                                ident.2 = condition_reference, 
                                min.pct = 0.10, 
                                max.cells.per.ident = 5e3) %>% 
  tibble::rownames_to_column("gene")

geneset_oi = DE_table_receiver %>% 
  filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% 
  pull(gene) %>% 
  .[. %in% rownames(ligand_target_matrix)]

#' Define genes: ligands
ligands = lr_network %>% 
  pull(from) %>% 
  unique()
receptors = lr_network %>% 
  pull(to) %>% 
  unique()

expressed_ligands = intersect(ligands, expressed_genes_sender)
expressed_receptors = intersect(receptors, expressed_genes_receiver)

potential_ligands = lr_network %>% 
  filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% 
  pull(from) %>% 
  unique()

##### Step 4. Perform NicheNet ##### 
#' Perform NicheNet
ligand_activities = predict_ligand_activities(geneset = geneset_oi, 
                                              background_expressed_genes = background_expressed_genes, 
                                              ligand_target_matrix = ligand_target_matrix, 
                                              potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% 
  arrange(-pearson) %>% 
  mutate(rank = rank(desc(pearson)))

# Export results
write.csv(ligand_activities,
          file.path(DIR_RES, "objects", "A_NMF30", 
                    paste0("hs_visium_A_nmf_1-30_NicheNet_", fname_nn_comparison, ".ligand_activities.csv")), 
          row.names = F)

ligand_activities <- read.csv(file.path(DIR_RES, "objects", "A_NMF30", 
                                        paste0("hs_visium_A_nmf_1-30_NicheNet_", fname_nn_comparison, ".ligand_activities.csv"))
                              )

# Look at results
best_upstream_ligands = ligand_activities %>% 
  top_n(40, pearson) %>% 
  arrange(-pearson) %>% 
  pull(test_ligand) %>% 
  unique()


##### Step 5-6 ##### 
# Infer R-L activity
# Active target gene inference
active_ligand_target_links_df = best_upstream_ligands %>% 
  lapply(get_weighted_ligand_target_links, geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% 
  bind_rows() %>% 
  tidyr::drop_na()
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, 
                                                                 ligand_target_matrix = ligand_target_matrix, 
                                                                 cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% 
  rev() %>% 
  make.names()
order_targets = active_ligand_target_links_df$target %>% 
  unique() %>% 
  intersect(rownames(active_ligand_target_links)) %>% 
  make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

#' Export results
write.csv(active_ligand_target_links_df,
          file.path(DIR_RES, "objects", "A_NMF30", 
                    paste0("hs_visium_A_nmf_1-30_NicheNet_", fname_nn_comparison, ".active_ligand_target_links.csv")),
          row.names = F)

active_ligand_target_links_df <- read.csv(file.path(DIR_RES, "objects", "A_NMF30", 
                                                    paste0("hs_visium_A_nmf_1-30_NicheNet_", fname_nn_comparison, ".active_ligand_target_links.csv"))
                                          )

#' Include LogFC
DE_table_all <- read.csv(file.path(DIR_RES, "objects", "A_NMF30", 
                                   "hs_visium_A_nmf_1-30_f14high_C0_nbs_clusters_NicheNet.DE_table_all.csv"))

# Combine ligand activities with DE information
ligand_activities_de = ligand_activities %>% 
  select(test_ligand, pearson) %>% 
  rename(ligand = test_ligand) %>% 
  left_join(DE_table_all %>% 
              rename(ligand = gene))
ligand_activities_de[is.na(ligand_activities_de)] = 0


##### Step 7. Visualize all ##### 
#' Summarize visualizations of results
# A: ligand activity heatmap
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% 
  make_heatmap_ggplot("Prioritized ligands", 
                      "Ligand activity", 
                      color = col_orange,
                      legend_position = "left", 
                      x_axis_position = "top", 
                      # x_axis = F,
                      legend_title = "Pearson correlation \ncoefficient (target gene\nprediction ability") + 
  theme(legend.text = element_text(size = 10),
        axis.title.y = element_text(size = 10, colour = "black"))

# B: ligand expression Seurat dotplot
order_ligands_adapted = order_ligands
rotated_dotplot = DotPlot(se.subset, features = order_ligands_adapted, cols = "RdYlBu") + 
  theme_dotplot +
  scale_y_discrete(position = "right") +
  coord_flip() +
  ylab("Expression") +
  theme(legend.position = "bottom",
        legend.box="vertical", 
        legend.key.size = unit(10, "pt"),
        legend.text = element_text(size = 9, colour = "grey30"),
        legend.title = element_text(size = 9, colour = "grey30"),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=0, vjust=1, size = 9, colour = "grey30"),
        axis.title.x = element_text(size = 10, colour = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic", size = 9, colour = "grey30"))

# C: LFC heatmap
lfc_matrix = ligand_activities_de %>% 
  select(-ligand, -pearson) %>% 
  as.matrix() %>% 
  magrittr::set_rownames(ligand_activities_de$ligand)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
vis_ligand_lfc = lfc_matrix[order_ligands,]

colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()

p_ligand_lfc = vis_ligand_lfc %>% 
  make_threecolor_heatmap_ggplot("Prioritized ligands",
                                 "Log2 Fold-Change", 
                                 low_color = col_blue,
                                 mid_color = "grey95", 
                                 mid = median(vis_ligand_lfc), 
                                 high_color = col_orange,
                                 legend_position = "right", 
                                 x_axis_position = "top", 
                                 legend_title = "Log2FC")
p_ligand_lfc <- p_ligand_lfc + theme(axis.text.y = element_text(face = "italic", size = 9, colour = "grey30"))
p_ligand_lfc2 <- p_ligand_lfc + theme(axis.title.y = element_blank())

# Patch together and save
p_patch <- p_ligand_pearson + rotated_dotplot + p_ligand_lfc2 + patchwork::plot_layout(nrow = 1, widths = c(1, 4, 4))
p_patch <- p_patch & theme(axis.title.x = element_text(size = 10, colour = "black"), 
                           axis.ticks = element_blank(),
                           legend.text = element_text(size = 9, colour = "grey30"),
                           legend.title = element_text(size = 9, colour = "grey30"))

pdf(file = file.path(DIR_NN, paste0("hs_visium_A_nmf_1-30_NicheNet_", fname_nn_comparison, ".plotgrid.pdf")), 
    width = 8, height = 8)
p_patch + patchwork::plot_annotation(title = fname_nn_comparison)
dev.off()


#### NicheNet Round 3 ####
#' receiver: nb. cluster 1 ("F14C0_nbs_C1") (Fib cells)
#' sender: hsNMF-F14hi-C0

fname_nn_comparison <- "senderF14-C0_receiverF14-C0-Nbs-C1"
se.subset <- SetIdent(se.subset, value = "f14_nbs_clusters2")

##### Step 1 ##### 
## receiver
receiver = "F14C0_nbs_C1"
expressed_genes_receiver = get_expressed_genes(receiver, se.subset, pct = 0.10, assay_oi = "SCT")
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## sender
sender_celltypes = "F14_C0"
list_expressed_genes_sender = sender_celltypes %>% 
  unique() %>% 
  lapply(get_expressed_genes, se.subset, pct = 0.10, assay_oi = "SCT") # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

##### Step 2-3 ##### 
#' Define genes: receivers
condition_oi = sender_celltypes
condition_reference = "other" 
DE_table_receiver = FindMarkers(object = se.subset, 
                                ident.1 = condition_oi, 
                                ident.2 = condition_reference, 
                                min.pct = 0.10, 
                                max.cells.per.ident = 5e3) %>% 
  tibble::rownames_to_column("gene")

geneset_oi = DE_table_receiver %>% 
  filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% 
  pull(gene) %>% 
  .[. %in% rownames(ligand_target_matrix)]

#' Define genes: ligands
ligands = lr_network %>% 
  pull(from) %>% 
  unique()
receptors = lr_network %>% 
  pull(to) %>% 
  unique()

expressed_ligands = intersect(ligands, expressed_genes_sender)
expressed_receptors = intersect(receptors, expressed_genes_receiver)

potential_ligands = lr_network %>% 
  filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% 
  pull(from) %>% 
  unique()


##### Step 4. Perform NicheNet ##### 
#' Perform NicheNet
ligand_activities = predict_ligand_activities(geneset = geneset_oi, 
                                              background_expressed_genes = background_expressed_genes, 
                                              ligand_target_matrix = ligand_target_matrix, 
                                              potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% 
  arrange(-pearson) %>% 
  mutate(rank = rank(desc(pearson)))

# Export results
write.csv(ligand_activities,
          file.path(DIR_RES, "objects", "A_NMF30", 
                    paste0("hs_visium_A_nmf_1-30_NicheNet_", fname_nn_comparison, ".ligand_activities.csv")), 
          row.names = F)

ligand_activities <- read.csv(file.path(DIR_RES, "objects", "A_NMF30", 
                                        paste0("hs_visium_A_nmf_1-30_NicheNet_", fname_nn_comparison, ".ligand_activities.csv"))
                              )

# Look at results
best_upstream_ligands = ligand_activities %>% 
  top_n(40, pearson) %>% 
  arrange(-pearson) %>% 
  pull(test_ligand) %>% 
  unique()
best_upstream_ligands


##### Step 5-6 ##### 
# Infer R-L activity
# Active target gene inference
active_ligand_target_links_df = best_upstream_ligands %>% 
  lapply(get_weighted_ligand_target_links, geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% 
  bind_rows() %>% 
  tidyr::drop_na()
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, 
                                                                 ligand_target_matrix = ligand_target_matrix, 
                                                                 cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% 
  rev() %>% 
  make.names()
order_targets = active_ligand_target_links_df$target %>% 
  unique() %>% 
  intersect(rownames(active_ligand_target_links)) %>% 
  make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

#' Export results
write.csv(active_ligand_target_links_df,
          file.path(DIR_RES, "objects", "A_NMF30", 
                    paste0("hs_visium_A_nmf_1-30_NicheNet_", fname_nn_comparison, ".active_ligand_target_links.csv")),
          row.names = F)

active_ligand_target_links_df <- read.csv(file.path(DIR_RES, "objects", "A_NMF30", 
                                                    paste0("hs_visium_A_nmf_1-30_NicheNet_", fname_nn_comparison, ".active_ligand_target_links.csv"))
                                          )

#' Include LogFC
DE_table_all <- read.csv(file.path(DIR_RES, "objects", "A_NMF30", 
                                   "hs_visium_A_nmf_1-30_f14high_C0_nbs_clusters_NicheNet.DE_table_all.csv"))

#' Combine ligand activities with DE information
ligand_activities_de = ligand_activities %>% 
  select(test_ligand, pearson) %>% 
  rename(ligand = test_ligand) %>% 
  left_join(DE_table_all %>% 
              rename(ligand = gene))
ligand_activities_de[is.na(ligand_activities_de)] = 0


##### Step 7. Visualize all ##### 
#' Summarize visualizations of results
# A: ligand activity heatmap
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% 
  make_heatmap_ggplot("Prioritized ligands", 
                      "Ligand activity", 
                      color = col_orange,
                      legend_position = "left", 
                      x_axis_position = "top", 
                      # x_axis = F,
                      legend_title = "Pearson correlation \ncoefficient (target gene\nprediction ability") + 
  theme(legend.text = element_text(size = 10),
        axis.title.y = element_text(size = 10, colour = "black"))

# B: ligand expression Seurat dotplot
order_ligands_adapted = order_ligands
rotated_dotplot = DotPlot(se.subset, features = order_ligands_adapted, cols = "RdYlBu") + 
  theme_dotplot +
  scale_y_discrete(position = "right") +
  coord_flip() +
  ylab("Expression") +
  theme(legend.position = "bottom",
        legend.box="vertical", 
        legend.key.size = unit(10, "pt"),
        legend.text = element_text(size = 9, colour = "grey30"),
        legend.title = element_text(size = 9, colour = "grey30"),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=0, vjust=1, size = 9, colour = "grey30"),
        axis.title.x = element_text(size = 10, colour = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic", size = 9, colour = "grey30"))

# C: LFC heatmap
lfc_matrix = ligand_activities_de %>% 
  select(-ligand, -pearson) %>% 
  as.matrix() %>% 
  magrittr::set_rownames(ligand_activities_de$ligand)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
vis_ligand_lfc = lfc_matrix[order_ligands,]

colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()

p_ligand_lfc = vis_ligand_lfc %>% 
  make_threecolor_heatmap_ggplot("Prioritized ligands",
                                 "Log2 Fold-Change", 
                                 low_color = col_blue,
                                 mid_color = "grey95", 
                                 mid = median(vis_ligand_lfc), 
                                 high_color = col_orange,
                                 legend_position = "right", 
                                 x_axis_position = "top", 
                                 legend_title = "Log2FC")
p_ligand_lfc <- p_ligand_lfc + theme(axis.text.y = element_text(face = "italic", size = 9, colour = "grey30"))
p_ligand_lfc2 <- p_ligand_lfc + theme(axis.title.y = element_blank())

# Patch together and save
p_patch <- p_ligand_pearson + rotated_dotplot + p_ligand_lfc2 + patchwork::plot_layout(nrow = 1, widths = c(1, 4, 4))
p_patch <- p_patch & theme(axis.title.x = element_text(size = 10, colour = "black"), 
                           axis.ticks = element_blank(),
                           legend.text = element_text(size = 9, colour = "grey30"),
                           legend.title = element_text(size = 9, colour = "grey30"))

pdf(file = file.path(DIR_NN, paste0("hs_visium_A_nmf_1-30_NicheNet_", fname_nn_comparison, ".plotgrid.pdf")), 
    width = 8, height = 8)
p_patch + patchwork::plot_annotation(title = fname_nn_comparison)
dev.off()


#### NicheNet Round 4 ####
#' receiver: nb. cluster 2 ("F14C0_nbs_C2") (Mac cells)
#' sender: hsNMF-F14hi-C0

fname_nn_comparison <- "senderF14-C0_receiverF14-C0-Nbs-C2"
se.subset <- SetIdent(se.subset, value = "f14_nbs_clusters2")

##### Step 1 ##### 
## receiver
receiver = "F14C0_nbs_C2"
expressed_genes_receiver = get_expressed_genes(receiver, se.subset, pct = 0.10, assay_oi = "SCT")
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## sender
sender_celltypes = "F14_C0"
list_expressed_genes_sender = sender_celltypes %>% 
  unique() %>% 
  lapply(get_expressed_genes, se.subset, pct = 0.10, assay_oi = "SCT") # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

##### Step 2-3 ##### 
#' Define genes: receivers
condition_oi = sender_celltypes
condition_reference = "other" 
DE_table_receiver = FindMarkers(object = se.subset, 
                                ident.1 = condition_oi, 
                                ident.2 = condition_reference, 
                                min.pct = 0.10, 
                                max.cells.per.ident = 5e3) %>% 
  tibble::rownames_to_column("gene")

geneset_oi = DE_table_receiver %>% 
  filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% 
  pull(gene) %>% 
  .[. %in% rownames(ligand_target_matrix)]

#' Define genes: ligands
ligands = lr_network %>% 
  pull(from) %>% 
  unique()
receptors = lr_network %>% 
  pull(to) %>% 
  unique()

expressed_ligands = intersect(ligands, expressed_genes_sender)
expressed_receptors = intersect(receptors, expressed_genes_receiver)

potential_ligands = lr_network %>% 
  filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% 
  pull(from) %>% 
  unique()


##### Step 4. Perform NicheNet ##### 
#' Perform NicheNet
ligand_activities = predict_ligand_activities(geneset = geneset_oi, 
                                              background_expressed_genes = background_expressed_genes, 
                                              ligand_target_matrix = ligand_target_matrix, 
                                              potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% 
  arrange(-pearson) %>% 
  mutate(rank = rank(desc(pearson)))

# Save results
write.csv(ligand_activities,
          file.path(DIR_RES, "objects", "A_NMF30", 
                    paste0("hs_visium_A_nmf_1-30_NicheNet_", fname_nn_comparison, ".ligand_activities.csv")), 
          row.names = F)

ligand_activities <- read.csv(file.path(DIR_RES, "objects", "A_NMF30", 
                                        paste0("hs_visium_A_nmf_1-30_NicheNet_", fname_nn_comparison, ".ligand_activities.csv")))

# Look at results
best_upstream_ligands = ligand_activities %>% 
  top_n(40, pearson) %>% 
  arrange(-pearson) %>% 
  pull(test_ligand) %>% 
  unique()
best_upstream_ligands


# p <- DotPlot(se.subset, features = best_upstream_ligands, cols = "RdYlBu") + 
#   RotatedAxis() + 
#   theme_dotplot +
#   scale_y_discrete(position = "right") +
#   coord_flip() +
#   labs(title=fname_nn_comparison) +
#   theme(axis.text.x = element_text(angle=45, hjust=0));p
# pdf(file = file.path(DIR_NN, "hs_visium_A_NN_F14-C0_vs_nbs_upsteam_ligands_dotplot.pdf"), width = 4, height = 8);p;dev.off()


##### Step 5-6 ##### 
# Infer R-L activity
# Active target gene inference
active_ligand_target_links_df = best_upstream_ligands %>% 
  lapply(get_weighted_ligand_target_links, geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% 
  bind_rows() %>% 
  tidyr::drop_na()
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, 
                                                                 ligand_target_matrix = ligand_target_matrix, 
                                                                 cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% 
  rev() %>% 
  make.names()
order_targets = active_ligand_target_links_df$target %>% 
  unique() %>% 
  intersect(rownames(active_ligand_target_links)) %>% 
  make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

#' Export results
write.csv(active_ligand_target_links_df,
          file.path(DIR_RES, "objects", "A_NMF30", 
                    paste0("hs_visium_A_nmf_1-30_NicheNet_", fname_nn_comparison, ".active_ligand_target_links.csv")),
          row.names = F)
active_ligand_target_links_df <- read.csv(file.path(DIR_RES, "objects", "A_NMF30", 
                                                    paste0("hs_visium_A_nmf_1-30_NicheNet_", fname_nn_comparison, ".active_ligand_target_links.csv")))


#' Include LogFC
DE_table_all <- read.csv(file.path(DIR_RES, "objects", "A_NMF30", 
                                   "hs_visium_A_nmf_1-30_f14high_C0_nbs_clusters_NicheNet.DE_table_all.csv"))

# Combine ligand activities with DE information
ligand_activities_de = ligand_activities %>% 
  select(test_ligand, pearson) %>% 
  rename(ligand = test_ligand) %>% 
  left_join(DE_table_all %>% 
              rename(ligand = gene))
ligand_activities_de[is.na(ligand_activities_de)] = 0


##### Step 7. Visualize all ##### 
#' Summarize visualizations of results
# A: ligand activity heatmap
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% 
  make_heatmap_ggplot("Prioritized ligands", 
                      "Ligand activity", 
                      color = col_orange,
                      legend_position = "left", 
                      x_axis_position = "top", 
                      # x_axis = F,
                      legend_title = "Pearson correlation \ncoefficient (target gene\nprediction ability") + 
  theme(legend.text = element_text(size = 10),
        axis.title.y = element_text(size = 10, colour = "black"))

# B: ligand expression Seurat dotplot
order_ligands_adapted = order_ligands
rotated_dotplot = DotPlot(se.subset, features = order_ligands_adapted, cols = "RdYlBu") + 
  theme_dotplot +
  scale_y_discrete(position = "right") +
  coord_flip() +
  ylab("Expression") +
  theme(legend.position = "bottom",
        legend.box="vertical", 
        legend.key.size = unit(10, "pt"),
        legend.text = element_text(size = 9, colour = "grey30"),
        legend.title = element_text(size = 9, colour = "grey30"),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=0, vjust=1, size = 9, colour = "grey30"),
        axis.title.x = element_text(size = 10, colour = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic", size = 9, colour = "grey30"))

# C: LFC heatmap
lfc_matrix = ligand_activities_de %>% 
  select(-ligand, -pearson) %>% 
  as.matrix() %>% 
  magrittr::set_rownames(ligand_activities_de$ligand)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
vis_ligand_lfc = lfc_matrix[order_ligands,]

colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()

p_ligand_lfc = vis_ligand_lfc %>% 
  make_threecolor_heatmap_ggplot("Prioritized ligands",
                                 "Log2 Fold-Change", 
                                 low_color = col_blue,
                                 mid_color = "grey95", 
                                 mid = median(vis_ligand_lfc), 
                                 high_color = col_orange,
                                 legend_position = "right", 
                                 x_axis_position = "top", 
                                 legend_title = "Log2FC")
p_ligand_lfc <- p_ligand_lfc + theme(axis.text.y = element_text(face = "italic", size = 9, colour = "grey30"))
p_ligand_lfc2 <- p_ligand_lfc + theme(axis.title.y = element_blank())

# Patch together and save
p_patch <- p_ligand_pearson + rotated_dotplot + p_ligand_lfc2 + patchwork::plot_layout(nrow = 1, widths = c(1, 4, 4))
p_patch <- p_patch & theme(axis.title.x = element_text(size = 10, colour = "black"), 
                           axis.ticks = element_blank(),
                           legend.text = element_text(size = 9, colour = "grey30"),
                           legend.title = element_text(size = 9, colour = "grey30"))

pdf(file = file.path(DIR_NN, paste0("hs_visium_A_nmf_1-30_NicheNet_", fname_nn_comparison, ".plotgrid.pdf")), 
    width = 8, height = 8)
p_patch + patchwork::plot_annotation(title = fname_nn_comparison)
dev.off()


#### Combine results ####
#' Comparison: F14-C0 (AbBa) vs F14-C0-Nbs C0-C2
#' 
#' Code used to generate figures shown in article

#' Fetch results
ligand_act_list <- list(
  "recF14C0_sendF14C0nbs" = read.csv(file.path(DIR_RES, "objects", "A_NMF30", "hs_visium_A_nmf_1-30_f14high_C0_vs_nbs_NicheNet.ligand_activities.csv")),
  "recF14C0nbsC0_sendF14C0" = read.csv(file.path(DIR_RES, "objects", "A_NMF30", 
                                                 paste0("hs_visium_A_nmf_1-30_NicheNet_", "senderF14-C0_receiverF14-C0-Nbs-C0", ".ligand_activities.csv"))
  ),
  "recF14C0nbsC1_sendF14C0" = read.csv(file.path(DIR_RES, "objects", "A_NMF30", 
                                                 paste0("hs_visium_A_nmf_1-30_NicheNet_", "senderF14-C0_receiverF14-C0-Nbs-C1", ".ligand_activities.csv"))
  ),
  "recF14C0nbsC2_sendF14C0" = read.csv(file.path(DIR_RES, "objects", "A_NMF30", 
                                                 paste0("hs_visium_A_nmf_1-30_NicheNet_", "senderF14-C0_receiverF14-C0-Nbs-C2", ".ligand_activities.csv"))
  )
)

ligand_act_list <- lapply(names(ligand_act_list), function(nn_comp){
  d <- ligand_act_list[[nn_comp]]
  comp_names <- strsplit(nn_comp, split = "_") %>% unlist()
  d$nn_comp <- nn_comp
  d$receiver <- grep("rec", comp_names, value = T) %>% gsub("rec", "", x=.)
  d$sender <- grep("send", comp_names, value = T) %>% gsub("send", "", x=.)
  return(d)
})

ligand_act_df <- bind_rows(ligand_act_list)

#' Export
write.csv(ligand_act_df, file.path(DIR_NN, "hs_visium_A_nmf_1-30_NicheNet_ligand_comparison_pearson.heatmap_unfiltered.csv"))

#' Filter
ligand_act_df <- ligand_act_df %>% filter(pearson>0)
head(ligand_act_df)

##### Plot Pearson cor ##### 
#' Unfiltered Pearson cor results
pdf(file = file.path(DIR_NN, paste0("hs_visium_A_nmf_1-30_NicheNet_ligand_comparison_pearson.heatmap.pdf")), 
    width = 3, height = 22)
ggplot(ligand_act_df, aes(y=reorder(test_ligand, desc(test_ligand)), x=nn_comp, fill=pearson)) + 
  geom_tile(color="white", width=0.95) +
  scale_fill_distiller(palette = "RdPu", direction = 1) +
  labs(x="Comparison", "Predicted ligand") +
  scale_x_discrete(position = "top") +
  theme_dotplot +
  theme(
    axis.text.y = element_text(face = "plain", size=7),
    axis.text.x = element_text(angle=45, hjust=0, face="plain", size=10),
    panel.grid = element_blank(),
    legend.position = "right")
dev.off()


#' Filtered Pearson cor results
ligands_filt <- ligand_act_df %>% 
  group_by(test_ligand) %>% 
  summarise(mean_p = mean(pearson)) %>%
  arrange(desc(mean_p)) %>%
  filter(mean_p>0.075)

ligand_act_df_filt <- ligand_act_df %>% filter(test_ligand %in% ligands_filt$test_ligand)

ligand_act_df_filt$nn_comp <- ligand_act_df_filt$nn_comp %>% 
  as.factor() %>% 
  recode_factor(
    recF14C0_sendF14C0nbs = "F14.C0.nbs_to_F14.C0",
    recF14C0nbsC0_sendF14C0 = "F14.C0_to_nbC0", 
    recF14C0nbsC1_sendF14C0 = "F14.C0_to_nbC1",
    recF14C0nbsC2_sendF14C0 = "F14.C0_to_nbC2"
  )

ligand_act_df_filt$test_ligand <- factor(ligand_act_df_filt$test_ligand,
                                         levels = arrange(ligand_act_df_filt, desc(pearson))$test_ligand %>% unique() %>% rev())

pdf(file = file.path(DIR_NN, paste0("hs_visium_A_nmf_1-30_NicheNet_ligand_comparison_pearson.heatmap_filtered2.pdf")), 
    width = 2.5, height = 4.5)
ggplot(ligand_act_df_filt, aes(y=test_ligand, x=nn_comp, fill=pearson)) +  # reorder(test_ligand, pearson)
  geom_tile(width=0.9, height=0.9) +
  scale_fill_scico(palette = "acton", direction = -1) +
  labs(x="Comparison", "Predicted ligand") +
  scale_x_discrete(position = "top") +
  theme_dotplot +
  theme(
    axis.text.y = element_text(face = "plain", size=8),
    axis.text.x = element_text(angle=45, hjust=0, face="plain", size=8),
    panel.grid = element_blank(),
    legend.text = element_text(size=8),
    legend.position = "right")
dev.off()


##### Plot gene expression ##### 
#' Gene expression dotplot for filtered ligands
genes_plot <- levels(ligand_act_df_filt$test_ligand)

pdf(file = file.path(DIR_NN, paste0("hs_visium_A_nmf_1-30_NicheNet_ligand_comparison_geneexpr.heatmap_filtered.pdf")), 
    width = 3.5, height = 4.5)
DotPlot(se.subset, features = genes_plot) +  # , cols = "RdYlBu"
  scale_color_gradientn(colours = col_scale_div_expr, limits = c(-2.5, 2.5)) +
  theme_dotplot +
  scale_y_discrete(position = "right") +
  coord_flip() +
  theme(legend.position = "right",
        legend.key.size = unit(10, "pt"),
        legend.text = element_text(size = 8, colour = "black"),
        legend.title = element_text(size = 8, colour = "black"),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=0, size = 8, colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic", size = 8, colour = "black"))
dev.off()


##### Plot gene Log2FC -  not included in article ##### 
#' Log2FC heatmap for filtered ligands
DE_table_all <- read.csv(file.path(DIR_RES, "objects", "A_NMF30", 
                                   "hs_visium_A_nmf_1-30_f14high_C0_nbs_clusters_NicheNet.DE_table_all.csv"))

genes_plot <- levels(ligand_act_df_filt$test_ligand)
DE_table_all[DE_table_all$gene %in% genes_plot,]

# Combine ligand activities with DE information
ligand_activities_de = ligand_act_df_filt %>% 
  select(test_ligand, pearson) %>% 
  rename(ligand = test_ligand) %>% 
  left_join(DE_table_all %>% 
              rename(ligand = gene))
ligand_activities_de[is.na(ligand_activities_de)] = 0

lfc_matrix = ligand_activities_de %>% 
  select(-ligand, -pearson) %>% 
  as.matrix() %>% 
  magrittr::set_rownames(ligand_activities_de$ligand)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
vis_ligand_lfc = lfc_matrix[order_ligands,]

colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()

p_ligand_lfc = vis_ligand_lfc %>% 
  make_threecolor_heatmap_ggplot("Prioritized ligands",
                                 "Log2 Fold-Change", 
                                 low_color = col_blue,
                                 mid_color = "grey95", 
                                 mid = median(vis_ligand_lfc), 
                                 high_color = col_orange,
                                 legend_position = "right", 
                                 x_axis_position = "top", 
                                 legend_title = "Log2FC")
p_ligand_lfc <- p_ligand_lfc + theme(axis.text.y = element_text(face = "italic", size = 9, colour = "grey30"))
p_ligand_lfc2 <- p_ligand_lfc + theme(axis.title.y = element_blank())


#### Other ####
##### Selected ligand R-L gene correlations ##### 
se_d_subset <- SubsetSTData(se.subset, d_c0 != "NA")

se_ipf_bg <- SubsetSTData(se.subset, 
                          spots = se.subset[[]] %>% filter(condition == "IPF" & is.na(d_c0)) %>% rownames())
se_ipf_bg$d_c0 <- 4

# Define gene of interest
gene_oi <- "APOE"
# gene_oi <- "IL1B"
# gene_oi <- "TGFB1"
# gene_oi <- "SPP1"

gene_oi_rl <- lr_network %>% filter(from==gene_oi) %>% select(to) %>% unique() %>% c(., gene_oi) %>% unlist() %>% as.character()

gene_data <- se_d_subset[[]] %>% 
  bind_cols(FetchData(se_d_subset, vars = gene_oi_rl))
gene_oi_rl <- intersect(gene_oi_rl, colnames(gene_data))

bg_data <- se_ipf_bg[[]] %>% 
  bind_cols(FetchData(se_ipf_bg, vars = gene_oi_rl)) %>% 
  select(matches(c("sample_name", "d_c0", gene_oi_rl))) %>% 
  pivot_longer(all_of(gene_oi_rl), 
               names_to = "variable", 
               values_to = "value") %>% 
  group_by(variable, d_c0) %>% 
  summarise(expr_sum = sum(value),
            expr_avg = mean(value))

gene_cor <- setNames(lapply(gene_oi_rl, function(g){
  cor(x = gene_data[["d_c0"]], y = gene_data[[g]])
}), nm = gene_oi_rl)

gene_cor_df <- data.frame(unlist(gene_cor) %>% sort())
colnames(gene_cor_df) <- "pearson_cor"

gene_plot_order <- gene_cor_df %>% arrange(-abs(pearson_cor)) %>% rownames()
gene_plot_order <- c(gene_oi, gene_plot_order[-match(gene_oi, gene_plot_order)])
gene_cor_df$variable <- rownames(gene_cor_df) %>% factor(levels = gene_plot_order)

p_gene_cor_data <- gene_data %>% 
  pivot_longer(all_of(gene_oi_rl), names_to = "variable", values_to = "value")
p_gene_cor_data$variable <- factor(p_gene_cor_data$variable, levels = gene_plot_order)

# heatmap
p_gene_cor_data_hm <- p_gene_cor_data %>% 
  group_by(variable, d_c0) %>% 
  summarise(expr_sum = sum(value),
            expr_avg = mean(value))

d_plot <- bind_rows(p_gene_cor_data_hm, bg_data)
d_plot$d_c0 <- recode_factor(factor(d_plot$d_c0), `4` = "n") %>% factor(levels = c(0:3, "n"))
d_plot$variable <- factor(d_plot$variable, levels=gene_plot_order)

p1 <- ggplot(d_plot, aes(y=reorder(variable, desc(variable)), x=d_c0, fill=expr_avg, size=expr_avg)) + 
  geom_vline(xintercept = "0", color = "#DBC8DB", linewidth = 5) +
  geom_vline(xintercept = "n", color = "grey90", linewidth = 5) +
  # geom_tile(color=NA, width=1, height=1) +
  geom_point(shape=21) +
  # geom_vline(xintercept = 0.5, color = "black", linetype = "dashed", linewidth = 0.25) +
  # lims(x=c(-0.5, 3.5)) +
  # lims(x=c(-0.5, 4.5)) +
  scale_fill_distiller(palette = "Oranges", direction = 1) +
  theme_dotplot +
  theme(
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(0,0,0,0),
    legend.position = "bottom")

cor_lims <- max(abs(gene_cor_df$pearson_cor)) %>% ifelse(. < 0.15, 0.15, .)
gene_cor_df$pt_shape <- ifelse(gene_cor_df$pearson_cor>0, "up", "down") %>% as.factor()
p2 <- ggplot(gene_cor_df, aes(y=reorder(variable, desc(variable)), x="x", fill=pearson_cor)) + 
  geom_tile(color=NA, width=1, height=1) +
  geom_point(data = subset(gene_cor_df, pt_shape=="up"), shape=24, fill="white", color="white", size=1.5) +
  geom_point(data = subset(gene_cor_df, pt_shape=="down"), shape=25, fill="white", color="white", size=1.5) +
  scale_fill_gradientn(colours = col_scale_div_custom2, limits=c(-cor_lims, cor_lims)) +
  theme_dotplot +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(), plot.margin = margin(0,0,0,0),
    legend.position = "right")


pdf(file = file.path(DIR_NN, paste0("hs_visium_A_NN_F14-C0_distance_", gene_oi,"_receptors_expr_correlation.pdf")), width = 3, height = 4)
(p1|p2) + patchwork::plot_layout(width = c(5,1))
dev.off()

#### Session Info ####
sessionInfo()
