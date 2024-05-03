#' [mm_visium_cell2location_proc_data.R]
#'
#' Run, analyse, and plot NMF results 
#'
#'
#' Jan 2023, L. Franzén [lovisa.franzen@scilifelab.se]

#### Set up ####
##### Define params. ####
set.seed(1)
SPECIES <- "mouse"
DIR_ROOT <- getwd()
DIR_DATA <- file.path(DIR_ROOT, "data", SPECIES, "sc_deconvolution_strunz")
DIR_RES <- file.path(DIR_ROOT, "results", SPECIES)
DIR_FIG_C2L <- file.path(DIR_RES, "figures", "cell2location_strunz")
dir.create(DIR_FIG_C2L)
fig_res <- 500

##### Load libs ####
library(STutility)
# library(harmony)
library(RColorBrewer)
library(patchwork)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(ggpmisc)
library(ggsankey)

##### Other ####
source(file.path(DIR_ROOT, "scripts", "custom_functions.R"))
source(file.path(DIR_ROOT, "scripts", "custom_colors.R"))
theme_custom <- theme(axis.title.x = element_blank())


##### Read se object ####
# fname <- paste0("mm_visium_preproc_se_obj.rds")
fname <- paste0("mm_visium_preproc_se_obj_nmf30_c2l.rds")
se <- readRDS(file = file.path(DIR_RES, "objects", fname))


##### Read tables ####
metadata <- read.table(file.path(DIR_DATA, "../visium/mm_visium_metadata.tsv"), sep = "\t", header = T)
rownames(metadata) <- metadata$sample_id
metadata$sample_n <- 1:nrow(metadata)


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
c2l_colnames <- colnames(c2l_all)

##### Add metadata to se obj ####
se <- AddMetaData(se, c2l_all)


##### Define variables ####
# factor_names <- paste0("factor_", 1:20)
factor_names <- paste0("factor_", 1:30)
c2l_names <- colnames(c2l_all)
new_c2l_names <- gsub("c2l_", "", gsub("[.]$", "", gsub("..", "_", c2l_names, fixed = TRUE)))


##### Violin & Spatial plots ####
png(file.path(DIR_FIG_C2L, paste0("mm_visium_c2l_res_violin_group.png")), width = 30*fig_res, height = 16*fig_res, res = fig_res)
VlnPlot(se, features = c2l_names, pt.size = 0, ncol = 10, group.by = "group", cols = cols_group) & theme_custom
dev.off()

png(file.path(DIR_FIG_C2L, paste0("mm_visium_c2l_res_violin_histpath.png")), width = 30*fig_res, height = 16*fig_res, res = fig_res)
VlnPlot(se, features = c2l_names, ncol = 10, group.by = "annotation", pt.size = 0) & theme_custom
dev.off()

png(file.path(DIR_FIG_C2L, paste0("mm_visium_c2l_res_violin_cluster.png")), width = 30*fig_res, height = 16*fig_res, res = fig_res)
VlnPlot(se, features = c2l_names, ncol = 10, group.by = "seurat_clusters", pt.size = 0) & theme_custom
dev.off()

png(file.path(DIR_FIG_C2L, paste0("mm_visium_c2l_res_spatial_Krt8ADI.png")), width = 14*fig_res, height = 14*fig_res, res = fig_res)
ST.FeaturePlot(se, features = grep("Krt8", c2l_names, value = T),
               ncol = 5, cols = col_scale_mako, pt.border = F, 
               min.cutoff = 1,
               label.by = "sample_name")
dev.off()


c2l_plot <- grep("c2l_", colnames(se@meta.data), value = T)
c2l_plot <- c("c2l_AT1.cells",
              "c2l_AT2.cells",
              "c2l_Ciliated.cells",
              "c2l_Fibroblasts",
              "c2l_Myofibroblasts",
              "c2l_B.lymphocytes",
              "c2l_T.lymphocytes",
              "c2l_Krt8.ADI")
pdf(file = file.path(DIR_FIG_C2L, paste0("mm_visium_c2l_res_spatial_selected_cells.pdf")), 
    width = 6, height = 6, useDingbats = F)
for(c in c2l_plot){ # grep("Krt8", c2l_plot, value = T)
  message(c)
  p <- ST.FeaturePlot(se, 
                      features = c, 
                      ncol = 5, 
                      cols = c("grey90", hcl.colors(6, palette = "BrwnYl", rev = T), "black"), 
                      pt.border = F, highlight.edges = F,
                      max.cutoff = 6, 
                      # min.cutoff = 1,
                      pt.size = 0.5,
                      label.by = "sample_name") &
    theme(aspect.ratio = 1)
  print(p)
}
dev.off()


pdf(file = file.path(DIR_FIG_C2L, paste0("mm_visium_c2l_res_spatial_selected_cells_samples.pdf")), 
    width = 7, height = 3, useDingbats = F)
ST.FeaturePlot(se, 
               indices = c(1,4,19),
               features = "c2l_AT2.cells",  # c2l_Myofibroblasts
               ncol = 3, 
               cols = c("grey95", col_scale_rocket), 
               pt.border = F, 
               min.cutoff = 0, 
               max.cutoff = 3,
               label.by = "sample_name") &
  theme(aspect.ratio = 1)
ST.FeaturePlot(se, 
               indices = c(1,4,19),
               features = "c2l_Fibroblasts",  # c2l_Myofibroblasts
               ncol = 3, 
               cols = c("grey95", col_scale_rocket),
               # cols = c("grey90", col_scale_mako), 
               pt.border = F, 
               min.cutoff = 0, 
               max.cutoff = 3,
               label.by = "sample_name") &
  theme(aspect.ratio = 1)
dev.off()

##### Boxplots selected cells ####
c2l_mdata <- se@meta.data
c2l_select <- c(
              "c2l_AT1.cells",
              "c2l_AT2.cells",
              "c2l_Ciliated.cells",
              "c2l_Fibroblasts",
              "c2l_Myofibroblasts",
              "c2l_Neutrophils",
              "c2l_B.lymphocytes",
              "c2l_T.lymphocytes",
              "c2l_Krt8.ADI")

###### Per section stats ######
c2l_cell_stats <- c2l_mdata %>% 
  group_by(sample_name, group, day, condition) %>% 
  summarise(.groups = "keep", 
            # n_sum_total = sum(colSums(across(c2l_colnames))),
            n_AT1 = mean(c2l_AT1.cells),
            n_AT2 = mean(c2l_AT2.cells),
            # n_Goblet = sum(c2l_Goblet.cells),
            # n_Club = mean(c2l_Club.cells),
            n_Ciliated = mean(c2l_Ciliated.cells),
            n_Fibroblast = mean(c2l_Fibroblasts),
            n_Myofib = mean(c2l_Myofibroblasts),
            n_Krt8ADI = mean(c2l_Krt8.ADI),
            # n_Neutrophil = mean(c2l_Neutrophils),
            n_RecMac = mean(c2l_Recruited.macrophages),
            n_ResMac = mean(c2l_Resolution.macrophages),
            n_Bcell = mean(c2l_B.lymphocytes),
            n_Tcell = mean(c2l_T.lymphocytes),
            )

c2l_cell_stats_long <- pivot_longer(c2l_cell_stats, 
                                    cols = grep("n_", colnames(c2l_cell_stats), value = T),
                                    names_to = "cell_type", 
                                    values_to = "n"
                                    )

c2l_cell_stats_long$day <- gsub("d", "", c2l_cell_stats_long$day)
c2l_cell_stats_long$day <- factor(c2l_cell_stats_long$day, levels = c("7", "21"))
c2l_cell_stats_long$condition <- factor(c2l_cell_stats_long$condition, levels = c("control", "bleomycin"))
c2l_cell_stats_long$cell_type <- gsub("n_", "", c2l_cell_stats_long$cell_type)
c2l_cell_stats_long$cell_type <- factor(c2l_cell_stats_long$cell_type,
                                        levels = c("AT1", "AT2",
                                                   "Fibroblast", "Myofib", "Krt8ADI",
                                                   # "Club", 
                                                   "Ciliated",
                                                   # "Neutrophil", 
                                                   "RecMac", "ResMac",
                                                   "Bcell", "Tcell"))
write.csv(c2l_cell_stats_long, file.path(DIR_FIG_C2L, paste0("mm_visium_c2l_res_nsum_boxplot_daata.csv")), row.names = F)
c2l_cell_stats_long <- read.csv(file.path(DIR_FIG_C2L, paste0("mm_visium_c2l_res_nsum_boxplot_daata.csv")))

# boxplot
pdf(file.path(DIR_FIG_C2L, paste0("mm_visium_c2l_res_nsum_boxplot2.pdf")), width = 4, height = 4)
ggplot(c2l_cell_stats_long, aes(x=day, y=n, fill=condition)) +
  # geom_point(size=2) +
  geom_boxplot(width=0.8, outlier.size = 0.25) +
  facet_wrap(~cell_type, nrow = 2) +
  scale_fill_manual(values = cols_cond) +
  scale_y_log10() +
  # ylim(0,3500) +
  labs(y="Avg. inferred cell density per section") +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5), 
        panel.spacing = unit(0, "lines"),
        panel.grid = element_blank(),
        # axis.title.x = element_blank(),
        axis.text = element_text(size=10, color="black"),
        legend.position = "top", 
        legend.title = element_blank())
dev.off()

# boxplot small
cells_select <- c2l_cell_stats_long$cell_type %>% levels()
cells_select <- cells_select[c(2,4,5,7)]

pdf(file.path(DIR_FIG_C2L, paste0("mm_visium_c2l_res_nsum_boxplot_small.pdf")), width = 8, height = 4)
ggplot(subset(c2l_cell_stats_long, cell_type %in% cells_select), 
       aes(x=day, y=n, fill=condition)) +
  geom_boxplot(width=0.8, outlier.size = 0.25) +
  facet_wrap(~cell_type, nrow = 1) +
  scale_fill_manual(values = cols_cond) +
  scale_y_log10() +
  labs(y="Avg. inferred cell\ndensity per section") +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5), 
        panel.spacing = unit(0, "lines"),
        panel.grid = element_blank(),
        # axis.title.x = element_blank(),
        axis.text = element_text(size=10, color="black"),
        legend.position = "right", 
        legend.title = element_blank(), aspect.ratio = 1)
dev.off()

# boxplot small +1
cells_select <- c2l_cell_stats_long$cell_type %>% levels()
cells_select <- cells_select[c(2,4,5,7,8)]

pdf(file.path(DIR_FIG_C2L, paste0("mm_visium_c2l_res_nsum_boxplot_small_v2.pdf")), width = 10, height = 4)
ggplot(subset(c2l_cell_stats_long, cell_type %in% cells_select), 
       aes(x=day, y=n, fill=condition)) +
  geom_boxplot(width=0.8, outlier.size = 0.25) +
  facet_wrap(~cell_type, nrow = 1) +
  scale_fill_manual(values = cols_cond) +
  scale_y_log10() +
  labs(y="Avg. inferred cell\ndensity per section") +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5), 
        panel.spacing = unit(0, "lines"),
        panel.grid = element_blank(),
        axis.text = element_text(size=10, color="black"),
        legend.position = "right", 
        legend.title = element_blank(), aspect.ratio = 1)
dev.off()


# Stat test for cells_select
subset(c2l_cell_stats_long, cell_type %in% cells_select)
cells_select <- c2l_cell_stats_long$cell_type %>% levels()
  
file_stat_test_res <- file.path(DIR_FIG_C2L, "mm_visium_c2l_res_nsum_boxplot_small_ttest_231129.txt")
file.create(file_stat_test_res)

cat(c("Welch Two Sample t-test performed between BLM and Veh for cell types included in the file 'mm_visium_c2l_res_nsum_boxplot_small.pdf'\n"), 
    file=file_stat_test_res, append = T, sep="\n")

for(c in cells_select){
  cat(paste0("# Testing for cell ", c), file=file_stat_test_res, sep="\n", append = T)
  message(c)
  for(d in c(7, 21)){
    cat(paste0("## Results for day ", d, " data:"), file=file_stat_test_res, sep="\n", append = T)
    test_data <- subset(c2l_cell_stats_long, cell_type %in% c & day == d)
    
    ttest_res <- t.test(subset(test_data, condition == "control")$n, subset(test_data, condition == "bleomycin")$n)
    res_formated <- paste0(names(unlist(ttest_res)), ": ", unlist(ttest_res)) %>% paste(collapse = "\n")
    cat(res_formated, 
        file=file_stat_test_res, sep="\n", append = T)
    cat("\n", 
        file=file_stat_test_res, sep="\n", append = T)
  }
}


# close(file_stat_test_res)


# dotplot
pdf(file.path(DIR_FIG_C2L, paste0("mm_visium_c2l_res_nsum_dotplot.pdf")), width = 8, height = 4)
ggplot(c2l_cell_stats_long) +
  geom_point(aes(x=sample_name, 
                 y=reorder(cell_type, as.character(cell_type)), 
                 fill=log10(n),
                 size=log10(n)),
             shape=21) +
  facet_wrap(~day, ncol = 2, scales = "free_x") +
  # scale_y_discrete(limits=rev) +
  scale_fill_gradientn(colours = col_scale_mako) +
  labs(fill="Sum of cell density") +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5), 
        # panel.spacing = unit(0, "lines"),
        panel.grid = element_blank(), 
        panel.background = element_rect(color=NA),
        strip.background = element_rect(color=NA),
        axis.title = element_blank(),
        axis.text = element_text(size=10, color="black"),
        axis.text.x = element_text(size=10, color="black", angle=45, hjust=1),
        legend.position = "right", 
        legend.title = element_blank())
dev.off()


###### Per animal stats ######
c2l_cell_stats_animal <- c2l_mdata %>% 
  group_by(animal, group, day, condition) %>% 
  summarise(.groups = "keep", 
            n_AT1 = mean(c2l_AT1.cells),
            n_AT2 = mean(c2l_AT2.cells),
            n_Ciliated = mean(c2l_Ciliated.cells),
            n_Fibroblast = mean(c2l_Fibroblasts),
            n_Myofib = mean(c2l_Myofibroblasts),
            n_Krt8ADI = mean(c2l_Krt8.ADI),
            n_RecMac = mean(c2l_Recruited.macrophages),
            n_ResMac = mean(c2l_Resolution.macrophages),
            n_Bcell = mean(c2l_B.lymphocytes),
            n_Tcell = mean(c2l_T.lymphocytes),
  )

c2l_cell_stats_animal_long <- pivot_longer(c2l_cell_stats_animal, 
                                    cols = grep("n_", colnames(c2l_cell_stats_animal), value = T),
                                    names_to = "cell_type", 
                                    values_to = "n"
)

c2l_cell_stats_animal_long$day <- gsub("d", "", c2l_cell_stats_animal_long$day)
c2l_cell_stats_animal_long$day <- factor(c2l_cell_stats_animal_long$day, levels = c("7", "21"))
c2l_cell_stats_animal_long$condition <- factor(c2l_cell_stats_animal_long$condition, levels = c("control", "bleomycin"))
c2l_cell_stats_animal_long$cell_type <- gsub("n_", "", c2l_cell_stats_animal_long$cell_type)
c2l_cell_stats_animal_long$cell_type <- factor(c2l_cell_stats_animal_long$cell_type,
                                        levels = c("AT1", "AT2",
                                                   "Fibroblast", "Myofib", "Krt8ADI",
                                                   "Ciliated",
                                                   "RecMac", "ResMac",
                                                   "Bcell", "Tcell"))

write.csv(c2l_cell_stats_animal_long, file.path(DIR_FIG_C2L, paste0("mm_visium_c2l_res_nsum_peranimal_boxplot_data.csv")), row.names = F)
# c2l_cell_stats_animal_long <- read.csv(file.path(DIR_FIG_C2L, paste0("mm_visium_c2l_res_nsum_peranimal_boxplot_data.csv")))


# boxplot small +1
cells_select <- c2l_cell_stats_animal_long$cell_type %>% levels()
cells_select <- cells_select[c(2,4,5,7,8)]

pdf(file.path(DIR_FIG_C2L, paste0("mm_visium_c2l_res_nsum_peranimal_boxplot_small_v2.pdf")), width = 8, height = 2.2)
ggplot(subset(c2l_cell_stats_animal_long, cell_type %in% cells_select), 
       aes(x=day, y=n, fill=condition)) +
  geom_boxplot(width=0.8, outlier.size = 0.2, color="black") +
  geom_point(data = subset(c2l_cell_stats_animal_long, cell_type %in% cells_select), 
             mapping = aes(x=day, y=n), 
             position = position_jitterdodge(jitter.width = 0), 
             size = 0.75) +
  facet_wrap(~cell_type, nrow = 1) +
  scale_fill_manual(values = cols_cond) +
  scale_y_log10() +
  labs(y="Avg. inferred cell\ndensity per animal") +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5), 
        strip.background = element_blank(),
        panel.spacing = unit(0, "lines"),
        panel.grid = element_blank(),
        axis.text = element_text(size=10, color="black"),
        legend.position = "right", 
        legend.title = element_blank())
dev.off()


# Stat test for cells_select
subset(c2l_cell_stats_animal_long, cell_type %in% cells_select)
cells_select <- c2l_cell_stats_animal_long$cell_type %>% levels()

file_stat_test_res <- file.path(DIR_FIG_C2L, "mm_visium_c2l_res_nsum_peranimal_boxplot_small_ttest_2023-12-08.txt")
file.create(file_stat_test_res)

cat(c("Welch Two Sample t-test performed between BLM and Veh for cell types included in the file 'mm_visium_c2l_res_nsum_boxplot_small.pdf'\n"), 
    file=file_stat_test_res, append = T, sep="\n")

for(c in cells_select){
  cat(paste0("# Testing for cell ", c), file=file_stat_test_res, sep="\n", append = T)
  message(c)
  for(d in c(7, 21)){
    cat(paste0("## Results for day ", d, " data:"), file=file_stat_test_res, sep="\n", append = T)
    test_data <- subset(c2l_cell_stats_animal_long, cell_type %in% c & day == d)
    
    ttest_res <- t.test(subset(test_data, condition == "control")$n, subset(test_data, condition == "bleomycin")$n)
    res_formated <- paste0(names(unlist(ttest_res)), ": ", unlist(ttest_res)) %>% paste(collapse = "\n")
    cat(res_formated, 
        file=file_stat_test_res, sep="\n", append = T)
    cat("\n", 
        file=file_stat_test_res, sep="\n", append = T)
  }
}



##### Zoom-in spatial HE overlay ####
# sections to select?
se_b21 <- SubsetSTData(se, sample_name == "d21_bleo_6b")
se_b21 <- LoadImages(se_b21, xdim = 1e3)

pdf(file = file.path(DIR_FIG_C2L, "mm_visium_c2l_res_spatial_selected_cells_HE_d21b6b.pdf"), 
    width = 6, height = 6, useDingbats = F)
FeatureOverlay(se_b21, 
               features = c("c2l_Fibroblasts"),
               add.alpha = T,
               cols = c("grey90", hcl.colors(6, palette = "OrYel", rev = T)), 
               max.cutoff = 6, 
               show.sb = T) &
  theme(aspect.ratio = 1)
FeatureOverlay(se_b21, 
               features = c("c2l_B.lymphocytes"),
               add.alpha = T,
               cols = c("grey90", hcl.colors(6, palette = "DarkMint", rev = T)), 
               max.cutoff = 6, 
               show.sb = T) &
  theme(aspect.ratio = 1)
dev.off()

#### All cell spot correlation ####
factor_cell_cor_data <- bind_cols(se@meta.data %>% select(sample_name, group, animal, day, condition, matches(c2l_names), matches(new_c2l_names)),
                                  se@reductions$NMF@cell.embeddings[, factor_names])

colnames(factor_cell_cor_data)[colnames(factor_cell_cor_data) %in% c2l_names] <- new_c2l_names

# write.csv(factor_cell_cor_data, file.path(DIR_RES, "objects", "NMF30", "mm_visium_nmf30_factorweight_cell2location_habermann_metadata.csv"))
factor_cell_cor_data <- read.csv(file = file.path(DIR_RES, "objects", "NMF30", "mm_visium_nmf30_factorweight_cell2location_habermann_metadata.csv"), 
                                 row.names = 1)

cell_anno <- read.csv(file = file.path(DIR_ROOT, "data", "misc", "strunz_cell_type_groups.csv"), header = T, sep = ";")


##### Cell~cell correlation #####
cell_cor_data <- factor_cell_cor_data %>% select(sample_name, day, animal, condition, group, matches(c2l_names), matches(new_c2l_names))
colnames(cell_cor_data)[colnames(cell_cor_data) %in% c2l_names] <- new_c2l_names

# Filter cells to show
new_c2l_names_filtered <- new_c2l_names[!new_c2l_names %in% c("NA", "low.quality.cells")]

cell_anno_filt <- cell_anno[cell_anno$cell_name%in%new_c2l_names_filtered, ]


###### All samples ######
cell_density_hm <- cor(cell_cor_data[,new_c2l_names_filtered])
diag(cell_density_hm) <- 0
cell_density_hm <- as.data.frame(cell_density_hm)


#' plot pheatmap
# max(cell_density_hm)
cell_density_hm_plot <- cell_density_hm[,new_c2l_names_filtered]

rownames(cell_density_hm_plot) <- cell_anno_filt$cell_name_special
colnames(cell_density_hm_plot) <- cell_anno_filt$cell_name_special

cell_group_colors <- data.frame(group = cell_anno_filt$cell_group, 
                                compartment = cell_anno_filt$lung_compartment_healthy,
                                row.names = cell_anno_filt$cell_name_special)

cell_colors <- list(group = c(brewer.pal(length(unique(cell_group_colors$group))-1, "Pastel1"), "grey90"),
                    compartment = c("#357BA2FF", "#49C1ADFF", "#DEF5E5FF", "grey90")
                      #c(brewer.pal(length(unique(cell_group_colors$compartment))-1, "Dark2"), "grey90")
                      )
names(cell_colors$group) <- unique(cell_group_colors$group)
names(cell_colors$compartment) <- unique(cell_group_colors$compartment)

hclusters <- hclust(dist(cell_density_hm_plot, ))
plot(hclusters)
cluster_cut_h1_nclust <- cutree(hclusters, h = 1.25) %>% max()
cluster_cut_h1_nclust <- cutree(hclusters, h = 1.5) %>% max()

# pal_length <- 11
# ph_colors <- RColorBrewer::brewer.pal(pal_length, "RdBu") %>% rev()
pal_length <- length(col_scale_div_custom2)
ph_colors <- col_scale_div_custom2
abs_max_val <- max(abs(cell_density_hm_plot))
ph_breaks <- c(seq(-abs_max_val, 0, length.out=ceiling(pal_length/2) + 1), 
               seq(abs_max_val/pal_length, abs_max_val, length.out=floor(pal_length/2)))

pdf(file = file.path(DIR_FIG_C2L, "mm_visium_c2l_res_cell_type_cor_heatmap.pdf"), width = 10, height = 10, useDingbats = F)
pheatmap::pheatmap(cell_density_hm_plot,
                   cellwidth = 10, 
                   cellheight = 10, 
                   annotation_names_row = F, show_rownames = F,
                   annotation_col = cell_group_colors,
                   # annotation_row = cell_group_colors,
                   annotation_colors = cell_colors, 
                   treeheight_col = 10,
                   treeheight_row = 0,
                   cutree_cols = cluster_cut_h1_nclust,
                   cutree_rows = cluster_cut_h1_nclust,
                   color = ph_colors, 
                   breaks = ph_breaks)
dev.off()

# ggplot
cell_density_hm$cell_1 <- rownames(cell_density_hm)
cell_density_hm_long <- pivot_longer(cell_density_hm, cols = all_of(new_c2l_names_filtered), names_to = "cell_2", values_to = "cor")

p_hm <- ggplot(cell_density_hm_long, aes(cell_1, cell_2, fill= cor)) + 
  geom_tile(color = "white",
            lwd = .5,
            linetype = 1) +
  scale_fill_gradient2(low = "#0474BA",
                       mid = "white",
                       high = "#F17720") +
  labs(x="", y="", title="", fill="correlation") +
  theme_bw() +
  coord_fixed() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust=0.5, size=12, face = "bold"),
        axis.text = element_text(size=10),
        legend.title = element_text(size=10), 
        legend.text = element_text(size=10))

pdf(file = file.path(DIR_FIG_C2L, "mm_visium_c2l_res_cell_type_cor_heatmap2.pdf"), width = 10, height = 10, useDingbats = F)
pheatmap::pheatmap(cell_density_hm[,new_c2l_names_filtered], cellwidth = 10, cellheight = 10, color = ph_colors, breaks = ph_breaks)
print(p_hm)
dev.off()

###### Fig 4e: Per group, Compartment analysis ###### 
####### > Heatmap ######
groups_use <- c("CTRL", "d7_BLM", "d21_BLM")

cell_cor_data_group <- setNames(lapply(groups_use, function(g){
  if(g=="CTRL"){
    cell_cor_data_subset <- subset(cell_cor_data, condition == "control")
  } else {
    cell_cor_data_subset <- subset(cell_cor_data, group == g)
  }
  cell_density_hm <- cor(cell_cor_data_subset[,new_c2l_names_filtered])
  diag(cell_density_hm) <- 0
  cell_density_hm <- as.data.frame(cell_density_hm)
  }), 
  nm = groups_use)

# Plot pheatmap
ph_colors <- hcl.colors(n = pal_length, palette = "PuOr") %>% rev()
ph_colors <- hcl.colors(n = pal_length, palette = "PiYG") %>% rev()
ph_colors <- hcl.colors(n = pal_length, palette = "RdBu") %>% rev()

ph_colors <- col_scale_div_custom2
pal_length <- length(ph_colors)
abs_max_val <- max(abs(cell_density_hm[,new_c2l_names_filtered]))
ph_breaks <- c(seq(-abs_max_val, 0, length.out=ceiling(pal_length/2) + 1), 
               seq(abs_max_val/pal_length, abs_max_val, length.out=floor(pal_length/2)))

pdf(file = file.path(DIR_FIG_C2L, "mm_visium_c2l_res_cell_type_cor_heatmap_per_group_new3.pdf"), width = 10, height = 10, useDingbats = F)
for(g in groups_use){
  cell_density_hm_plot <- cell_cor_data_group[[g]][,new_c2l_names_filtered]
  cell_anno_filt <- cell_anno[cell_anno$cell_name%in%rownames(cell_density_hm_plot), ]
  rownames(cell_density_hm_plot) <- cell_anno_filt$cell_name_special
  colnames(cell_density_hm_plot) <- cell_anno_filt$cell_name_special

#   abs_max_val <- max(abs(cell_density_hm_plot))
#   ph_breaks <- c(seq(-abs_max_val, 0, length.out=ceiling(pal_length/2) + 1),
#                  seq(abs_max_val/pal_length, abs_max_val, length.out=floor(pal_length/2)))
  
  h_cutoff <- 1.5
  n_clust <- hclust(dist(x = cell_density_hm_plot)) %>% cutree(h = h_cutoff) %>% max()
  
  print(
    pheatmap::pheatmap(cell_density_hm_plot,
                       annotation_names_row = F, show_rownames = F,
                       annotation_col = cell_group_colors, # %>% select(compartment),
                       # annotation_row = cell_group_colors,
                       annotation_colors = cell_colors, 
                       treeheight_col = 16,
                       treeheight_row = 0, 
                       cutree_rows = n_clust,
                       cutree_cols = n_clust,
                       cellwidth = 10, 
                       cellheight = 10, border_color = NA,
                       color = ph_colors, 
                       breaks = ph_breaks, 
                       main = g)
  )
}
dev.off()


####### > Sankey #######
#' Define tree cutoff
h_cutoff <- 1.5
hclusters_group <- setNames(lapply(groups_use, function(g){
  cell_density_hm_plot <- cell_cor_data_group[[g]][,new_c2l_names_filtered]
  rownames(cell_density_hm_plot) <- cell_anno_filt$cell_name_special
  colnames(cell_density_hm_plot) <- cell_anno_filt$cell_name_special
  plot( hclust(dist(x = cell_density_hm_plot)))
  hclusters <- hclust(dist(x = cell_density_hm_plot)) %>% cutree(h = h_cutoff)
  hclusters <- paste0(hclusters, ".", g)
  return(hclusters)
  }), nm = groups_use) %>% 
  bind_cols()
hclusters_group$cell_type <- cell_anno_filt$cell_name_special

# max(hclusters_group$CTRL);max(hclusters_group$d7_BLM);max(hclusters_group$d21_BLM)

hclusters_group <- hclusters_group %>% 
  mutate_if(is.character, as.factor) %>% 
  mutate_if(is.integer, as.factor)

hclusters_group_long <- pivot_longer(hclusters_group, cols = c("CTRL", "d7_BLM", "d21_BLM"), names_to = "group", values_to = "cluster")
hclusters_group_long <- merge(x=hclusters_group_long, y=cell_anno_filt, by.x= "cell_type", by.y = "cell_name_special")

write.csv(hclusters_group, file = file.path(DIR_FIG_C2L, "mm_visium_c2l_res_cell_type_cor_heatmap_per_group_hclust_data.csv"))
write.csv(hclusters_group_long, file = file.path(DIR_FIG_C2L, "mm_visium_c2l_res_cell_type_cor_heatmap_per_group_hclust_data_long.csv"))

hclusters_group <- read.csv(file.path(DIR_FIG_C2L, "mm_visium_c2l_res_cell_type_cor_heatmap_per_group_hclust_data.csv"))

sank_plot_df <- hclusters_group %>% ggsankey::make_long(CTRL, d7_BLM, d21_BLM)

node_colors <- setNames(
  object = 
    # c("#357BA2", "#dea4c1",  "#38629DFF", "#49C1AD",  "#27436b",  "#926390",  "#B1E4C2FF",   "#78D6AEFF"),
    c("#5DA6C6", "#E391BA",  "#3873B9", "#E7B922",  "#215189",  "#8D71A6",  "#F4D8A8",   "#EDC962"),
  nm = 
    c("1.CTRL", "1.d21_BLM", "1.d7_BLM", "2.CTRL", "2.d21_BLM", "2.d7_BLM", "3.d21_BLM", "3.d7_BLM"))

p_sank <- ggplot(sank_plot_df, aes(x = x, 
               next_x = next_x, 
               node = node, 
               next_node = next_node,
               fill = factor(node))) +
  geom_sankey(flow.alpha = 0.6, node.color = 1, width = 0.25) +
  scale_fill_manual(values = t(node_colors)) +
  # scale_fill_viridis_d() +
  # scale_fill_viridis_d(option = "F", direction = -1) +
  theme_sankey(base_size = 12, base_line_size = 0.1) +
  theme(legend.position = "bottom", axis.title = element_blank());p_sank

pdf(file = file.path(DIR_FIG_C2L, "mm_visium_c2l_res_cell_type_cor_heatmap_per_group_hclust_sankey.pdf"), width = 6, height = 4, useDingbats = F)
p_sank
dev.off()


####### > Spatial compartment scores ######
hclusters_group_long <- read.csv(file = file.path(DIR_FIG_C2L, "mm_visium_c2l_res_cell_type_cor_heatmap_per_group_hclust_data_long.csv"))

cols_cluster <- 
  # setNames(lapply(unique(hclusters_group_long$group), function(g){
  # df <- subset(hclusters_group_long, group == g)
  # cluster_cells <- 
  setNames(lapply(unique(hclusters_group_long$cluster), function(c){
    subset(hclusters_group_long, cluster == c)$c2l_col_name
  }), nm=unique(hclusters_group_long$cluster))
  # }), nm=unique(hclusters_group_long$group))

mdat <- se@meta.data

mdat %>% select(matches(cols_cluster$`1.CTRL`)) %>% head()

for(n in unique(names(cols_cluster))){
  se <- AddMetaData(object = se, 
                    metadata = (mdat %>% mutate(sum = rowSums(across(cols_cluster[[n]]))) %>% select(sum)),
                    col.name = paste0("comp_", n))
}


col_plot <- grep("comp_", colnames(se@meta.data), value = T)
write.csv(se@meta.data[, col_plot], file = file.path(DIR_FIG_C2L, paste0("mm_visium_c2l_res_cell_type_cor_heatmap_per_group_hclust_score_metadata.csv")), row.names = T)


pdf(file = file.path(DIR_FIG_C2L, paste0("mm_visium_c2l_res_cell_type_cor_heatmap_per_group_hclust_score_spatial.pdf")), 
    width = 7, height = 7, useDingbats = F)
for(c in col_plot[3]){
  message(c)
  p <- ST.FeaturePlot(se, 
                      features = c, 
                      ncol = 5, 
                      cols = c("grey90", hcl.colors(6, palette = "BrwnYl", rev = T)), 
                      max.cutoff = 15, 
                      min.cutoff = 1,
                      pt.size = 0.5,
                      label.by = "sample_name") &
    theme(aspect.ratio = 1)
  print(p)
}
dev.off()


pdf(file = file.path(DIR_FIG_C2L, paste0("mm_visium_c2l_res_cell_type_cor_heatmap_per_group_hclust_score_spatial_d21_b4a.pdf")), 
    width = 7, height = 4, useDingbats = F)
ST.FeaturePlot(se, 
               features = c("comp_1.d21_BLM", "comp_2.d21_BLM", "comp_3.d21_BLM"), 
               ncol = 3,
               indices = 19,
               # cols = col_scale_spec %>% rev(),
               cols = hcl.colors(6, palette = "RdYlBu", rev = T),
               max.cutoff = 10, 
               min.cutoff = 1,
               pt.size = 1.1,
               label.by = "sample_name") &
  theme(aspect.ratio = 1, legend.position = "bottom")
dev.off()


ST.FeaturePlot(se, 
               features = c("comp_1.d7_BLM", "comp_2.d7_BLM", "comp_3.d7_BLM"), 
               ncol = 3,
               indices = 19,
               # cols = col_scale_spec %>% rev(),
               cols = hcl.colors(6, palette = "RdYlBu", rev = T),
               max.cutoff = 10, 
               min.cutoff = 2,
               pt.size = 1.1,
               label.by = "sample_name") &
  theme(aspect.ratio = 1, legend.position = "bottom")


###### Per animal (heatmap) ######
cell_cor_data$animal <- factor(cell_cor_data$animal, levels = c(paste0("d7_c", 1:3), paste0("d7_b", 1:6), paste0("d21_c", 1:3), paste0("d21_b", 1:6)))
animals_use <- levels(cell_cor_data$animal)
cell_cor_data_animal <- setNames(lapply(animals_use, function(a){
  cell_cor_data_subset <- subset(cell_cor_data, animal == a)
  cell_density_hm <- cor(cell_cor_data_subset[,new_c2l_names_filtered])
  diag(cell_density_hm) <- 0
  cell_density_hm <- as.data.frame(cell_density_hm)
}), 
nm = animals_use)

cell_cor_data_animal_all <- do.call(cbind, cell_cor_data_animal)


##### Krt8 ADI ~ c2l corr #####
cell_cor_data_group_abba <- do.call(cbind, cell_cor_data_animal) %>% select(contains("Krt8"))
cell_cor_data_group_abba %>% head()
cell_cor_data_group_abba <- cell_cor_data_group_abba[rownames(cell_cor_data_group_abba)!="Krt8.ADI", ]

mdat <- metadata %>% group_by(animal, group) %>% summarise(n=n())
plot_groups <- data.frame(row.names = paste0(mdat$animal, ".Krt8.ADI"),
                          group = mdat$group)
group_colors = list(
  group = cols_group
)
pal_length <- length(col_scale_div_custom2)
ph_colors <- col_scale_div_custom2
abs_max_val <- max(cell_cor_data_group_abba)
ph_breaks <- c(seq(-abs_max_val, 0, length.out=ceiling(pal_length/2) + 1), 
               seq(abs_max_val/pal_length, abs_max_val, length.out=floor(pal_length/2)))

pdf(file = file.path(DIR_FIG_C2L, paste0("mm_visium_c2l_res_cell_type_cor_heatmap_per_animal_Krt8ADI.pdf")), 
    width = 8, height = 8, useDingbats = F)
pheatmap::pheatmap(cell_cor_data_group_abba, 
                   cellwidth = 12, 
                   cellheight = 10, 
                   treeheight_row = 10,  
                   treeheight_col = 4, 
                   annotation_col = plot_groups, 
                   annotation_colors = group_colors, 
                   # cluster_cols = F,
                   color = ph_colors,
                   breaks = ph_breaks
)
dev.off()

pdf(file = file.path(DIR_FIG_C2L, paste0("mm_visium_c2l_res_cell_type_cor_heatmap_per_animal_Krt8ADI_onlyBLM.pdf")), 
    width = 8, height = 8, useDingbats = F)
pheatmap::pheatmap(cell_cor_data_group_abba[, grep("_b", colnames(cell_cor_data_group_abba))], 
                   cellwidth = 12, 
                   cellheight = 10, 
                   treeheight_row = 10,  
                   treeheight_col = 4, 
                   annotation_col = plot_groups,
                   annotation_colors = group_colors, 
                   # cluster_cols = F,
                   color = ph_colors,
                   breaks = ph_breaks
)
dev.off()


##### Fig 4g: cell~factor d7/d21 correlation #####
# Read in se objects
se_d7 <- readRDS(file = file.path(DIR_RES, "objects", "mm_visium_preproc_se_obj_subset_d7_nmf30_c2l.rds"))
se_d21 <- readRDS(file = file.path(DIR_RES, "objects", "mm_visium_preproc_se_obj_subset_d21_nmf30_c2l.rds"))

# prepare data
factor_cell_cor_data_d7 <- bind_cols(
  se_d7@meta.data %>% select(sample_name, group, animal, day, condition, matches(c2l_names), matches(new_c2l_names)),
  se_d7@reductions$NMF@cell.embeddings[, factor_names])

factor_cell_cor_data_d21 <- bind_cols(
  se_d21@meta.data %>% select(sample_name, group, animal, day, condition, matches(c2l_names), matches(new_c2l_names)),
  se_d21@reductions$NMF@cell.embeddings[, factor_names])


# compute correlations
groups_use <- c("d7", "d21")
factor_cell_cor_data_group <- setNames(lapply(groups_use, function(g){
  if (g=="d7") {
    factor_cell_density_hm <- cor(factor_cell_cor_data_d7[,c(c2l_names, factor_names)])
  } else if (g=="d21"){
    factor_cell_density_hm <- cor(factor_cell_cor_data_d21[,c(c2l_names, factor_names)])
  }
  # format output
  diag(factor_cell_density_hm) <- 0
  factor_cell_density_hm <- as.data.frame(factor_cell_density_hm)
  factor_cell_density_hm <- factor_cell_density_hm[factor_names, c2l_names]
  
  # rename factor names
  rownames(factor_cell_density_hm) <- gsub("factor_", paste0(g, "-F"), rownames(factor_cell_density_hm))
  
  # filter and rename cell type names
  colnames(factor_cell_density_hm) <- new_c2l_names
  factor_cell_density_hm <- factor_cell_density_hm[, new_c2l_names_filtered]
  cell_anno_filt <- cell_anno[cell_anno$cell_name %in% colnames(factor_cell_density_hm), ]
  colnames(factor_cell_density_hm) <- cell_anno_filt$cell_name_special
  
  return(t(factor_cell_density_hm))
  }),
  nm = groups_use)


# Cell group annotations from previous analysis (per group hclust / sankey)
# group_annotation_color <- setNames(
#   object = 
#     c("#5DA6C6", "#E391BA",  "#3873B9", "#E7B922",  "#215189",  "#8D71A6",  "#F4D8A8",   "#EDC962"),
#   nm = 
#     c("1.CTRL", "1.d21_BLM", "1.d7_BLM", "2.CTRL", "2.d21_BLM", "2.d7_BLM", "3.d21_BLM", "3.d7_BLM"))

group_annotation_color <- list(
  CTRL = setNames(c("#5DA6C6", "#E7B922"), nm=c("1.CTRL", "2.CTRL")),
  d7_BLM = setNames(c("#3873B9", "#8D71A6", "#EDC962"), nm = c("1.d7_BLM", "2.d7_BLM", "3.d7_BLM")),
  d21_BLM = setNames(c("#E391BA", "#215189", "#F4D8A8"), nm = c("1.d21_BLM", "2.d21_BLM", "3.d21_BLM"))
)

group_annotation_cells <- hclusters_group %>% as.data.frame()
rownames(group_annotation_cells) <- group_annotation_cells$cell_type
group_annotation_cells$cell_type <- NULL
group_annotation_cells <- group_annotation_cells[,rev(colnames(group_annotation_cells))]
# test heatmap
# plot_groups <- data.frame(row.names = paste0(mdat$animal, ".Krt8.ADI"),
#                           group = mdat$group)
# group_colors = list(
#   group = cols_group
# )

ph_colors <- col_scale_div_custom2
pal_length <- length(ph_colors)

abs_max_val <- max(cbind(factor_cell_cor_data_group[[1]], factor_cell_cor_data_group[[2]]))
ph_breaks <- c(seq(-abs_max_val, 0, length.out=ceiling(pal_length/2) + 1), 
               seq(abs_max_val/pal_length, abs_max_val, length.out=floor(pal_length/2)))

pdf(file = file.path(DIR_FIG_C2L, paste0("mm_visium_c2l_res_cell_type_dayNMF_cor_heatmap_grouped_raw.pdf")), 
    width = 10, height = 8, useDingbats = F)
pheatmap::pheatmap(factor_cell_cor_data_group[[1]], 
                   cellwidth = 10, 
                   cellheight = 10, 
                   treeheight_row = 10,  
                   treeheight_col = 4, 
                   annotation_row = group_annotation_cells,
                   annotation_colors = group_annotation_color,
                   # cluster_cols = F,
                   color = ph_colors,
                   breaks = ph_breaks
                   )
pheatmap::pheatmap(factor_cell_cor_data_group[[2]], 
                   cellwidth = 10, 
                   cellheight = 10, 
                   treeheight_row = 10,  
                   treeheight_col = 4, 
                   annotation_row = group_annotation_cells,
                   annotation_colors = group_annotation_color,
                   # cluster_cols = F,
                   color = ph_colors,
                   breaks = ph_breaks
)
dev.off()

pdf(file = file.path(DIR_FIG_C2L, paste0("mm_visium_c2l_res_cell_type_dayNMF_cor_heatmap_grouped_formated.pdf")), 
    width = 10, height = 4, useDingbats = F)
pheatmap::pheatmap(factor_cell_cor_data_group[[1]], 
                   cellwidth = 10, 
                   cellheight = 3, 
                   border_color = NA,
                   treeheight_row = 6,  
                   treeheight_col = 10, 
                   show_rownames = F,
                   annotation_row = group_annotation_cells,
                   annotation_colors = group_annotation_color,
                   color = ph_colors,
                   breaks = ph_breaks
                   )
pheatmap::pheatmap(factor_cell_cor_data_group[[2]], 
                   cellwidth = 10,
                   cellheight = 3, 
                   border_color = NA,
                   treeheight_row = 6,  
                   treeheight_col = 10, 
                   show_rownames = F,
                   annotation_row = group_annotation_cells,
                   annotation_colors = group_annotation_color,
                   color = ph_colors,
                   breaks = ph_breaks
)
dev.off()


##### cell~factor correlation #####
# factor_cell_cor_data <- bind_cols(se@meta.data %>% select(sample_name, day, condition, group, matches(c2l_names)),
#                                   se@reductions$NMF@cell.embeddings[, factor_names])
# 
# colnames(factor_cell_cor_data)[colnames(factor_cell_cor_data) %in% c2l_names] <- new_c2l_names
# 
# 
# # Filter cells to show
# new_c2l_names_filtered <- new_c2l_names[!new_c2l_names %in% c("NA", "low.quality.cells")]
# 
# 
# groups_use <- c("all", "CTRL", "d7_BLM", "d21_BLM")
# 
# factor_cell_cor_data_group <- setNames(lapply(groups_use, function(g){
#   if (g=="all") {
#     factor_cell_cor_subset <- factor_cell_cor_data
#   } else if (g=="CTRL") {
#     factor_cell_cor_subset <- subset(factor_cell_cor_data, condition == "control")
#   } else {
#     factor_cell_cor_subset <- subset(factor_cell_cor_data, group == g)
#   }
#   factor_cell_density_hm <- cor(factor_cell_cor_subset[,c(new_c2l_names_filtered, factor_names)])
#   diag(factor_cell_density_hm) <- 0
#   factor_cell_density_hm <- as.data.frame(factor_cell_density_hm)
#   factor_cell_density_hm <- factor_cell_density_hm[factor_names, new_c2l_names_filtered]
#   }), 
#   nm = groups_use)
# 
# 
# # Plot pheatmap
# pal_length <- 11
# ph_colors <- RColorBrewer::brewer.pal(pal_length, "RdBu") %>% rev()
# 
# pdf(file = file.path(DIR_FIG_C2L, "mm_visium_c2l_res_cell_type_NMF20_factor_cor_heatmap_per_group.pdf"), width = 10, height = 10, useDingbats = F)
# for(g in groups_use){
#   factor_cell_density_hm <- factor_cell_cor_data_group[[g]]
#   
#   abs_max_val <- max(abs(factor_cell_density_hm))
#   ph_breaks <- c(seq(-abs_max_val, 0, length.out=ceiling(pal_length/2) + 1), 
#                  seq(abs_max_val/pal_length, abs_max_val, length.out=floor(pal_length/2)))
#   
#   print(
#     pheatmap::pheatmap(t(factor_cell_density_hm), 
#                        cellwidth = 10, 
#                        cellheight = 10, 
#                        color = ph_colors, 
#                        breaks = ph_breaks, 
#                        main = g)
#   )
# }
# dev.off()


#### Extra: Krt8 ADI ####

##### Count cells per sample #####
c2l_mdata <- se@meta.data

c2l_krt8adi_stats <- c2l_mdata %>% 
  group_by(sample_name, group, day, condition) %>% 
  summarise(.groups = "keep", 
            n_cells_krt8adi = sum(c2l_Krt8.ADI),
            mean_density_krt8ad = mean(c2l_Krt8.ADI))
c2l_krt8adi_stats$day <- factor(c2l_krt8adi_stats$day, levels = c("d7", "d21"))

ggplot(c2l_krt8adi_stats, aes(x=day, y=n_cells_krt8adi, color=condition)) +
  # geom_point(size=2) +
  geom_boxplot(width=0.4) +
  scale_color_manual(values = cols_cond) +
  labs(y="Sum of cell density", title="Krt8+ADI density\nper sample") +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5))

t.test(x = subset(c2l_krt8adi_stats, day=="d7" & condition == "control")$n_cells_krt8adi, 
       y = subset(c2l_krt8adi_stats, day=="d7" & condition == "bleomycin")$n_cells_krt8adi, 
       alternative = "less")

t.test(x = subset(c2l_krt8adi_stats, day=="d21" & condition == "control")$n_cells_krt8adi, 
       y = subset(c2l_krt8adi_stats, day=="d21" & condition == "bleomycin")$n_cells_krt8adi, 
       alternative = "less")

t.test(x = subset(c2l_krt8adi_stats, day=="d7" & condition == "bleomycin")$n_cells_krt8adi, 
       y = subset(c2l_krt8adi_stats, day=="d21" & condition == "bleomycin")$n_cells_krt8adi, 
       alternative = "greater")



##### H&E Spatial #####
# se_he <- SubsetSTData(se, sample_name %in% c("d7_bleo_3", "d7_bleo_5b", "d21_bleo_3", "d21_bleo_4a")) %>%
#   LoadImages(se_he, xdim = 1e3)

se_he_d7 <- SubsetSTData(se, sample_name %in% c("d7_bleo_1", "d7_bleo_2", "d7_bleo_3", "d7_bleo_6b")) %>% 
  LoadImages(xdim = 1e3)

# se_he_d21 <- SubsetSTData(se, sample_name %in% c("d21_bleo_1", "d21_bleo_2", "d21_bleo_5b", "d21_bleo_6b")) %>% 
#   LoadImages(xdim = 1e3)
se_he_d21 <- SubsetSTData(se, sample_name %in% c("d21_bleo_1", "d21_bleo_3", "d21_bleo_4a", "d21_bleo_5b")) %>%
  LoadImages(xdim = 1e3)

# se_he <- se_he_d7
se_he <- se_he_d21

# ImagePlot(se_he)
p <- FeatureOverlay(se_he, 
                    sampleids = 1:4, features = "c2l_Krt8.ADI", 
                    label.by = "sample_name", value.scale = "all", 
                    min.cutoff = 0.5, pt.size = 1.8,
                    ncols = 2, cols = (col_scale_mako), add.alpha = T) &
  theme(aspect.ratio = 1)
p <- ggrastr::rasterize(p, layers = "Point", dpi = fig_res)

png(file.path(DIR_FIG_C2L, paste0("mm_visium_c2l_res_spatial_Krt8ADI_d21_he.png")), width = 14*fig_res, height = 14*fig_res, res = fig_res);p;dev.off()


p1 <- FeatureOverlay(se_he, sampleids = 1:4, 
                     features = "factor_17", 
                     label.by = "sample_name", value.scale = "all", 
                     min.cutoff = 0.5, pt.size = 1.8,
                     ncols = 2, cols = (col_scale_mako), add.alpha = T) &
  theme(aspect.ratio = 1)
p1 <- ggrastr::rasterize(p1, layers = "Point", dpi = fig_res)

p1.2 <- FeatureOverlay(se_he, sampleids = 1:4, 
                     features = "factor_13", 
                     label.by = "sample_name", value.scale = "all", 
                     min.cutoff = 0.5, pt.size = 1.8,
                     ncols = 2, cols = (col_scale_mako), add.alpha = T) &
  theme(aspect.ratio = 1)
p1.2 <- ggrastr::rasterize(p1.2, layers = "Point", dpi = fig_res)

p1.3 <- FeatureOverlay(se_he, sampleids = 1:4, 
                       features = "factor_9", 
                       label.by = "sample_name", value.scale = "all", 
                       min.cutoff = 0.25, pt.size = 1.8,
                       ncols = 2, cols = (col_scale_mako), add.alpha = T) &
  theme(aspect.ratio = 1)
p1.3 <- ggrastr::rasterize(p1.3, layers = "Point", dpi = fig_res)


p_apoe <- FeatureOverlay(se_he, sampleids = 1:4, 
                       features = "Apoe", 
                       label.by = "sample_name", value.scale = "all", 
                       min.cutoff = 3, 
                       pt.size = 1.8,
                       ncols = 2, cols = (col_scale_mako), add.alpha = T) &
  theme(aspect.ratio = 1)
p_apoe <- ggrastr::rasterize(p_apoe, layers = "Point", dpi = fig_res)


p2 <- FeatureOverlay(se_he, sampleids = 1:4, 
                     features = "annotation", 
                     label.by = "sample_name", 
                     # value.scale = "all", 
                     # min.cutoff = 0.5, 
                     pt.size = 1.8, 
                     cols = cols_annotation,
                     ncols = 2) &
  theme(aspect.ratio = 1)
p2 <- ggrastr::rasterize(p2, layers = "Point", dpi = fig_res)


pdf(file.path(DIR_FIG_C2L, paste0("mm_visium_c2l_res_spatial_Krt8ADI_he_selection_d7.pdf")), width = 14, height = 14)
# pdf(file.path(DIR_FIG_C2L, paste0("mm_visium_c2l_res_spatial_Krt8ADI_he_selection_d21.pdf")), width = 14, height = 14)
p
p1
p1.2
p1.3
p_apoe
p2
dev.off()





##### Krt8 ADI ~ NMF corr #####

#' Scatter plot with correlations per sample
pdf(file = file.path(DIR_FIG_C2L, "mm_visium_c2l_res_krt8adi_NMF_cor_per_sample.pdf"), width = 14, height = 10, useDingbats = F)
for(s_name in metadata$sample_name){
  message(s_name)
  spots_select <- rownames(subset(se@meta.data, sample_name == s_name))
  dat_nmf <- se@reductions$NMF@cell.embeddings[spots_select, ]
  dat_nmf <- as.data.frame(dat_nmf)
  dat_nmf$Krt8.ADI <- se@meta.data[spots_select, grep("Krt8", c2l_names, value = T)]
  dat_nmf <- dat_nmf %>% 
    mutate_if(is.character, as.numeric) %>% 
    mutate_if(is.factor, as.numeric)
  
  p_list <- lapply(factor_names, function(f){
    p <- ggplot(dat_nmf[dat_nmf[[f]]>0,], aes_string(x="Krt8.ADI", y=f)) +
      geom_point(size=0.5, alpha = 0.5) +
      labs(title=f) +
      scale_y_log10() +
      scale_x_log10() +
      stat_poly_line() +
      stat_poly_eq(color="grey70") +
      theme_bw()
    p <- ggrastr::rasterize(p, layers = "Point", dpi = fig_res)
  })
  print(wrap_plots(p_list, ncol = 5) & patchwork::plot_annotation(title = s_name))
}
dev.off()

#' Pearson's Correlation r values per sample
cor_list <- setNames(
  lapply(metadata$sample_name, function(s_name){
    message(s_name)
    spots_select <- rownames(subset(se@meta.data, sample_name == s_name))
    dat_nmf <- se@reductions$NMF@cell.embeddings[spots_select, ]
    dat_nmf <- as.data.frame(dat_nmf)
    dat_nmf$Krt8.ADI <- se@meta.data[spots_select, grep("Krt8", c2l_names, value = T)]
    dat_nmf <- dat_nmf %>% 
      mutate_if(is.character, as.numeric) %>% 
      mutate_if(is.factor, as.numeric)
    lapply(factor_names, function(f){cor(x = dat_nmf[["Krt8.ADI"]], y = dat_nmf[[f]])}) %>% unlist()
  }), nm = metadata$sample_name)

cor_df <- bind_cols(cor_list)
cor_df <- cor_df %>% t()
colnames(cor_df) <- factor_names
cor_df <- as.data.frame(cor_df) 
cor_df$sample_name <- rownames(cor_df)
cor_df <- merge(x = cor_df, y = metadata, by = "sample_name")

cor_df_long <- tidyr::pivot_longer(cor_df, cols = factor_names, names_to = "factor", values_to = "value")
cor_df_long$group <- factor(cor_df_long$group, levels = c("d7_CTRL","d7_BLM","d21_CTRL","d21_BLM"))
cor_df_long$factor <- factor(cor_df_long$factor, levels = factor_names)

#' Boxplot grid
p_cor_grid <- ggplot(cor_df_long, aes(x = group, y = value, fill = group)) +
  geom_hline(yintercept = 0, linewidth=0.5) +
  geom_boxplot() +
  scale_fill_manual(values = cols_group) +
  facet_wrap(~factor) +
  labs(y="r (Krt8+ADI ~ factor)") +
  theme_bw() +
  theme(axis.text.x = element_blank(), legend.position = "bottom")

png(file.path(DIR_FIG_C2L, paste0("mm_visium_c2l_res_Krt8ADI_NMF_corr_boxplot.png")), width = 10*fig_res, height = 8*fig_res, res = fig_res)
p_cor_grid 
dev.off()

#' Dotplot
p_cor_p <- ggplot(cor_df_long, aes(y=reorder(factor, desc(factor)), x=reorder(sample_name, desc(sample_name)), color=value, size=abs(value))) +
  geom_point() +
  scale_color_gradient2(low = "steelblue", mid = "grey90", high = "orange") +
  scale_x_discrete(position = "top") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title = element_blank(),
        axis.text.x = element_text(angle=45, hjust=0, vjust=1));p_cor_p

png(file.path(DIR_FIG_C2L, paste0("mm_visium_c2l_res_Krt8ADI_NMF_corr_dotplot.png")), width = 8*fig_res, height = 4*fig_res, res = fig_res)
p_cor_p
dev.off()


##### Krt8.ADI-high spots and neighbors
cell_column_name <- "c2l_Krt8.ADI"
x_var <- sym(cell_column_name) # enquo() !!x_var

###### Find cut-off ###### 
####### Histogram #######
mdata_hist <- c2l_mdata %>% 
  select(sample_name, group, day, condition, matches(cell_column_name))

mean_vlines <- mdata_hist %>% 
  group_by(group) %>% 
  summarise(mean_density = mean(!!x_var),
            median_density = median(!!x_var))

ggplot(mdata_hist, aes(x = !!x_var, fill = group, color=group)) +
  geom_histogram(position = "identity", alpha = 0.4, bins = 60) +
  geom_vline(data = mean_vlines, aes(xintercept = median_density, color = group),
             linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = 2.5, linetype = "dotted") +
  scale_fill_manual(values = cols_group) +
  scale_color_manual(values = cols_group) +
  facet_wrap(~group) +
  scale_x_log10() +
  theme_bw()

c2l_mdata %>% 
  filter(!!x_var > 2.5) %>%
  group_by(group) %>% 
  summarise(
    n_spots = n(),
    sum_density = sum(!!x_var),
    max_density = max(!!x_var),
    mean_density = mean(!!x_var),
    median_density = median(!!x_var)) %>% 
  View()


###### Define high spots ###### 
krt8adi_cutoff <- 2.5
se$krt8adi_high <- ifelse(se$c2l_Krt8.ADI > krt8adi_cutoff, "high", "low")

ST.FeaturePlot(se, features = "krt8adi_high", ncol = 6, label.by = "sample_name", cols = c("grey90", "red"))

###### Neighbors ###### 
# Find nbs
se_krt8adi <- SubsetSTData(se, condition == "bleomycin")
se_krt8adi <- LoadImages(se_krt8adi, xdim = 10)

se_krt8adi <- SetIdent(se_krt8adi, value = "krt8adi_high")
se_krt8adi <- RegionNeighbours(se_krt8adi, id = "high")

# ST.FeaturePlot(se_krt8adi, features = "nbs_high", ncol = 6)

#' subcluster nbs spots
se_krt8adi_nbs <- subset(se_krt8adi, nbs_high == "nbs_high")

se_krt8adi_nbs <- se_krt8adi_nbs %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30, n.neighbors = 50)

se_krt8adi_nbs <- RunHarmony(se_krt8adi_nbs, group.by.vars = "day", 
                             assay.use = "SCT", plot_convergence = T, 
                             reduction.save = "pca_harmony")

se_krt8adi_nbs <- RunUMAP(se_krt8adi_nbs, 
                          reduction = "pca_harmony", 
                          dims = 1:15, 
                          reduction.name = "umap_harmony")


DimPlot(se_krt8adi_nbs, group.by = "group", cols = cols_group)
DimPlot(se_krt8adi_nbs, reduction = "umap_harmony", group.by = "group", cols = cols_group)


se_krt8adi_nbs <- FindNeighbors(se_krt8adi_nbs, reduction = "pca_harmony", dims = 1:15)
se_krt8adi_nbs <- FindClusters(se_krt8adi_nbs, resolution = 0.4)
se_krt8adi_nbs$krt8adi_nbs_clusters <- factor(paste0("NbC_", as.numeric(se_krt8adi_nbs$seurat_clusters)))
se_krt8adi_nbs <- SetIdent(se_krt8adi_nbs, value = "krt8adi_nbs_clusters")

saveRDS(se_krt8adi_nbs, file = file.path(DIR_RES, "objects", "mm_visium_se_obj_krt8adi_nbs.rds"))


# plot umap
p_umap <- DimPlot(se_krt8adi_nbs, 
        reduction = "umap_harmony", 
        group.by = "krt8adi_nbs_clusters", 
        pt.size = 0.5,
        cols = cols_d3_20) /
  DimPlot(se_krt8adi_nbs, 
          reduction = "umap_harmony", 
          pt.size = 0.5,
          group.by = "group", 
          cols = cols_group)
p_umap <- p_umap & theme(legend.position = "right", 
                         aspect.ratio = 1, 
                         axis.text = element_text(size=8),
                         axis.title = element_text(size=8),
                         legend.text = element_text(size=8),
                         plot.title = element_text(hjust=0.5, size=11, face = "plain"))

cells_plot <- c("c2l_Krt8.ADI", 
                "c2l_Myofibroblasts",
                "c2l_Activated.AT2.cells",
                "c2l_M2.macrophages",
                "c2l_AM..Bleo.",
                "c2l_Plasma.cells")
p_umap_cell <- FeaturePlot(se_krt8adi_nbs, 
                           features = cells_plot,
                           reduction = "umap_harmony", 
                           cols = col_scale_rocket, 
                           pt.size = 0.5, 
                           ncol = length(cells_plot)) & 
  theme_void() &
  theme(legend.position = "none", 
        aspect.ratio = 1, 
        plot.title = element_text(size=10, hjust=0.5))


# sample proportions
anno_df_sample <- se_krt8adi_nbs[[]] %>%
  group_by(sample_name, krt8adi_nbs_clusters) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n),
         pct = round(100 * n/sum(n), 0),
         rel.freq = paste0(round(100 * n/sum(n), 0), "%")) %>%
  as.data.frame() %>%
  arrange(sample_name, desc(freq))
anno_df_sample$sample_name <- factor(anno_df_sample$sample_name, levels = rev(unique(anno_df_sample$sample_name)))

anno_df_group <- se_krt8adi_nbs[[]] %>%
  group_by(group, krt8adi_nbs_clusters) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n),
         pct = round(100 * n/sum(n), 0),
         rel.freq = paste0(round(100 * n/sum(n), 0), "%")) %>%
  as.data.frame() %>%
  arrange(group, desc(freq))
anno_df_group$group <- factor(anno_df_group$group, levels = rev(unique(anno_df_group$group)))


p_prop_sample <- ggplot(anno_df_sample, aes(x=sample_name, y=freq, fill=krt8adi_nbs_clusters)) +
  geom_bar(stat = 'identity', colour=NA, position = "stack", width = 0.7) +
  scale_fill_manual(values = cols_d3_20) +
  labs(x="", y="Frequency", fill="Cluster", title="Cluster proportions in samples") +
  theme_linedraw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right", 
    text = element_text(size=8),
    plot.title = element_text(hjust=0.5)) +
  coord_flip()

p_prop_group <- ggplot(anno_df_group, aes(x=group, y=freq, fill=krt8adi_nbs_clusters)) +
  geom_bar(stat = 'identity', colour=NA, position = "dodge") +
  scale_fill_manual(values = cols_d3_20) +
  labs(x="", y="Frequency", fill="Cluster", title="Cluster proportions in groups") +
  theme_linedraw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none", 
    text = element_text(size=8),
    plot.title = element_text(hjust=0.5)) +
  coord_flip()

p_prop <- (p_prop_group / p_prop_sample) + patchwork::plot_layout(heights = c(1,2))


# marker genes
nbs_cluster_markers <- FindAllMarkers(se_krt8adi_nbs)

top_markers <- nbs_cluster_markers %>% 
  group_by(cluster) %>%
  filter(p_val_adj < 0.01 & avg_log2FC > 0) %>%
  # arrange(desc(avg_log2FC)) %>% 
  top_n(n = 12, wt = avg_log2FC) %>% 
  arrange(cluster)

p_marker_dp <- DotPlot(se_krt8adi_nbs, features = unique(top_markers$gene), cols="RdYlBu") +
  theme_dotplot +
  theme(axis.title = element_blank(), 
        legend.position = "none", 
        plot.title = element_text(face = "bold", size = 11, colour = "black"),
        axis.text.y = element_text(face = "plain", size = 10, colour = "black"),
        axis.text.x = element_text(angle=45, hjust=1, vjust=1, size = 8))


# save plots
pdf(file = file.path(DIR_FIG_C2L, "mm_visium_c2l_res_krt8adi_high_nbs_subclusters.pdf"), width = 11, height = 10, useDingbats = F)
wrap_plots(list(p_prop, p_umap), ncol = 2, widths = c(3,2)) / wrap_plots(list(p_umap_cell, p_marker_dp), ncol = 1)
dev.off()


###### Add to se blm metadata ######
# Add to blm subset
se_blm <- SubsetSTData(se, condition == "bleomycin")
se_blm$krt8adi_high_nbs <- ifelse(se_blm$krt8adi_high == "high", "krt8adi_high", "NA")
se_blm$krt8adi_high_nbs[rownames(se_krt8adi_nbs[[]])] <- as.character(se_krt8adi_nbs$krt8adi_nbs_clusters)
se_blm$krt8adi_high_nbs <- factor(se_blm$krt8adi_high_nbs)

ST.FeaturePlot(se_blm, features = "krt8adi_high_nbs", ncol = 6, 
               cols = c("black", "grey90", cols_d3_20), label.by = "sample_name") &
  theme(aspect.ratio = 1)

plot_colors <- setNames(c("black", "grey90", cols_d3_20), nm = levels(se_blm$krt8adi_high_nbs))

png(file.path(DIR_FIG_C2L, paste0("mm_visium_c2l_res_krt8adi_high_nbs_subclusters_cell_density_violin.png")), width = 20*fig_res, height = 10*fig_res, res = fig_res)
VlnPlot(se_blm, features = c2l_names, ncol = 10, group.by = "krt8adi_high_nbs", pt.size = 0, cols = plot_colors) & 
  theme_custom +
  theme(legend.position = "none", 
        plot.title = element_text(face = "plain", size = 10, colour = "black"),
        axis.text.y = element_text(face = "plain", size = 10, colour = "black"),
        axis.text.x = element_text(angle=45, hjust=1, vjust=1, size = 8))
dev.off()


# Plot cell density in clusters
se_blm_cell_cluster_plot <- se_blm@meta.data %>% 
  select(sample_name, group, day, krt8adi_high_nbs, matches(c2l_names)) %>% 
  filter(krt8adi_high_nbs != "NA")

se_blm_cell_cluster_plot_long <- pivot_longer(se_blm_cell_cluster_plot, 
                                              cols = starts_with("c2l_"),
                                              names_to = "cell_type",
                                              values_to = "cell_density")
se_blm_cell_cluster_plot_long$cell_type <- gsub("c2l_", "", gsub("[.]$", "", gsub("..", "_", se_blm_cell_cluster_plot_long$cell_type, fixed = TRUE)))
celltype_names <- unique(se_blm_cell_cluster_plot_long$cell_type)

celltype_names_filtered <- celltype_names[!celltype_names %in% c("NA", "low.quality.cells")]
se_blm_cell_cluster_plot_long <- se_blm_cell_cluster_plot_long %>% filter(cell_type %in% celltype_names_filtered)


c_type_order <- se_blm_cell_cluster_plot_long %>% 
  group_by(cell_type) %>% 
  summarise(mean_density = mean(cell_density)) %>% 
  arrange(desc(mean_density))

se_blm_cell_cluster_plot_long$cell_type <- factor(se_blm_cell_cluster_plot_long$cell_type, levels = as.character(c_type_order$cell_type))

# unique(se_blm_cell_cluster_plot_long$cell_type)
# dim(se_blm_cell_cluster_plot_long)
theme_density_grid <- theme(panel.grid = element_blank(),
                            panel.background = element_rect(fill = NA, color = NA),
                            panel.border = element_rect(fill = NA, color = NA),
                            legend.position = "none",
                            panel.spacing=unit(0, "lines"), 
                            axis.ticks = element_line(linewidth = 0.25, color = "black"), 
                            axis.line = element_blank(),
                            axis.text = element_text(size=6, color="black"),
                            strip.background = element_rect(fill=NA, color = NA),
                            strip.text = element_text(color="black"),
                            strip.text.y.right = element_text(angle = 0, hjust=0))

pdf(file = file.path(DIR_FIG_C2L, "mm_visium_c2l_res_krt8adi_high_nbs_subclusters_cell_density_grid_log10.pdf"), width = 8, height = 18, useDingbats = F)
ggplot(se_blm_cell_cluster_plot_long, aes(x=cell_density, fill=cell_type)) +
  geom_density(color=NA) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "black") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "black", linetype = "dashed") +
  scale_x_log10() +
  facet_grid(cell_type ~ krt8adi_high_nbs, scales = "free_y") +
  theme_bw() +
  theme_density_grid
dev.off()

pdf(file = file.path(DIR_FIG_C2L, "mm_visium_c2l_res_krt8adi_high_nbs_subclusters_cell_density_grid.pdf"), width = 8, height = 18, useDingbats = F)
ggplot(se_blm_cell_cluster_plot_long, aes(x=cell_density, fill=cell_type)) +
  geom_density(color=NA) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "black") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "black", linetype = "dashed") +
  facet_grid(cell_type ~ krt8adi_high_nbs, scales = "free_y") +
  theme_bw() +
  theme_density_grid
dev.off()


# plot on HE
sample_select <- "d21_bleo_4a"
# sample_select <- "d7_bleo_3"
se_blm_he <- SubsetSTData(se_blm, sample_name == sample_select & krt8adi_high_nbs != "NA")
se_blm_he <- LoadImages(se_blm_he, xdim = 1e3)

plot_colors <- setNames(c("black", NA, cols_d3_20), nm = levels(se_blm_he$krt8adi_high_nbs))

png(file.path(DIR_FIG_C2L, paste0("mm_visium_c2l_res_krt8adi_high_nbs_subclusters_", sample_select, ".png")), width = 10*fig_res, height = 10*fig_res, res = fig_res)
FeatureOverlay(se_blm_he, features = "krt8adi_high_nbs", 
               cols = plot_colors, pt.alpha = 0.6, pt.size = 2,
               label.by = "sample_name") &
  theme(aspect.ratio = 1)
dev.off()



