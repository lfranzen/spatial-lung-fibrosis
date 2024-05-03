#' [mm_visium_plots.R]
#'
#' Generation of other plots outside of main analyses 
#'
#'
#' Aug 2022, L. Franzén [lovisa.franzen@scilifelab.se]

#### Set up ####
##### Define params. ####
set.seed(1)
SPECIES <- "mouse"
DIR_ROOT <- "/home/st-analysis_home/lovisa.franzen/analysis/lung/spatial-lung-fibrosis"  #getwd()
DIR_DATA <- file.path(DIR_ROOT, "data", SPECIES, "visium")
DIR_RES <- file.path(DIR_ROOT, "results", SPECIES)
DIR_FIG <- file.path(DIR_RES, "figures")
fig_res <- 500

##### Load libs ####
library(STutility)
library(harmony)
library(RColorBrewer)
library(patchwork)
library(dplyr)
library(magrittr)


##### Other ####
# source(file.path(DIR_ROOT, "scripts", "colors.R"))
source(file.path(DIR_ROOT, "scripts", "custom_functions.R"))
source(file.path(DIR_ROOT, "scripts", "custom_colors.R"))
theme_custom <- theme(axis.title.x = element_blank())

##### Read objects ####
fname <- paste0("mm_visium_preproc_se_obj.rds")
se <- readRDS(file = file.path(DIR_RES, "objects", fname))

# Processs annotation labels
se$annotation2 <- ifelse(se$annotation == "Inflammation", paste0(se$annotation, "_", se$day), se$annotation)
se$annotation3 <- ifelse(se$annotation == "Inflammation" | se$annotation == "Suspect Fibrosis/Fibroplasia", 
                         "Diseased", se$annotation)


#### Plot ####

##### Histopath props ####
metadata <- se[[]]

anno_df <- metadata %>%
  group_by(group, sample_name, annotation) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n),
         pct = round(100 * n/sum(n), 0),
         rel.freq = paste0(round(100 * n/sum(n), 0), "%")) %>%
  as.data.frame() %>%
  arrange(annotation, desc(freq))

# anno_df$sample_name <- factor(anno_df$sample_name, levels = c())
# subset(anno_df, annotation=="Normal Alveolar and Other")
# anno_df$rank <- 1:nrow(anno_df)

p <- ggplot(subset(anno_df, !is.na(annotation))) +
  geom_segment(mapping = aes(x=sample_name, xend=sample_name, y=0, yend=pct), color="black") +
  geom_point(mapping = aes(x=sample_name, y=pct, color=group)) +
  geom_hline(yintercept = 0, color="black", size=.2) +
  ylim(0,100) + 
  labs(y="% of total spots", x="") +
  # geom_bar(position="dodge", stat="identity") %>%
  scale_color_manual(values = cols_group) +
  facet_grid(~annotation, scales = "free_y", space='free') +
  coord_flip() +
  theme_linedraw() +
  theme(
    panel.grid.major.x = element_line(color="black"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(), 
    panel.grid.minor.y = element_blank(),
    legend.position = "none")

png(filename = file.path(DIR_FIG, "mm_visium_plots_histopath_lollipop.png"), 
    width = 12*fig_res, height = 3*fig_res, res = fig_res);p;dev.off()


##### Histopath spatial ####
# Spatial clusters
se_annotated <- SubsetSTData(se, annotation != "NA")

p_spatial <- ST.FeaturePlot(se_annotated, 
                            features = "annotation",
                            label.by = "sample_name",
                            ncol = 3, 
                            pt.size = 0.7, 
                            pt.border = F, 
                            cols = cols_annotation,
                            show.sb = F) &
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5), legend.position = "top") &
  guides(fill = guide_legend(override.aes = list(size = 2), ncol = 2))

png(file = file.path(DIR_RES, "figures", "mm_visium_histopath_spatial.png"), width = 5*fig_res, height = 12*fig_res, res = fig_res);p_spatial;dev.off()



##### Fig Ext Data 4 b-c: Stats ####
#' 1: n spots per animal
#' 2: % annotation as diseased tissue (inflammation+fibrosis) per animal
#' 3: pseudo-bulk DEA per t.p
#' 4. qPCR overlap with DEA results (optional - excluded)

# Read data for 1+2
metadata <- se[[]]
anno_df_add <- metadata %>%
  group_by(animal) %>% 
  summarise(.groups = "keep", n_spots = n())

anno_df_plot <- metadata %>%
  group_by(condition, day, group, animal, annotation3) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n),
         pct = round(100 * n/sum(n), 0),
         rel.freq = paste0(round(100 * n/sum(n), 0), "%")) %>%
  as.data.frame() %>%
  arrange(annotation3, desc(freq)) %>% 
  filter(!is.na(annotation3))
anno_df_plot <- merge(anno_df_plot, anno_df_add, by="animal")
colnames(anno_df_plot)[colnames(anno_df_plot)=="annotation3"] <- "region"

write.csv(anno_df_plot, 
          file.path(DIR_RES, "objects", "mm_visium_preproc_animal_summary_histopath_fig4b.csv"), 
          row.names = F)
anno_df_plot <- read.csv(file.path(DIR_RES, "objects", "mm_visium_preproc_animal_summary_histopath_fig4b.csv"))

sum(subset(anno_df_plot, day=="d7")[,"n_spots"])
sum(subset(anno_df_plot, day=="d21")[,"n_spots"])
sum(subset(anno_df_plot, day=="d7" & condition=="bleomycin")[,"n_spots"])
sum(subset(anno_df_plot, day=="d21" & condition=="bleomycin")[,"n_spots"])

# Plot n spots per sample
anno_df_plot$condition <- factor(anno_df_plot$condition, levels = c("bleomycin", "control"))
anno_df_plot$day <- factor(anno_df_plot$day, levels = c("d7", "d21"))

p1 <- ggplot(anno_df_plot[!duplicated(anno_df_plot$animal),], aes(x="spots", y=n_spots, fill=condition)) +
  geom_bar(position="fill", stat="identity", color="black") +
  facet_wrap(~day) +
  scale_fill_manual(values = cols_cond) +
  coord_flip() +
  theme_minimal() +
  theme(axis.title = element_blank(), 
        axis.text.x = element_blank(),
        # axis.text = element_blank(),
        panel.grid = element_blank());p1

p2 <- ggplot(subset(anno_df_plot), aes(x=region, y=pct, fill=condition)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~day) +
  scale_fill_manual(values = cols_cond) +
  ylim(0,100) +
  coord_flip() +
  theme_minimal() +
  theme(axis.title = element_blank(), 
        strip.text = element_blank(),
        # axis.text = element_blank(),
        panel.grid = element_blank());p2

p1/p2 + plot_layout(heights = c(1,4))


#' DEA volcano
dea_d7 <- read.csv(file.path(DIR_RES, "objects", "DEA", "mm_visium_pseudobulk_anno_all_day7_deseq2_res.csv"), row.names = 1)
dea_d21 <- read.csv(file.path(DIR_RES, "objects", "DEA", "mm_visium_pseudobulk_anno_all_day21_deseq2_res.csv"), row.names = 1)

dea_d7$day <- "d7"
dea_d21$day <- "d21"
dea_d7$gene <- rownames(dea_d7)
dea_d21$gene <- rownames(dea_d21)
dea_res <- rbind(dea_d7, dea_d21)
dea_res <- dea_res %>% filter(padj<0.01)

dea_res$day <- factor(dea_res$day, levels = c("d7", "d21"))
dea_res$fc <- 2^dea_res$log2FoldChange
dea_res$log10FC <- log10(dea_res$fc)

#' qPCR
qpcr <- read.csv(file.path(DIR_DATA, "..", "qpcr", "Log10FC_heatmap.csv"), sep = ";") # Optional
colnames(qpcr)[1] <- "gene"
rownames(qpcr) <- qpcr$gene

# qpcr_dea_res_add_d7 <- subset(dea_res, day=="d7")[unique(qpcr$gene), ] %>% filter(!is.na(gene))
# qpcr_dea_res_add_d21 <- subset(dea_res, day=="d21")[unique(qpcr$gene), ] %>% filter(!is.na(gene))
# 
# qpcr$visium_d7 <- ifelse(rownames(qpcr) %in% qpcr_dea_res_add_d7$gene, qpcr_dea_res_add_d7$log2FoldChange, 0)
# qpcr$visium_d21 <- ifelse(rownames(qpcr) %in% qpcr_dea_res_add_d21$gene, qpcr_dea_res_add_d21$log2FoldChange, 0)
# 
# 
# qpcr_long <- pivot_longer(qpcr, cols = grep("^d[1-9]", colnames(qpcr), value = T), names_to = "animal", values_to = "log10FC")
# qpcr_long$fc <- 10^qpcr_long$log10FC
# qpcr_long$log2FoldChange <- log2(qpcr_long$fc)
# 
# qpcr_long$animal <- gsub("ctr_", "c", qpcr_long$animal)
# qpcr_long$animal <- gsub("bleo_", "b", qpcr_long$animal)
# qpcr_long <- merge(qpcr_long, mdat, by="animal")
# qpcr_long$day <- factor(qpcr_long$day, levels = c("d7", "d21"))
# qpcr_long$animal2 <- gsub("d7_|d21_", "", qpcr_long$animal)
# 
# qpcr_long$visium_log2fc <- ifelse(qpcr_long$day=="d7", qpcr_long$visium_d7, qpcr_long$visium_d21)
# qpcr_long$visium_dir <- ifelse(qpcr_long$visium_log2fc>0, "pos", "neg")
# qpcr_long$qpcr_dir <- ifelse(qpcr_long$log2FoldChange>0, "pos", "neg")

# ifelse(qpcr_dea_res_add$gene %in% qpcr$gene, qpcr$gene, qpcr$gene)

# ylims <- max(qpcr_long$visium_log2fc)
# 
# p4 <- ggplot(qpcr_long, aes(x=animal2, y=visium_log2fc, color=log2FoldChange)) +
#   geom_point(shape="|", size=3) +
#   facet_wrap(~day) +
#   labs(y="log2 FC") +
#   # scale_color_manual(values = cols_cond) +
#   # scale_color_manual(values = c("lightblue", "hotpink")) +
#   scale_color_gradientn(
#   # colors = col_scale_div_custom2,
#     colors = hcl.colors(9, palette = "Tropic"),
#     limits = c(-ylims, ylims)) +
#   ylim(-ylims,ylims) +
#   coord_flip() +
#   geom_hline(yintercept = 0, linewidth=0.1) +
#   theme_minimal() +
#   theme(axis.title.y = element_blank(), 
#         # strip.text = element_blank(),
#         # axis.text = element_blank(),
#         panel.grid = element_blank(),
#         legend.position = "bottom");p4
# 
# (p1/p2/p4) + plot_layout(heights = c(1,4,2.5))


#' DEA volcano + qpcr
# dea_res$qpcr_gene <- ifelse(dea_res$gene %in% qpcr$gene, "qpcr", "no") # Optional
dea_res$direction <- ifelse(dea_res$log2FoldChange>0, "bleomycin", "control")

top_genes <- dea_res %>% 
  group_by(day, direction) %>% 
  slice_max(n = 4, order_by = abs(log2FoldChange)) %>% 
  as.data.frame() %>% 
  select(gene)
top_genes2 <- dea_res %>% 
  group_by(day, direction) %>% 
  slice_max(n = 4, order_by = -padj) %>% 
  as.data.frame() %>% 
  select(gene)
top_genes <- unique(c(top_genes$gene, top_genes2$gene))

dea_res$top_genes <- ifelse(dea_res$gene %in% top_genes, dea_res$gene, "x")

cols_qpcr <- setNames(c("black", "grey"),
                      nm = c("qpcr", "no"))
xylims <- max(dea_res$log2FoldChange)

p3 <- ggplot(dea_res, aes(x = log2FoldChange, y = -log10(padj), color = direction)) +
  geom_vline(xintercept = 0, linetype="solid", color = "black", linewidth=0.1) +
  geom_point(shape=16, size = 1, alpha = 0.5) +
  # ggrepel::geom_text_repel(data = subset(dea_res, qpcr_gene == "qpcr"), 
  #                          mapping = aes(x = log2FoldChange, y = -log10(padj), label = gene), 
  #                          color = "black", 
  #                          size = 3, 
  #                          max.overlaps = 50) +
  ggrepel::geom_text_repel(data = subset(dea_res, top_genes != "x"),
                           mapping = aes(x = log2FoldChange, y = -log10(padj), label = gene),
                           color = "black",
                           size = 3,
                           max.overlaps = 50) +
  scale_color_manual(values = cols_cond) +
  xlim(c(-xylims, xylims)) +
  facet_wrap(~day) +
  labs(x="log2 FC") +
  theme_linedraw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        plot.title = element_text(size = rel(1), hjust = 0.5),
        axis.title = element_text(size = rel(1)),
        axis.line = element_line(colour = "black"));p3

# Export images
p_up <- (p1/p2) + plot_layout(heights = c(1,4)) & theme(legend.position = "none")

pdf(file = file.path(DIR_FIG, "mm_visium_plots_fig4b_p1_stats.pdf"), width = 5, height = 2) # Ext Data Fig 4b
p_up
dev.off()

pdf(file = file.path(DIR_FIG, "mm_visium_plots_fig4b_p2_dea.pdf"), width = 5, height = 1.5) # Ext Data Fig 4c
p3
dev.off()

# pdf(file = file.path(DIR_FIG, "mm_visium_plots_fig4b_p2_dea_large.pdf"), width = 10, height = 6)
# ggplot(dea_res, aes(x = log2FoldChange, y = -log10(padj), color = direction)) +
#   geom_vline(xintercept = 0, linetype="solid", color = "black", linewidth=0.1) +
#   geom_point(shape=16, size = 1, alpha = 0.5) +
#   geom_point(data = subset(dea_res, top_genes != "x"),
#              mapping = aes(x = log2FoldChange, y = -log10(padj)),
#              shape=18, size = 1, color="black") +
#   ggrepel::geom_text_repel(data = subset(dea_res, top_genes != "x"),
#                            mapping = aes(x = log2FoldChange, y = -log10(padj), label = gene),
#                            color = "black",
#                            size = 3,
#                            max.overlaps = 50) +
#   scale_color_manual(values = cols_cond) +
#   xlim(c(-xylims, xylims)) +
#   facet_wrap(~day) +
#   labs(x="log2 FC") +
#   theme_linedraw() +
#   theme(legend.position = "none",
#         panel.grid = element_blank(),
#         plot.title = element_text(size = rel(1), hjust = 0.5),
#         axis.title = element_text(size = rel(1)),
#         axis.line = element_line(colour = "black"))
# dev.off()
