#' [custom_functions.R]
#'
#' Collection of custom functions used within the project
#'
#' L. Franzén [lovisa.franzen@scilifelab.se]


#' Create an InfoTable of paths and metadata for STutility
#' 
#' @param data.dir Parent path to spaceranger output
#' @param sample.select.ids List of samples to include corrsponding to data folder ids
#' @param metadata.df A dataframe with sample metadata to be added to STutility object
#' @param use.raw Use raw unfiltered spot data (default: FALSE)
#' @param subfolder.name Specify name of parent subfolder if needed (default: NA)
#' @return an InfoTable data frame
#' 
#' @export
#' 
createInfoTable <- function(data.dir = DIR_DATA, 
                            sample.select.ids, 
                            metadata.df, 
                            use.raw = FALSE,
                            subfolder.name=NA) {
  sample_dirs_all <- c()
  for (s in sample.select.ids){
    message("Fetching data path for sample ", s)
    sample_dirs <- grep(pattern = s, x = list.dirs(path = data.dir, recursive = F, full.names = T), value = T)
    sample_dirs_all <- c(sample_dirs_all, sample_dirs)
  }
  
  if(!is.na(subfolder.name)){
    sample_dirs_all <- paste0(sample_dirs_all, subfolder.name)
  }
  
  if(use.raw){
    data.dirs <- paste0(sample_dirs_all, "/raw_feature_bc_matrix.h5")
  } else {
    data.dirs <- paste0(sample_dirs_all, "/filtered_feature_bc_matrix.h5")
  }
  spot_dirs <- paste0(sample_dirs_all, "/spatial/tissue_positions_list.csv")
  img_dirs <- paste0(sample_dirs_all, "/spatial/tissue_hires_image.png")
  img_scale_dirs <- paste0(sample_dirs_all, "/spatial/scalefactors_json.json")
  
  infoTable <- data.frame(samples = data.dirs,
                          spotfiles = spot_dirs,
                          imgs = img_dirs,
                          json = img_scale_dirs,
                          metadata.df, 
                          stringsAsFactors = F)
  return(infoTable)
}


#' Plot QC stats using STutility
#' 
#' @param seurat.object 
#' @param color.group 
#' @param fill.by 
#' @param label.by 
#' @param extra.feat.plot 
#' @return plot object
#' 
#' @export
#' 
plotVisiumQC <- function(seurat.object, 
                         color.group = NA, 
                         fill.by = "group", 
                         label.by = "sample_name",
                         extra.feat.plot = NA) {
  require(STutility)
  require(ggplot2)
  
  theme_custom <- theme(axis.title.x = element_blank())
  if(is.na(color.group[1])){
    color.group <- RColorBrewer::brewer.pal(8, "Set2")
  }
  se_staffli <- GetStaffli(seurat.object)@meta.data
  p_dat <- cbind(seurat.object@meta.data, se_staffli)
  
  p1 <- ggplot(p_dat, aes_string(x="nFeature_RNA")) +
    geom_histogram(aes(y=..density..), colour=NA, fill="black", bins = 75)+
    geom_density(data = p_dat, aes_string(x="nFeature_RNA", y="..density..", fill = fill.by), alpha=.5) +
    geom_vline(aes(xintercept=mean(nFeature_RNA)), color="orange", linetype="dashed", linewidth=1) +
    scale_fill_manual(values = color.group) +
    ggtitle("Unique genes per spot") +
    theme_classic() + 
    theme(legend.position = "bottom", 
          axis.text = element_text(colour = "black"), 
          plot.title = element_text(hjust=0.5, face = "bold"))
  
  p2 <- ggplot(p_dat, aes_string(x="nCount_RNA")) +
    geom_histogram(aes(y=..density..), colour=NA, fill="black", bins = 75)+
    geom_density(data = p_dat, aes_string(x="nCount_RNA", y="..density..", fill = fill.by), alpha=.5) +
    geom_vline(aes(xintercept=mean(nFeature_RNA)), color="orange", linetype="dashed", linewidth=1) +
    scale_fill_manual(values = color.group) +
    scale_x_log10() +
    ggtitle("Total counts per spots") +
    theme_classic() + 
    theme(legend.position = "bottom", 
          axis.text = element_text(colour = "black"), 
          plot.title = element_text(hjust=0.5, face = "bold"))
  
  p3 <- VlnPlot(seurat.object, 
                features = c("nFeature_RNA", "nCount_RNA"), 
                pt.size=0, 
                group.by = label.by, 
                split.by = fill.by, 
                cols = color.group) & theme_custom
  
  if(!is.na(extra.feat.plot[1])){
    p4 <- VlnPlot(seurat.object, 
                  features = extra.feat.plot, 
                  pt.size=0, 
                  group.by = label.by, 
                  split.by = fill.by, 
                  cols = color.group) & theme_custom
    } else {
      p4=NULL
    }
  
  p <- (p1+p2)/p3/p4
  return(p)
}


#' Plot NMF spatial and gene loading
#' 
#' @param seurat.object 
#' @param factor 
#' @param col.scale 
#' @param labels 
#' @param indices 
#' @return plot object
#' 
#' @export
#' 
FactorPlot <- function(seurat.object, 
                       factor = 1, 
                       col.scale = c("darkblue", "cyan", "yellow", "red", "darkred"),
                       labels = "sample_name",
                       indices = c(1,2),
                       n.columns = NULL
) {
  # Check if NMF has been computed
  if (!"NMF" %in% names(seurat.object@reductions)) stop("NMF has not been computed ... \n", call. = FALSE)
  if(missing(n.columns)){
    n.columns <- length(indices)
  }
  
  p1 <- ST.DimPlot(seurat.object,
                   dims = factor, 
                   indices = indices,
                   label.by = labels,
                   center.zero = F, 
                   reduction = "NMF", 
                   ncol = n.columns, 
                   pt.size = 1, 
                   cols = col.scale, 
                   show.sb = F) & 
    labs(fill = "factor\nactivity") &
    theme(aspect.ratio = 1)
  p2 <- FactorGeneLoadingPlot(seurat.object, 
                              factor = factor, 
                              dark.theme = F) + 
    labs(title = "") +
    theme(axis.title.y = element_blank()) & labs(y = "weight")
  cowplot::plot_grid(p1, p2, ncol = 2, rel_widths = c(4, 1))
}


#' Plot NMF gene loadings horizontally
#' 
#' @param seurat.object 
#' @param factor.n 
#' @param topn 
#' @return plot object
#' 
#' @export
#' 
FactorGeneLoadingPlotHorizontal <- function (
    seurat.object,
    factor = 1,
    topn = 20
) {
  # Check if NMF has been computed
  if (!"NMF" %in% names(seurat.object@reductions)) stop("NMF has not been computed ... \n", call. = FALSE)
  
  ftr <- paste0("factor_", factor)
  nmf <- seurat.object@reductions$NMF@feature.loadings[, ftr]
  gene <- names(nmf)
  df <- data.frame(gene, val = nmf, stringsAsFactors = F)
  df <- df[order(df$val, decreasing = T), ]
  df <- df[1:topn, ]
  df$gene <- factor(df$gene, levels = df$gene)
  # CI95 <- mean(nmf) + (1.960 * sd(nmf) / sqrt(length(nmf)))
  
  p <- ggplot(df[1:topn, ], aes(reorder(gene, -val), val)) +
    geom_bar(stat = "identity", 
             fill = "grey20", 
             color = NA,
             width = 0.7) +
    labs(x = "", y = "value") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  # if(add.CI95.hline) p <- p + geom_hline(yintercept = CI95, color="orange")
  
  return(p)
}


#' Plot NMF factor spot ranking
#' 
#' @param seurat.object 
#' @param factor 
#' @param top.n.spots
#' @param include.zoom
#' @param split.by
#' @param group.by
#' @param color.group
#' @param line.alpha
#' @param plot.title
#' @return plot object
#' 
#' @export
#' 
FactorRankSpotPlot <- function(seurat.object,
                               factor = 1,
                               top.n.spots = 500,
                               include.zoom = TRUE,
                               split.by = "sample_name",
                               group.by = "group",
                               color.group = NA,
                               line.alpha = 1,
                               plot.title = NA){
  require(ggplot2)
  require(patchwork)
  
  factor_x <- factor
  factor_plot <- as.numeric(gsub("factor_", "", factor_x))
  rank_dat <- cbind(seurat.object@reductions$NMF@cell.embeddings[,factor_x], seurat.object@meta.data)
  colnames(rank_dat)[1] <- "nmf_factor"
  rank_dat$group_data_by <- rank_dat[ ,split.by]
  
  if(is.na(color.group[1])){
    n_groups <- length(unique(rank_dat[ ,group.by]))
    if(n_groups<3){
      n_cols <- 4
      color.group <- RColorBrewer::brewer.pal(n_cols, "Spectral")[c(1,4)]
    } else {
      color.group <- RColorBrewer::brewer.pal(n_groups, "Spectral")
    }
  }
  if(is.na(plot.title[1])){
    plot.title <- paste("Factor", factor_x)
  }
  
  rank_dat_plot <- rank_dat %>%
    dplyr::group_by(group_data_by) %>%
    dplyr::slice_max(order_by = nmf_factor, n = top.n.spots) %>%
    mutate(rank = dense_rank(-nmf_factor)) %>% 
    mutate_at(.vars = c(split.by, group.by), ~ as.factor(.))
  
  theme_plot <- theme(legend.justification = "center", 
                      plot.title = element_blank(), 
                      panel.grid = element_blank(), 
                      text = element_text(color = "black"), 
                      panel.border = element_rect(colour = "black"))
  p <- ggplot(rank_dat_plot, aes_string(x="rank", y="nmf_factor", group=split.by, color=group.by)) +
    geom_line(size=.5, alpha = line.alpha) +
    labs(y="", x="rank", color="") +
    scale_color_manual(values = color.group) +
    theme_bw() +
    theme_plot +
    theme(legend.position = "right")
  
  if(include.zoom){
    p_zoom <- ggplot(subset(rank_dat_plot, rank<11), aes_string(x="rank", y="nmf_factor", group=split.by, color=group.by)) +
      geom_line(linewidth=.5) +
      labs(y="factor activity", x="rank") +
      scale_color_manual(values = color.group) +
      scale_x_continuous(breaks = c(1,5,10)) +
      theme_bw() +
      theme_plot +
      theme(legend.position = "none")
    p <- p + theme(axis.title.y = element_blank())
    p_rank <- p_zoom + p + plot_layout(nrow = 1, widths = c(1,3))
  } else {
    p_zoom <- NULL
    p_rank <- p & labs(y="factor activity")
  }
  p_rank <- p_rank + plot_annotation(title = plot.title)
  return(p_rank)
}


#' Identify components of connected spots based on NMF factor embedding
#' 
#' @param seurat.object 
#' @param nmf.factor 
#' @param factor.embedding.cutoff
#' @param min.component.size
#' @param output.column.name
#' @return seurat.object
#' 
#' @export
#' 
IdentifyComponentsFromFactor <- function (
    seurat.object,
    nmf.factor = "factor_1",
    factor.embedding.cutoff = 0.1,
    min.component.size = 2,
    output.column.name = "graph_comp"
) {
  require(igraph)
  require(STutility)
  
  #' Prep
  if(sum(colnames(seurat.object@meta.data) %in% output.column.name)>0){
    seurat.object[[output.column.name]] <- NULL
  }
  
  se_nmf_emb <- seurat.object@reductions$NMF@cell.embeddings
  spots_f_cutoff <- intersect(colnames(seurat.object), 
                              rownames(se_nmf_emb[se_nmf_emb[, nmf.factor] > factor.embedding.cutoff, ])
  )
  se_f <- SubsetSTData(seurat.object, spots = spots_f_cutoff) # expression = factor_7>factor_cutoff
  
  #' Generate graph
  spatnet_init <- as.data.frame(GetSpatNet(se_f))
  
  el <- as.matrix(spatnet_init[,c("from", "to")])
  g <- graph_from_edgelist(el, directed = F)
  g <- simplify(g)  # remove double edges between nodes
  c <- components(g)
  comps <- c$membership
  
  #' Trim graph
  v_keep <- V(g)[comps %in% which(table(c$membership) >= min.component.size)]
  g_trim <- induced_subgraph(g, v_keep)
  
  spots_keep <- v_keep$name
  spot_comps <- comps[spots_keep]
  
  #' Add new labels to seurat metadata
  seurat.object$graph_comp <- NA
  seurat.object@meta.data[spots_keep, "graph_comp"] <- as.character(spot_comps)
  colnames(seurat.object@meta.data)[colnames(seurat.object@meta.data)=="graph_comp"] <- output.column.name
  
  message(paste0("Storing results in metadata column ", output.column.name))
  return(seurat.object)
}


#' Identify n nearest neighbors from a given center of origin
#' 
#' @param seurat.object 
#' @param feature.column.name 
#' @param center.feature
#' @param n.neighbors
#' @return seurat.object
#' 
#' @export
#' 
IdentifyNNeighbors <- function (
    seurat.object,
    feature.column.name,
    center.feature,
    n.neighbors = 4
) 
{
  message(paste0("Identifying neighbors for selected center ", center.feature))
  seurat.object$temp <- seurat.object@meta.data[, feature.column.name]
  seurat.object@meta.data[seurat.object$temp %in% center.feature, "temp"] <- "x"
  
  for(i in 1:n.neighbors){
    message(paste0("Finding spots of distance ", i))
    if(i==1){
      seurat.object <- SetIdent(seurat.object, value = "temp")
      seurat.object <- RegionNeighbours(seurat.object, id = "x", keep.within.id = T, verbose = TRUE)
      seurat.object$d <- 999
      seurat.object@meta.data[seurat.object$temp %in% "x", "d"] <- 0
      seurat.object@meta.data[rownames(subset(seurat.object@meta.data, nbs_x=="nbs_x")), "d"] <- i
    } else {
      seurat.object$temp <- "NA"
      seurat.object@meta.data[seurat.object$d <= i, "temp"] <- "x"
      seurat.object <- SetIdent(seurat.object, value = "temp")
      seurat.object <- RegionNeighbours(seurat.object, id = "x", keep.within.id = T, verbose = TRUE)
      seurat.object@meta.data[rownames(subset(seurat.object[[]], nbs_x=="nbs_x")), "d"] <- i
    }
  }
  seurat.object$d <- as.numeric(seurat.object$d)
  seurat.object@meta.data[seurat.object$d==999, "d"] <- NA
  
  colnames(seurat.object@meta.data)[colnames(seurat.object@meta.data)=="d"] <- paste0("d_c", center.feature)
  seurat.object$temp <- NULL
  seurat.object$nbs_x <- NULL
  
  return(seurat.object)
}



Scale01 <- function(x) {
  z <- (x - min(x)) / (max(x) - min(x))
  return(z)
}

ZScoreMartrix <- function(df, by_row = T){
  if (by_row){
    mar = 1
  } else {
    mar = 2
  }
  m = apply(df, mar, mean, na.rm = T)
  s = apply(df, mar, sd, na.rm = T)
  return((df - m) / s)
}

#### ggplot custom themes ####
theme_custom <- theme(axis.title.x = element_blank())

theme_dotplot <- theme( 
  plot.title = element_text(face = "bold", size = 10, colour = "black"),
  # axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
  panel.background = element_rect(fill = NA, color = "black", size = .5),
  panel.grid = element_blank(),
  axis.title = element_blank(), 
  axis.ticks = element_line(color = "black", linewidth = .5),
  axis.line = element_blank(),
  axis.text.x = element_text(size = 10, colour = "black"),
  axis.text.y = element_text(face = "italic", size = 10, colour = "black"),
  legend.text = element_text(size = 10, colour = "black"),
  legend.title = element_text(size = 10, colour = "black"),
  legend.position = "right")

theme_umap <- theme( 
  panel.background = element_rect(fill = NA, color = "black", size = .5),
  plot.title = element_text(face = "bold", size = 10, colour = "black"),
  axis.title= element_text(size = 10, colour = "black"),
  axis.ticks = element_blank(),
  axis.line = element_blank(),
  axis.text = element_blank(),
  legend.text = element_text(size = 10, colour = "black"),
  legend.title = element_text(size = 10, colour = "black"))
