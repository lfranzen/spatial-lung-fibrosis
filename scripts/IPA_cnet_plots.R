# Load necessary libraries
library(multienrichjam)
library(jamba)
library(colorjam)
suppressPackageStartupMessages(library(ComplexHeatmap))
options("stringsAsFactors" = FALSE, "warn" = -1)
knitr::opts_chunk$set(
  fig.height = 10,
  fig.width = 10,
  fig.align = "center"
)
ragg_png = function(..., res = 192) {
  ragg::agg_png(..., res = res, units = "in")
}
knitr::opts_chunk$set(dev = "ragg_png", fig.ext = "png")

# Import and use IPA enrichment data
hs_F14C0 <- "human_f14_hi_C0_markers_vs_all_IPF_only.txt" #IPA output file
mm_F14C0 <- "mouse_f14_hi_C0_markers.txt" #IPA output file

# Import IPA data from text files
hs_F14C0_dfl <- importIPAenrichment(hs_F14C0)
mm_F14C0_dfl <- importIPAenrichment(mm_F14C0)

# Combine lists into one
ipa_l <- list(hs_F14C0 = hs_F14C0_dfl, mm_F14C0 = mm_F14C0_dfl)

# Check the dimensions within each list
ssdim(ipa_l)

# Take only the Ingenuity Canonical Pathways
enrichList_upstream <- lapply(ipa_l, function(i) {
  i[["Upstream Regulators"]]
})
sdim(enrichList_upstream)

# Convert data.frame to enrichResult
er_upstream <- lapply(enrichList_upstream, function(i) {
  enrichDF2enrichResult(
    i,
    keyColname = "Name",
    pvalueColname = "P-value",
    geneColname = "geneNames",
    pvalueCutoff = 0.05 #replace as necessary
  )
})

kable_coloring(
  head(as.data.frame(er_upstream[[1]])),
  caption = "Top 10 rows of enrichment data",
  row.names = FALSE
) %>%
  kableExtra::column_spec(
    column = seq_len(ncol(er_upstream[[1]])),
    border_left = "1px solid #DDDDDD",
    extra_css = "white-space: nowrap;"
  )

# Perform multi-enrichment mapping
mem_upstream <- multiEnrichMap(
  er_upstream,
  enrichBaseline = 1,
  cutoffRowMinP = 0.05,  #replace as necessary
  colorV = c("#B26694", "#517595"),
  topEnrichN = 20 #replace as necessary
)

# Display sdim of mem_upstream
kable_coloring(
  sdim(mem_upstream),
  caption = "sdim(mem_upstream)"
) %>%
  kableExtra::column_spec(
    column = seq_len(4),
    border_left = "1px solid #DDDDDD",
    extra_css = "white-space: nowrap;"
  )


dev.off()

# Create a PDF file for the upstream regulators
pdf(
  file = "hsF14C0_vs_mmF14C0_cnetUpstream_findMarkers.pdf",
  height = 10,
  width = 10
)

# Generate the cnet plots
mem_upstream_plots <- multienrichjam::mem_plot_folio(
  mem_upstream,
  pathway_column_split = 4,
  node_factor = 2,
  use_shadowText = TRUE,
  label_factor = 1.1,
  do_which = 5,
  verbose = TRUE,
  main = "Canonical Pathways"
)

dev.off()


# Extract and save selected columns for hs_F14C0
hs_obj <- mem_upstream[["enrichList"]][["hs_F14C0"]]
hs_result_df <- hs_obj@result
hs_selectedCols <- hs_result_df[, c("ID", "Upstream Regulator", "Molecule Type", "Predicted Activation State", "Activation z-score", "Bias-corrected z-score", "pvalue", "geneID", "Mechanistic Network", "p.adjust", "Count")]
write.csv(hs_selectedCols, file = "hs_us_findMarkers_posAndNeg_4clust_selectedCols.csv", row.names = FALSE)

# Extract and save selected columns for mm_F14C0
mm_obj <- mem_upstream[["enrichList"]][["mm_F14C0"]]
mm_result_df <- mm_obj@result
mm_selectedCols <- mm_result_df[, c("ID", "Upstream Regulator", "Molecule Type", "Predicted Activation State", "Activation z-score", "Bias-corrected z-score", "pvalue", "geneID", "Mechanistic Network", "p.adjust", "Count")]
write.csv(mm_selectedCols, file = "mm_us_findMarkers_posAndNeg_4clust_selectedCols.csv", row.names = FALSE)

# Extract and save cluster data
clust_obj <- mem_upstream_plots[["clusters_mem"]]

clust_df <- data.frame(
  A = c(clust_obj$A, rep(NA, max(lengths(clust_obj)) - length(clust_obj$A))),
  B = c(clust_obj$B, rep(NA, max(lengths(clust_obj)) - length(clust_obj$B))),
  C = c(clust_obj$C, rep(NA, max(lengths(clust_obj)) - length(clust_obj$C))),
  D = c(clust_obj$D, rep(NA, max(lengths(clust_obj)) - length(clust_obj$D)))
)
write.csv(clust_df, "clust_us_findMarkers_posAndNeg_4clust_output.csv", row.names = FALSE)