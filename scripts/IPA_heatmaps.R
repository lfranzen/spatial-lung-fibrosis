library(pheatmap)

# Custom breaks for heatmap
breaks_custom <- seq(-6, 6, length.out = 14)

# IPA comparison output (Upstream or pathways or diseases)
IPA_us <- read.table("human_F14hiC0_vs_allSubclusters_upstream_comparison.txt", sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE, na.strings = "N/A")

# Extract the specified columns
selected_cols <- IPA_us[, c(1, 8, 7, 6, 5)]
selected_rows <- selected_cols[1:20, ]
write.table(selected_rows, file = "selected_us_data.txt", sep = "\t", quote = FALSE, row.names = FALSE)

upstream <- read.table("selected_us_data.txt", sep = "\t", header = TRUE, row.names = 1)
data <- as.matrix(upstream)

# Create heatmap
  pheatmap(
    data,
    border_color = "gray100",
    na_col = "grey70",
    fontsize_row = 14,
    breaks = breaks_custom,
    fontsize_col = 11,
    filename = "IPA_us_scale.pdf",
    width = 10,
    height = 15,
    angle_col = 45,
    main = "Upstream Regulators",
    display_numbers = FALSE,
    color = col_scale_div_custom2,
    cellwidth = 20,
    cellheight = 20,
    cluster_cols = FALSE,
    cluster_rows = FALSE
  )
