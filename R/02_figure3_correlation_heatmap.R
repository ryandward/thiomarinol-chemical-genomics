# Figure 3: Pearson correlation between the chemical-gene S-score profiles of the 36 treatments,
# drawn as a heatmap with average-linkage dendrograms on both margins.

source("R/00_common.R")

library(pheatmap)
library(viridis)
library(ggplot2)
library(svglite)

median_chem <- load_median_s_scores()
score_matrix <- as_score_matrix(median_chem)

correlation_matrix <- cor(score_matrix)

# Self-correlations are exactly 1 and would otherwise dominate the scale; they get their own
# colour (white) by placing them above the top break.
breaks <- seq(min(correlation_matrix),
              max(correlation_matrix[correlation_matrix != 1]),
              length.out = 17)
colors <- c(viridis(length(breaks) - 1), "white")

# pheatmap prints the row/column names verbatim, so swap underscores for spaces for the figure
# only, then restore them so downstream joins keep working.
display_matrix <- correlation_matrix
dimnames(display_matrix) <- lapply(dimnames(display_matrix), gsub, pattern = "_", replacement = " ")

heatmap <- pheatmap(
  display_matrix,
  col        = colors,
  border_color = NA,
  cellwidth  = 12,
  cellheight = 12,
  breaks     = breaks,
  main       = "Average Correlation of S-Value by Condition",
  clustering_distance_cols = "correlation",
  clustering_distance_row  = "correlation",
  clustering_method        = "average",
  filename   = file.path(FIG_DIR, "MedianChemCorByCondition.png")
)

ggsave(
  file.path(FIG_DIR, "MedianChemCorByCondition.svg"),
  plot   = heatmap,
  width  = 0.2 * ncol(score_matrix) + 2,
  height = 0.2 * ncol(score_matrix) + 2,
  limitsize = FALSE,
  device = "svg"
)

# The correlation matrix in the row/column order chosen by the dendrograms.
ordered <- correlation_matrix[heatmap$tree_row$order, heatmap$tree_col$order]
fwrite(data.table(Condition = rownames(ordered), ordered),
       file.path(FIG_DIR, "MedianChemCorByCondition.tsv"), sep = "\t")

# Second-order structure: how similar two conditions' *correlation profiles* are to each other.
second_order <- cor(correlation_matrix)[heatmap$tree_row$order, heatmap$tree_row$order]
write.table(second_order,
            file.path(PROCESSED_DIR, "heatmap_correlation_second_order.tsv"), sep = "\t")

message("Figure 3 written to ", FIG_DIR)
