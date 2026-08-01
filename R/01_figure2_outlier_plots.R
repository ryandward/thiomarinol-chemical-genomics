# Figure 2: median S-score per gene, ordered by position in the B. subtilis genome, with the
# top and bottom 2.5% of genes labelled. Panel A is Mupirocin, panel B is PAC-holo.
# Also produces the Holomycin, Thiolutin and Gliotoxin panels used in the supplement.

source("R/00_common.R")

library(ggplot2)
library(ggrepel)
library(svglite)

peptidoglycan_process <- fread(PEPTIDOGLYCAN)
shape_process         <- fread(CELL_SHAPE)
bsu_to_gene           <- fread(BSU_TO_GENE)

median_chem     <- load_median_s_scores()
median_bsu_chem <- bsu_to_gene[median_chem, on = .(Gene)]
setorder(median_bsu_chem, +BSU)

PANEL_CONDITIONS <- c("Holomycin", "Thiolutin", "Gliotoxin", "PAC-holo", "Mupirocin")
CUTOFF <- 0.025

# The 2021 ggsave() calls passed no size, so the panels inherited whatever the interactive RStudio
# plot pane happened to be: 623.25 x 410.25 pt. Pinned here so a headless run matches the published
# canvas instead of falling back to the 7x7 inch default.
PANEL_WIDTH_IN  <- 623.25 / 72
PANEL_HEIGHT_IN <- 410.25 / 72

# Shared y-axis across the five panels. Wide enough for all of them (Holomycin is the widest at
# -7.03 to 4.82), but not for every condition in the screen: the full matrix reaches -9.05. Adding a
# condition to PANEL_CONDITIONS without widening this would silently drop points, so it is checked.
PANEL_YLIM <- c(-7.1, 4.9)

plot_outliers <- function(this_condition, cutoff_value) {
  d <- data.table(
    `S-Value`        = median_bsu_chem[[this_condition]],
    `BSU Gene Order` = median_bsu_chem[, .I],
    Gene             = median_bsu_chem$Gene,
    BSU              = median_bsu_chem$BSU
  )

  # Guide number and gene name are packed into one label, e.g. "murB2"; split them so the gene
  # name can be italicised while the guide number stays upright.
  d[, guide_number := gsub("[^0-9.]", "", Gene)]
  d[, gene_name    := gsub("[^A-z.]", "", Gene)]
  d[grep("[0-9]", guide_number), stylized_gene := paste0("italic('", gene_name, "')~", guide_number)]
  d[guide_number == "",          stylized_gene := paste0("italic('", gene_name, "')")]

  lower <- quantile(d$`S-Value`, cutoff_value)
  upper <- quantile(d$`S-Value`, 1 - cutoff_value)
  is_outlier <- d$`S-Value` > upper | d$`S-Value` < lower

  d[, Key := paste0("Inner ", 100 - cutoff_value * 100 * 2, "%")]
  d[is_outlier, Key := paste0("Outlier (Top or Bottom ", cutoff_value * 100, "%)")]

  if (this_condition == "PAC-holo") {
    d[is_outlier & Gene %in% peptidoglycan_process$Gene, Key := "Peptidoglycan Outlier"]
    d[is_outlier & Gene %in% shape_process$Gene,         Key := "Shape Outlier"]
    d[is_outlier & Gene %in% shape_process$Gene &
      Gene %in% peptidoglycan_process$Gene,              Key := "Peptidoglycan and Shape Outlier"]
  }

  d[Gene == "ileS", Key := "ileS"]

  p <- ggplot(d, aes(x = `BSU Gene Order`, y = `S-Value`)) +
    geom_point(aes(color = Key), size = 2) +
    scale_color_manual(values = c("#D81B60", "#9e9e9e", "#212121", "#1E88E5", "#FFC107")) +
    theme_bw(base_size = 12) +
    geom_label_repel(
      data = d[is_outlier | Key == "ileS"],
      aes(label = stylized_gene),
      size = 5, box.padding = unit(0.5, "lines"), point.padding = unit(0.5, "lines"),
      max.iter = 500000, parse = TRUE, seed = 1
    ) +
    ggtitle(this_condition) +
    theme(
      plot.title   = element_text(hjust = 0.5, size = 20),
      axis.text    = element_text(size = 14, color = "black"),
      axis.title   = element_text(size = 14, color = "black"),
      legend.text  = element_text(size = 8,  color = "black"),
      legend.title = element_text(size = 14, color = "black"),
      legend.position = "bottom"
    ) +
    ylim(PANEL_YLIM)

  clipped <- sum(d$`S-Value` < PANEL_YLIM[1] | d$`S-Value` > PANEL_YLIM[2])
  if (clipped > 0) {
    warning(sprintf("%s: %d point(s) fall outside PANEL_YLIM and will not be drawn (range %.2f to %.2f)",
                    this_condition, clipped, min(d$`S-Value`), max(d$`S-Value`)), call. = FALSE)
  }

  ggsave(file.path(FIG_DIR, paste0(this_condition, ".svg")), plot = p,
         width = PANEL_WIDTH_IN, height = PANEL_HEIGHT_IN)

  # Per-condition median S-score, one row per guide. The equivalent files shipped in 2021 were
  # collapsed to one row per locus tag, which drops guides where several target the same locus.
  fwrite(d[, .(BSU, Gene, `S-Value`)],
         file.path(FIG_DIR, paste0(this_condition, "_median_by_guide.tsv")), sep = "\t")

  invisible(p)
}

for (condition in PANEL_CONDITIONS) {
  message("Figure 2 panel: ", condition)
  plot_outliers(condition, CUTOFF)
}
