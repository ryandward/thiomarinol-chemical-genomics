# Shared setup: paths, packages, and the median S-score matrix used by both figures.
# Source this from the repository root, or let the figure scripts source it themselves.

library(data.table)

# All paths are relative to the repository root.
if (!file.exists("data/raw/Bsu_CRISPRi_S-scores_w_Bo_chemicals.tsv")) {
  stop("Run from the repository root: setwd() to the directory containing this repo's README.md")
}

RAW_SCORES  <- "data/raw/Bsu_CRISPRi_S-scores_w_Bo_chemicals.tsv"
CONDITIONS  <- "data/raw/conditions.tsv"
BSU_TO_GENE <- "data/raw/BSUtoGene.tsv"
PEPTIDOGLYCAN <- "data/raw/annotations/peptidoglycan.tsv"
CELL_SHAPE    <- "data/raw/annotations/cell_shape.tsv"
FIG_DIR       <- "figures"
PROCESSED_DIR <- "data/processed"

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PROCESSED_DIR, showWarnings = FALSE, recursive = TRUE)

# Conditions excluded from the published analysis: not chemical treatments.
NON_DRUG_CONDITIONS <- c("UV", "Minimal_S7-50")

# ftsA was excluded from the published analysis.
EXCLUDED_GENES <- c("ftsA")

# Rows of the raw table are one chemical at one concentration, labelled "Mupirocin-C [C]", where the
# bracketed letter indexes the dose. The published Data Set 1 writes the same row as
# "mupirocin [C] 0.03", with the concentration in ug/mL. Strip the tag to get the chemical.
condition_of <- function(x) sub("-[A-Za-z0-9]+ \\[[A-Za-z0-9]+\\]$", "", x)

# Genes x conditions matrix of S-scores, taking the median over each chemical's concentration
# series. Conditions assayed at a single concentration pass through unchanged.
load_median_s_scores <- function() {
  chem_raw   <- fread(RAW_SCORES, sep = "\t")
  conditions <- fread(CONDITIONS, sep = "\t")

  chemical <- condition_of(chem_raw$Condition)
  unknown <- setdiff(unique(chemical), conditions$Condition)
  if (length(unknown)) stop("Row labels not present in conditions.tsv: ", paste(unknown, collapse = ", "))

  median_chem <- melt(chem_raw, id.vars = "Condition", variable.name = "Gene",
                      value.name = "S-Value")[, .(Gene = unique(Gene))]

  for (condition in conditions$Condition) {
    median_chem[[condition]] <- melt(
      chem_raw[chemical == condition],
      id.vars = "Condition", variable.name = "Gene", value.name = "S-Value"
    )[, .(`S-Value` = median(`S-Value`)), by = .(Gene)]$`S-Value`
  }

  median_chem <- median_chem[, setdiff(names(median_chem), NON_DRUG_CONDITIONS), with = FALSE]
  median_chem <- median_chem[!Gene %in% EXCLUDED_GENES]
  na.omit(median_chem)
}

# Same table as a numeric matrix with gene names on the rows.
as_score_matrix <- function(median_chem) {
  m <- data.matrix(median_chem[, -1])
  rownames(m) <- median_chem$Gene
  m
}
