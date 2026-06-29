# Name: run_subtype_coupling.R
# Author: Andrew Willems <awillems@vols.utk.edu>
# Purpose: Subtype-coupling arm -- does the per-cell EMT axis track SCLC
#   neuroendocrine (NE) state and A/N/P subtype, as expected from SCLC biology
#   (Groves et al. 2023: SCLC-A is NE-high/epithelial-leaning, non-NE leans
#   mesenchymal; EMT should anticorrelate with the NE score)?
#
#   1. Load a single-cell counts matrix (HGNC symbols) + per-cell metadata.
#   2. Score each cell on the EMT axis (score_emt_singlecell, WS1) and the NE
#      axis (ne_score_singlecell, WS2; higher = more neuroendocrine).
#   3. Call the A/N/P subtype per cell (call_sclc_subtype, sparse-aware).
#   4. Cross-tabulate EMT vs NE / subtype (map_emt_to_subtype) and report the
#      EMT~NE correlation (expected NEGATIVE) overall and per group.
#
# Usage: source("R/run_subtype_coupling.R") from the project root, or
#   Rscript R/run_subtype_coupling.R. The big single-cell inputs are not in the
#   repo; produce them with R/fetch_cdx_data.R or R/fetch_chan_data.R.

# =========================================================================
# 0. CONFIGURATION
# =========================================================================

config <- list(
  # Genes-in-rows (HGNC symbols) single-cell counts matrix; .rds of a (sparse) matrix.
  counts_path = "Data/cdx_scrnaseq/cdx_counts.rds",
  # Per-cell metadata TSV with at least a `cell` column (+ any grouping columns
  # to break results down by, e.g. model / condition / group).
  cell_meta_path = "Data/cdx_scrnaseq/cell_metadata.tsv",
  group_col = "model",        # optional per-group breakdown column (NULL to skip)

  sc_method = "UCell",        # per-cell scorer for both EMT and NE
  output_dir = "Outputs/subtype_coupling"
)


# =========================================================================
# 1. SETUP
# =========================================================================

cat("============================================================\n")
cat("  SCLC subtype-coupling: EMT axis vs NE / A-N-P subtype\n")
cat("============================================================\n\n")

source("R/emt_scoring.R")
source("R/emt_subtype_map.R")

dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(config$counts_path) || !file.exists(config$cell_meta_path)) {
  stop(sprintf(paste0(
    "Single-cell inputs not found:\n  %s\n  %s\n",
    "Produce them with R/fetch_cdx_data.R or R/fetch_chan_data.R."),
    config$counts_path, config$cell_meta_path))
}


# =========================================================================
# 2. LOAD
# =========================================================================

cat("Loading single-cell data...\n")
counts <- readRDS(config$counts_path)
cell_meta <- utils::read.delim(config$cell_meta_path, stringsAsFactors = FALSE)
cat(sprintf("  %d genes x %d cells; %d metadata rows\n",
            nrow(counts), ncol(counts), nrow(cell_meta)))


# =========================================================================
# 3. SCORE EMT + NE PER CELL, CALL SUBTYPE
# =========================================================================

sigs <- load_emt_signatures()
cat("\nScoring per-cell EMT axis (higher = mesenchymal)...\n")
emt <- score_emt_singlecell(counts, epithelial_genes = sigs$ks_epithelial,
                            mesenchymal_genes = sigs$ks_mesenchymal,
                            method = config$sc_method)

cat("Scoring per-cell NE axis (higher = neuroendocrine)...\n")
ne <- ne_score_singlecell(counts, ne_template = load_ne_template(),
                          method = config$sc_method)

cat("Calling A/N/P subtype per cell...\n")
subt <- suppressWarnings(call_sclc_subtype(counts))


# =========================================================================
# 4. COUPLE: EMT vs NE / subtype
# =========================================================================

emt_df <- data.frame(sample = names(emt), consensus = as.numeric(emt),
                     stringsAsFactors = FALSE)
ne_df  <- data.frame(sample = names(ne), ne = as.numeric(ne), stringsAsFactors = FALSE)
coupled <- map_emt_to_subtype(emt_df, subt, ne = ne_df)

shared <- Reduce(intersect, list(names(emt), names(ne)))
overall_r <- stats::cor(emt[shared], ne[shared])
overall_rho <- stats::cor(emt[shared], ne[shared], method = "spearman")

cat("\n============================================================\n")
cat("  RESULT: does EMT track NE / subtype?\n")
cat("============================================================\n")
cat(sprintf("\n[overall] cor(EMT, NE) = %.3f (Pearson), %.3f (Spearman)  [expect < 0]\n",
            overall_r, overall_rho))
cat("\n[by subtype] EMT consensus + EMT~NE correlation:\n")
print(as.data.frame(coupled$by_subtype), row.names = FALSE)

# Optional per-group breakdown (e.g. per CDX model) of cor(EMT, NE).
if (!is.null(config$group_col) && config$group_col %in% names(cell_meta)) {
  gmap <- stats::setNames(cell_meta[[config$group_col]], cell_meta$cell)
  g <- gmap[shared]
  cat(sprintf("\n[per %s] cor(EMT, NE):\n", config$group_col))
  for (lev in sort(unique(g[!is.na(g)]))) {
    ix <- which(g == lev)
    if (length(ix) > 10) cat(sprintf("  %-12s n=%6d  r=%+.3f\n", lev, length(ix),
                                      stats::cor(emt[shared][ix], ne[shared][ix])))
  }
}


# =========================================================================
# 5. SAVE
# =========================================================================

per_cell <- merge(merge(emt_df, ne_df, by = "sample"),
                  data.frame(sample = subt$sample, subtype = subt$subtype,
                             stringsAsFactors = FALSE), by = "sample")
names(per_cell)[names(per_cell) == "sample"] <- "cell"
utils::write.csv(per_cell, file.path(config$output_dir, "per_cell_emt_ne_subtype.csv"), row.names = FALSE)
utils::write.csv(as.data.frame(coupled$by_subtype),
                 file.path(config$output_dir, "by_subtype.csv"), row.names = FALSE)

cat(sprintf("\nResults written under: %s/\n", config$output_dir))
cat("============================================================\n")
cat("  Subtype-coupling analysis complete.\n")
cat("============================================================\n")
