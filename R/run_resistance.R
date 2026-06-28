# Name: run_resistance.R
# Author: Andrew Willems <awillems@vols.utk.edu>
# Purpose: WS4 runner — does EMT arise in SCLC as a function of treatment?
#   1. Load the CDX single-cell counts + per-cell metadata (Stewart 2020,
#      GSE138267; see Data/cdx_scrnaseq/PROVENANCE.md).
#   2. Score every cell on the WS1 EMT axis (score_emt_singlecell()).
#   3. Assemble the paired sensitive-vs-resistant frame.
#   4. Test the three WS4 questions and write tidy results + per-model tables:
#        - location   : does the EMT axis shift up with resistance?
#        - dispersion : does EMT heterogeneity widen with resistance?
#        - composition: does the mesenchymal-cell fraction rise with resistance?
#
# Usage: source("R/run_resistance.R") from the project root, or
#   Rscript R/run_resistance.R. The big inputs are NOT in the repo; obtain them
#   per Data/cdx_scrnaseq/PROVENANCE.md and point the config below at them.

# =========================================================================
# 0. CONFIGURATION — edit these paths/params as needed
# =========================================================================

config <- list(
  # --- Inputs (see Data/cdx_scrnaseq/PROVENANCE.md to produce these) ---
  # Genes-in-rows (HGNC symbols) single-cell matrix; .rds of a (sparse) matrix.
  counts_path = "Data/cdx_scrnaseq/cdx_counts.rds",
  # Per-cell metadata TSV with at least: cell, model, condition.
  cell_meta_path = "Data/cdx_scrnaseq/cell_metadata.tsv",

  # If your GEO condition labels are study-specific (e.g. drug names for the
  # resistant arm), map them explicitly here; leave NULL to use keyword parsing.
  sensitive_labels = NULL,
  resistant_labels = NULL,

  # --- Scoring / analysis params ---
  sc_method      = "UCell",   # per-cell EMT scorer (UCell or AUCell)
  min_cells      = 20,        # min cells per model x condition group
  min_models     = 3,         # min paired models for the Wilcoxon test
  dispersion     = "sd",      # heterogeneity metric: sd / var / iqr / mad
  alternative    = "greater", # directional: resistance RAISES EMT

  output_dir = "Outputs/resistance"
)


# =========================================================================
# 1. SETUP
# =========================================================================

cat("============================================================\n")
cat("  SCLC WS4: EMT dynamics across treatment resistance\n")
cat("============================================================\n\n")

source("R/emt_scoring.R")
source("R/emt_resistance.R")

dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(config$counts_path) || !file.exists(config$cell_meta_path)) {
  stop(sprintf(paste0(
    "CDX inputs not found:\n  %s\n  %s\n",
    "Obtain GSE138267 and build them as described in ",
    "Data/cdx_scrnaseq/PROVENANCE.md."),
    config$counts_path, config$cell_meta_path))
}


# =========================================================================
# 2. LOAD DATA
# =========================================================================

cat("Loading CDX single-cell data...\n")
counts <- readRDS(config$counts_path)
cell_meta <- utils::read.delim(config$cell_meta_path, stringsAsFactors = FALSE)
cat(sprintf("  counts: %d genes x %d cells\n", nrow(counts), ncol(counts)))
cat(sprintf("  metadata: %d cells, %d model(s)\n",
            nrow(cell_meta), length(unique(cell_meta$model))))


# =========================================================================
# 3. SCORE THE EMT AXIS PER CELL (WS1)
# =========================================================================

cat("\nScoring per-cell EMT axis (higher = more mesenchymal)...\n")
sigs <- load_emt_signatures()
if (is.null(sigs$ks_epithelial) || is.null(sigs$ks_mesenchymal)) {
  stop("KS epithelial/mesenchymal gene sets missing. Run fetch_emt_signatures() first.")
}
emt_cells <- score_emt_singlecell(
  counts,
  epithelial_genes  = sigs$ks_epithelial,
  mesenchymal_genes = sigs$ks_mesenchymal,
  method = config$sc_method
)
cat(sprintf("  scored %d cells\n", length(emt_cells)))


# =========================================================================
# 4. ASSEMBLE THE PAIRED FRAME
# =========================================================================

cat("\nAssembling paired sensitive-vs-resistant frame...\n")
prepared <- prepare_resistance_emt(
  emt_cells, cell_meta,
  cell_col = "cell",                     # the metadata TSV carries an explicit cell column
  sensitive_labels = config$sensitive_labels,
  resistant_labels = config$resistant_labels,
  min_cells = config$min_cells
)
unpaired <- attr(prepared, "unpaired_models")
cat(sprintf("  %d cells across %d model(s); %d unpaired model(s)%s\n",
            nrow(prepared), length(unique(prepared$model)), length(unpaired),
            if (length(unpaired)) paste0(": ", paste(unpaired, collapse = ", ")) else ""))


# =========================================================================
# 5. THE THREE WS4 TESTS
# =========================================================================

cat("\n--- (a) LOCATION: EMT axis shift -------------------------\n")
loc <- emt_resistance_shift(prepared, alternative = config$alternative,
                            min_models = config$min_models)
print(loc$test)

cat("\n--- (b) DISPERSION: EMT heterogeneity --------------------\n")
het <- emt_heterogeneity_shift(prepared, dispersion = config$dispersion,
                               alternative = config$alternative,
                               min_models = config$min_models)
print(het$test)

cat("\n--- (c) COMPOSITION: mesenchymal-cell fraction -----------\n")
comp <- emt_state_composition(prepared, alternative = config$alternative,
                              min_models = config$min_models)
print(comp$test)


# =========================================================================
# 6. SAVE
# =========================================================================

summary_tbl <- do.call(rbind, list(loc$test, het$test, comp$test))
# Persist the mesenchymal-fraction threshold (an attribute on comp$test) as a
# column so ws4_resistance_summary.csv is unambiguously reproducible.
summary_tbl$mes_threshold <- NA_real_
summary_tbl$mes_threshold[summary_tbl$quantity == "mes_frac"] <- attr(comp$test, "threshold")
utils::write.csv(summary_tbl, file.path(config$output_dir, "ws4_resistance_summary.csv"),
                 row.names = FALSE)
utils::write.csv(loc$per_model,  file.path(config$output_dir, "per_model_location.csv"),  row.names = FALSE)
utils::write.csv(het$per_model,  file.path(config$output_dir, "per_model_dispersion.csv"), row.names = FALSE)
utils::write.csv(comp$per_model, file.path(config$output_dir, "per_model_composition.csv"), row.names = FALSE)
utils::write.csv(as.data.frame(prepared),
                 file.path(config$output_dir, "per_cell_emt.csv"), row.names = FALSE)

cat(sprintf("\nResults written under: %s/\n", config$output_dir))
cat("============================================================\n")
cat("  WS4 resistance analysis complete.\n")
cat("============================================================\n")
