# ============================================================================
# fetch_chan_data.R
#
# WS4 external validation loader: read the Chan et al. 2021 human SCLC atlas
# (Cancer Cell; CZ CELLxGENE collection 62e8f058-...) from its CELLxGENE .h5ad
# and produce the per-cell inputs the WS4 cross-sectional test consumes -- a
# genes-in-rows counts matrix (HGNC symbols) of the malignant SCLC cells plus a
# per-cell metadata table with (sample, donor, treatment group). See
# Data/chan_sclc/PROVENANCE.md for the source and citation.
#
# Why this cohort: independent HUMAN primary SCLC tumors (vs the Stewart CDX
# models), scRNA-seq, with a treatment-naive vs platinum-treated split -- a
# genuine cross-platform, between-patient replication of the EMT-heterogeneity
# finding. The design is CROSS-SECTIONAL (each tumor is wholly naive or treated),
# so the test is emt_dispersion_groupwise(), not the paired CDX tests.
#
# The .h5ad is read with {rhdf5} (no Python/basilisk). CELLxGENE schema details
# used here (verified against the file):
#   - raw/X : CSR (csr_matrix) integer counts, shape [n_cells, n_genes].
#   - var/feature_name : HGNC symbols (var/_index is Ensembl).
#   - obs/disease, obs/cell_type_general, obs/treatment, obs/donor_id,
#     obs/HTAN_Biospecimen_ID, obs/libsize : the columns we filter/group on.
#
# Conventions match the other modules: {cli} errors, tidy output, soft deps.
# ============================================================================

#' tibble if available else data.frame
#' @keywords internal
.chan_as_tbl <- function(df) {
  if (requireNamespace("tibble", quietly = TRUE)) tibble::as_tibble(df) else df
}

#' @keywords internal
.chan_require <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    how <- if (bioc) sprintf("BiocManager::install('%s')", pkg) else sprintf("install.packages('%s')", pkg)
    cli::cli_abort(c("x" = "Package {.pkg {pkg}} is required.", "i" = "Install with {.code {how}}."))
  }
}

#' Read one obs/var column from an .h5ad, decoding anndata categoricals
#'
#' Categorical columns are stored as a group with `categories` + `codes`; plain
#' columns are a dataset. Returns a character/numeric vector in cell order.
#' @keywords internal
.chan_read_col <- function(file, path) {
  x <- rhdf5::h5read(file, path)
  if (is.list(x) && all(c("categories", "codes") %in% names(x))) {
    cats <- as.character(x$categories)
    codes <- as.integer(x$codes)
    codes[codes < 0] <- NA_integer_           # anndata uses -1 for NA
    return(cats[codes + 1L])
  }
  as.vector(x)
}

#' Map the Chan `treatment` field to a two-level naive/treated factor
#'
#' Everything equal to `naive_label` is "naive"; every other (non-NA) value --
#' all of which are platinum-doublet-containing regimens in the SCLC arm -- is
#' "treated".
#'
#' @param treatment Character vector of the obs `treatment` values.
#' @param naive_label The value denoting treatment-naive. Default "Naive".
#' @return Factor with levels c("naive", "treated"); NA stays NA.
#' @export
chan_treatment_group <- function(treatment, naive_label = "Naive") {
  t <- as.character(treatment)
  out <- rep(NA_character_, length(t))
  out[!is.na(t) & t == naive_label] <- "naive"
  out[!is.na(t) & t != naive_label] <- "treated"
  factor(out, levels = c("naive", "treated"))
}

#' Read a CELLxGENE .h5ad count matrix as a symbols x cells sparse matrix
#'
#' Reads the CSR count layer and builds a genes(symbols)-in-rows `dgCMatrix`,
#' optionally restricted to a subset of cells (1-based positions in obs order).
#' Duplicate gene symbols are collapsed by summing their rows.
#'
#' @param file Path to the .h5ad.
#' @param keep_cells Integer positions of cells to keep (default: all).
#' @param layer "raw" (default; raw/X counts) or "X" (the main, normalized X).
#' @param symbol_path obs-relative path to gene symbols. Default "var/feature_name".
#' @return A sparse `dgCMatrix`, genes (symbols) in rows, kept cells in columns
#'   (column names from `obs/_index`).
#' @export
read_h5ad_counts <- function(file, keep_cells = NULL, layer = c("raw", "X"),
                             symbol_path = "var/feature_name") {
  .chan_require("rhdf5", bioc = TRUE)
  .chan_require("Matrix", bioc = FALSE)
  grp <- if (match.arg(layer) == "raw") "raw/X" else "X"
  attrs <- rhdf5::h5readAttributes(file, grp)
  if (!identical(as.character(attrs[["encoding-type"]]), "csr_matrix")) {
    cli::cli_abort("Expected a CSR matrix at {.val {grp}}; got {.val {attrs[['encoding-type']]}}.")
  }
  shape <- as.integer(attrs[["shape"]])          # [n_cells, n_genes]
  ncells <- shape[1]; ngenes <- shape[2]
  if (is.null(keep_cells)) keep_cells <- seq_len(ncells)
  keep_cells <- as.integer(keep_cells)
  if (any(keep_cells < 1 | keep_cells > ncells)) cli::cli_abort("{.arg keep_cells} out of range.")

  # as.double/as.integer also strip the 1-D array attribute h5read() returns,
  # which sparseMatrix() rejects for `x`.
  indptr <- as.double(rhdf5::h5read(file, paste0(grp, "/indptr")))  # 0-based, len ncells+1
  data <- as.double(rhdf5::h5read(file, paste0(grp, "/data")))      # nonzero values
  indices <- as.integer(rhdf5::h5read(file, paste0(grp, "/indices")))  # 0-based gene idx
  # Gather each kept cell's contiguous nonzero block (CSR row) -> CSC columns.
  starts <- indptr[keep_cells]; ends <- indptr[keep_cells + 1L]
  lens <- as.integer(ends - starts)
  sel <- sequence(lens, from = starts + 1)        # 1-based positions into data/indices
  m <- Matrix::sparseMatrix(
    i = indices[sel] + 1L, p = c(0L, cumsum(lens)), x = data[sel],
    dims = c(ngenes, length(keep_cells)), index1 = TRUE
  )
  rm(data, indices, sel)
  symbols <- as.character(.chan_read_col(file, symbol_path))
  cellids <- as.character(.chan_read_col(file, "obs/_index"))[keep_cells]
  dimnames(m) <- list(symbols, cellids)
  if (anyDuplicated(symbols)) {                   # collapse duplicate symbols by summing
    g <- factor(symbols, levels = unique(symbols))
    m <- methods::as(Matrix::fac2sparse(g), "CsparseMatrix") %*% m
    rownames(m) <- levels(g)
  }
  m
}

#' Load the malignant SCLC cells from the Chan atlas for WS4 validation
#'
#' Filters to SCLC tumor (epithelial) cells, reads their raw counts, builds the
#' per-cell metadata (sample = biospecimen, donor, naive/treated group), and
#' applies a minimum-UMI QC consistent with the CDX analysis.
#'
#' @param file Path to the Chan combined .h5ad.
#' @param min_umi Minimum library size (`obs/libsize`) per cell. Default 1000.
#' @param disease_value obs `disease` value selecting SCLC. Default
#'   "small cell lung carcinoma".
#' @param epithelial_value obs `cell_type_general` value selecting tumor cells.
#'   Default "Epithelial".
#' @param sample_col obs column used as the sample (replication) unit. Default
#'   "HTAN_Biospecimen_ID".
#' @return A list: `counts` (symbols x cells) and `cell_meta` (data.frame:
#'   `cell`, `sample`, `donor`, `group`, `libsize`).
#' @export
load_chan_sclc <- function(file, min_umi = 1000,
                           disease_value = "small cell lung carcinoma",
                           epithelial_value = "Epithelial",
                           sample_col = "HTAN_Biospecimen_ID") {
  .chan_require("rhdf5", bioc = TRUE)
  .chan_require("Matrix", bioc = FALSE)
  if (!file.exists(file)) cli::cli_abort("h5ad not found: {.path {file}}.")
  disease <- .chan_read_col(file, "obs/disease")
  ctype   <- .chan_read_col(file, "obs/cell_type_general")
  treat   <- .chan_read_col(file, "obs/treatment")
  donor   <- .chan_read_col(file, "obs/donor_id")
  sampid  <- .chan_read_col(file, paste0("obs/", sample_col))
  group   <- chan_treatment_group(treat)

  # Select SCLC malignant cells first, then QC on the actual library size
  # (obs/libsize is log10-transformed, so we recompute from raw counts to match
  # the CDX pipeline's min-UMI threshold exactly).
  keep0 <- which(disease == disease_value & ctype == epithelial_value & !is.na(group))
  if (length(keep0) == 0) cli::cli_abort("No SCLC tumor cells match the disease/cell-type filters.")
  counts0 <- read_h5ad_counts(file, keep_cells = keep0, layer = "raw")
  depth <- Matrix::colSums(counts0)
  pass <- depth >= min_umi
  if (sum(pass) == 0) cli::cli_abort("No SCLC tumor cells reach {min_umi} UMI.")
  if (sum(!pass) > 0) cli::cli_inform("UMI QC: dropping {sum(!pass)}/{length(pass)} SCLC cells with < {min_umi} counts.")
  counts <- counts0[, pass, drop = FALSE]
  keep <- keep0[pass]
  cli::cli_inform(c("v" = "Chan SCLC: {ncol(counts)} tumor cells across {length(unique(sampid[keep]))} samples ({length(unique(donor[keep]))} donors)."))

  cell_meta <- data.frame(
    cell = colnames(counts), sample = as.character(sampid[keep]),
    donor = as.character(donor[keep]), group = as.character(group[keep]),
    libsize = as.numeric(depth[pass]), row.names = NULL, stringsAsFactors = FALSE
  )
  list(counts = counts, cell_meta = cell_meta)
}
