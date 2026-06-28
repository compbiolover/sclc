# ============================================================================
# fetch_cdx_data.R
#
# WS4 data loader: turn the Stewart CDX single-cell series (GEO GSE138267) into
# the two inputs R/emt_resistance.R consumes -- a genes-in-rows counts matrix
# (HGNC symbols) and a per-cell metadata table with (model, condition). See
# Data/cdx_scrnaseq/PROVENANCE.md for the source and citation.
#
# Layout of the series (verified against GEO):
#   - 25 GSM samples, each a CDX library titled "<MODEL><variant>.<LIBRARY>"
#     (e.g. "SC4.LB17009", "SC4_Talazoparib.LB17004", "SC68_cr.LB19003",
#     "SC53cis", "SC55-3.LB17026"), with a `treatment` characteristic.
#   - Each GSM ships a 10x triplet (barcodes.tsv / genes.tsv / matrix.mtx) inside
#     "<GSM>_<title>.matrix.tar.gz". genes.tsv is `<Ensembl>\t<symbol>` on the
#     hg19 reference (human only); matrix.mtx is genes x cells integer counts.
#
# THE CRITICAL JOIN (get this wrong and the analysis is meaningless): the
# biological CDX model is the `SC<number>` stem, NOT the dotted prefix. The drug
# tags (`_LY2606368`, `_Talazoparib`, `_cr`, `cis`, `-1/-2/-3`) are
# treatment-derived VARIANTS of the same model. So SC4, SC4_LY2606368 and
# SC4_Talazoparib are one model (SC4) in three states; keying on the dotted
# prefix would shatter it into unpaired singletons and destroy the pairing.
#   -> base model  = regmatches(title, "^SC[0-9]+")
#   -> condition   : "untreated"/"vehicle-treated" = sensitive;
#                    anything "...relapsed"          = resistant;
#                    bare "cisplatin" (on-treatment, not yet relapsed) = NA,
#                    excluded by default as neither a naive baseline nor an
#                    established resistant line.
# This yields 4 paired models (SC4, SC53, SC55, SC68) and 4 sensitive-only
# models (SC16, SC39, SC49, SC75).
#
# The parsing (pure, testable offline) is separated from the network fetch.
# GEOquery / Matrix are soft deps used only by the fetching/reading functions.
# Conventions match the other modules: {cli} errors, tidy output.
# ============================================================================

#' tibble if available else data.frame
#' @keywords internal
.cdx_as_tbl <- function(df) {
  if (requireNamespace("tibble", quietly = TRUE)) tibble::as_tibble(df) else df
}

#' Require a soft dependency with an actionable error
#' @keywords internal
.cdx_require <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    how <- if (bioc) sprintf("BiocManager::install('%s')", pkg) else sprintf("install.packages('%s')", pkg)
    cli::cli_abort(c("x" = "Package {.pkg {pkg}} is required for this step.",
                     "i" = "Install it with {.code {how}}."))
  }
}

# ============================================================================
# 1. Pure parsing: title -> base model, treatment -> condition
# ============================================================================

#' Extract the base CDX model (the SC-number stem) from sample titles
#'
#' @param title Character vector of GEO sample titles.
#' @return Character vector of base models (e.g. "SC4", "SC55"); errors if any
#'   title does not start with an `SC<number>` stem.
#' @export
cdx_base_model <- function(title) {
  # Check the regexpr() match position directly (-1 / NA = no match) so a bad
  # title can't slip through; regmatches() silently DROPS non-matches, which
  # would otherwise misalign the result with `title`.
  pos <- regexpr("^SC[0-9]+", title)
  bad <- is.na(pos) | pos < 0
  if (any(bad)) {
    cli::cli_abort(c("x" = "{sum(bad)} title(s) lack an SC<number> stem: {.val {title[bad]}}.",
                     "i" = "Expected titles like {.val SC4.LB17009} or {.val SC53cis}."))
  }
  regmatches(title, pos)
}

#' Map a CDX `treatment` characteristic to a sensitive/resistant condition
#'
#' @param treatment Character vector of the GEO `treatment` characteristic.
#' @return Factor (levels sensitive < resistant). On-treatment-but-not-relapsed
#'   values (e.g. bare "cisplatin") and anything unrecognized map to `NA` --
#'   they are neither a naive baseline nor an established resistant line.
#' @export
cdx_condition <- function(treatment) {
  t <- tolower(trimws(as.character(treatment)))
  out <- rep(NA_character_, length(t))
  out[grepl("relaps", t)] <- "resistant"                       # "...relapsed"
  out[is.na(out) & grepl("untreated|vehicle", t)] <- "sensitive"
  factor(out, levels = c("sensitive", "resistant"))
}

#' Classify CDX samples into model + condition (pure; no network)
#'
#' The offline core of [cdx_sample_table()]: given a data.frame with sample
#' titles and treatments, derive the base model, the sensitive/resistant
#' condition, and per-model pairing.
#'
#' @param df A data.frame with a title column and a treatment column.
#' @param title_col,treatment_col Column names. Defaults "title","treatment".
#' @return The input as a tibble with added columns: `model`, `condition`
#'   (factor, NA where excluded), `excluded` (TRUE for NA-condition rows), and
#'   `paired` (TRUE if the row's model has >=1 sensitive AND >=1 resistant
#'   sample among non-excluded rows).
#' @export
cdx_classify_samples <- function(df, title_col = "title", treatment_col = "treatment") {
  miss <- setdiff(c(title_col, treatment_col), names(df))
  if (length(miss) > 0) cli::cli_abort("Missing column(s): {.val {miss}}.")
  out <- as.data.frame(df, stringsAsFactors = FALSE)
  out$model <- cdx_base_model(as.character(out[[title_col]]))
  out$condition <- cdx_condition(out[[treatment_col]])
  out$excluded <- is.na(out$condition)
  # A model is paired iff it has both states among its non-excluded samples.
  has_both <- tapply(as.character(out$condition), out$model, function(z) {
    z <- z[!is.na(z)]
    all(c("sensitive", "resistant") %in% z)
  })
  out$paired <- unname(has_both[out$model])
  .cdx_as_tbl(out)
}

# ============================================================================
# 2. Network: build the GEO sample table
# ============================================================================

#' Build the CDX sample table from GEO (needs network + GEOquery)
#'
#' @param gse GEO series accession. Default "GSE138267" (the scRNA-seq sub-series
#'   of the GSE138474 SuperSeries).
#' @param destdir Directory to cache the downloaded series matrix.
#'   Default `tempdir()`.
#' @return A tibble: `gsm`, `title`, `treatment`, `model`, `condition`,
#'   `excluded`, `paired`, `suppl_url` (the per-sample matrix.tar.gz). Emits a
#'   summary of paired/unpaired/excluded samples.
#' @export
cdx_sample_table <- function(gse = "GSE138267", destdir = tempdir()) {
  .cdx_require("GEOquery", bioc = TRUE)
  es <- GEOquery::getGEO(gse, GSEMatrix = TRUE, getGPL = FALSE, destdir = destdir)[[1]]
  pd <- Biobase::pData(es)
  treat_col <- if ("treatment:ch1" %in% names(pd)) "treatment:ch1" else
    grep("^treatment", names(pd), value = TRUE)[1]
  if (is.na(treat_col)) cli::cli_abort("No treatment characteristic found in {gse} pData.")
  base <- data.frame(
    gsm = pd$geo_accession,
    title = pd$title,
    treatment = pd[[treat_col]],
    suppl_url = as.character(pd$supplementary_file_1),
    stringsAsFactors = FALSE
  )
  tab <- cdx_classify_samples(base)
  n_pair_models <- length(unique(tab$model[tab$paired]))
  cli::cli_inform(c(
    "v" = "{gse}: {nrow(tab)} samples, {length(unique(tab$model))} models.",
    "i" = "Paired models (sensitive + resistant): {n_pair_models} -- {.val {sort(unique(tab$model[tab$paired]))}}.",
    "i" = "Excluded (on-treatment/ambiguous) samples: {sum(tab$excluded)}."
  ))
  tab
}

# ============================================================================
# 3. Read a 10x triplet -> symbols x cells sparse matrix
# ============================================================================

#' Read a 10x-style genes/barcodes/matrix bundle into a symbols x cells matrix
#'
#' Handles BOTH layouts present in GSE138267: the LB17 batch ships a sparse
#' MatrixMarket triplet (`genes.tsv`/`barcodes.tsv`/`matrix.mtx`), the LB19 batch
#' ships a DENSE tab-delimited table (`genes.txt`/`barcodes.txt`/`matrix.txt`,
#' genes x cells). The matrix file is detected by its first line.
#'
#' @param dir Directory containing the bundle (files optionally gzipped).
#' @param cell_prefix String prepended to each barcode (with "|") to keep cell
#'   ids unique across samples. Default "" (no prefix).
#' @param symbol_col Column of the genes file holding HGNC symbols. Default 2.
#' @return A sparse `dgCMatrix`, genes (symbols) in rows, cells in columns.
#'   Duplicate symbols (the 10x reference has some) are collapsed by summing
#'   their counts.
#' @export
read_10x_triplet <- function(dir, cell_prefix = "", symbol_col = 2) {
  .cdx_require("Matrix", bioc = FALSE)
  # Match by SUFFIX so both bare ("matrix.mtx") and sample-prefixed
  # ("SC53cis.matrix.mtx") names are found across the series' batches.
  pick <- function(suffixes) {
    for (s in suffixes) {
      hits <- list.files(dir, pattern = paste0(s, "(\\.gz)?$"), full.names = TRUE)
      if (length(hits) > 0) return(hits[1])
    }
    cli::cli_abort("No file matching {.val {suffixes}} found in {.path {dir}}.")
  }
  mfile <- pick(c("matrix\\.mtx", "matrix\\.txt"))
  genes <- utils::read.delim(pick(c("genes\\.tsv", "genes\\.txt", "features\\.tsv")),
                             header = FALSE, stringsAsFactors = FALSE)
  barcodes <- utils::read.delim(pick(c("barcodes\\.tsv", "barcodes\\.txt")),
                                header = FALSE, stringsAsFactors = FALSE)[[1]]
  ng <- nrow(genes); nc <- length(barcodes)
  # Open a fresh (possibly gzip) connection each time -- readMM()/scan() do not
  # transparently decompress a ".gz" *path*, so a connection is required.
  gz <- grepl("\\.gz$", mfile)
  mcon <- function() if (gz) gzfile(mfile) else mfile
  # Detect MatrixMarket vs dense table from the first non-comment line.
  con <- if (gz) gzfile(mfile) else file(mfile)
  first <- readLines(con, n = 1); close(con)
  if (grepl("^%%MatrixMarket", first)) {
    m <- methods::as(Matrix::readMM(mcon()), "CsparseMatrix")
  } else {
    # Dense genes x cells table: scan row-major (fast, C-level), then reshape.
    vals <- scan(mcon(), what = integer(), quiet = TRUE, sep = "\t")
    if (length(vals) != ng * nc) {
      cli::cli_abort(c("x" = "Dense matrix has {length(vals)} values, expected {ng}x{nc}={ng*nc}.",
                       "i" = "Check {.path {mfile}} against the genes/barcodes files."))
    }
    m <- methods::as(Matrix::Matrix(matrix(vals, nrow = ng, ncol = nc, byrow = TRUE),
                                    sparse = TRUE), "CsparseMatrix")
    rm(vals)
  }
  if (ng != nrow(m) || nc != ncol(m)) {
    cli::cli_abort(c("x" = "Dimensions disagree in {.path {dir}}.",
                     "i" = "matrix {nrow(m)}x{ncol(m)}, genes {ng}, barcodes {nc}."))
  }
  symbols <- as.character(genes[[symbol_col]])
  cells <- if (nzchar(cell_prefix)) paste0(cell_prefix, "|", barcodes) else barcodes
  dimnames(m) <- list(symbols, cells)
  # Collapse duplicate symbols by summing their rows (sparse-friendly).
  if (anyDuplicated(symbols)) {
    grp <- factor(symbols, levels = unique(symbols))
    agg <- methods::as(Matrix::fac2sparse(grp), "CsparseMatrix") %*% m
    rownames(agg) <- levels(grp)
    m <- agg
  }
  m
}

# ============================================================================
# 4. Download + assemble the combined matrix and per-cell metadata
# ============================================================================

#' Download (cached) and read one GSM's matrix tarball
#' @keywords internal
.cdx_fetch_one <- function(gsm, suppl_url, destdir, prefix) {
  .cdx_require("GEOquery", bioc = TRUE)
  gdir <- file.path(destdir, gsm)
  tar <- file.path(gdir, basename(sub("\\?.*$", "", suppl_url)))
  if (!file.exists(tar)) {
    GEOquery::getGEOSuppFiles(gsm, baseDir = destdir, fetch_files = TRUE)
    found <- list.files(gdir, pattern = "matrix\\.tar\\.gz$", full.names = TRUE)
    if (length(found) == 0) cli::cli_abort("No matrix tarball downloaded for {gsm}.")
    tar <- found[1]
  }
  # Clear any prior extraction so a recursive search can't pick up a stale
  # matrix file (e.g. if the tarball was refreshed) -- keeps caching deterministic.
  exdir <- file.path(gdir, "triplet")
  unlink(exdir, recursive = TRUE)
  dir.create(exdir, showWarnings = FALSE, recursive = TRUE)
  utils::untar(tar, exdir = exdir)
  # Some tarballs nest the bundle in a subfolder; find where the matrix landed
  # (LB17 ships matrix.mtx, LB19 ships matrix.txt).
  mtx <- list.files(exdir, pattern = "matrix\\.(mtx|txt)(\\.gz)?$",
                    recursive = TRUE, full.names = TRUE)
  if (length(mtx) == 0) cli::cli_abort("No matrix.(mtx|txt) found after untar of {gsm}.")
  read_10x_triplet(dirname(mtx[1]), cell_prefix = prefix)
}

#' Load CDX single-cell counts + per-cell metadata for the analysis
#'
#' Downloads each selected sample's matrix (cached under `destdir`), reads the
#' 10x triplet, prefixes barcodes by library title so cell ids stay unique,
#' aligns all samples on shared gene symbols, and column-binds them into one
#' sparse matrix. By default only samples with a resolved sensitive/resistant
#' condition are loaded (on-treatment/ambiguous samples are skipped).
#'
#' @param sample_table Output of [cdx_sample_table()] (or the offline
#'   [cdx_classify_samples()]), needing `gsm`,`title`,`model`,`condition`,
#'   `suppl_url`.
#' @param destdir Cache/download directory. Default "Data/cdx_scrnaseq/download".
#' @param models Optional character vector to restrict to specific base models
#'   (e.g. only the paired ones). Default NULL = all eligible.
#' @param include_excluded If TRUE, also load NA-condition samples. Default FALSE.
#' @return A list: `counts` (sparse symbols x cells) and `cell_meta` (data.frame:
#'   `cell`, `gsm`, `library`, `model`, `condition`).
#' @export
load_cdx_counts <- function(sample_table, destdir = "Data/cdx_scrnaseq/download",
                            models = NULL, include_excluded = FALSE) {
  need <- c("gsm", "title", "model", "condition", "suppl_url")
  miss <- setdiff(need, names(sample_table))
  if (length(miss) > 0) cli::cli_abort("{.arg sample_table} missing column(s): {.val {miss}}.")
  st <- as.data.frame(sample_table, stringsAsFactors = FALSE)
  if (!include_excluded) st <- st[!is.na(st$condition), , drop = FALSE]
  if (!is.null(models)) st <- st[st$model %in% models, , drop = FALSE]
  if (nrow(st) == 0) cli::cli_abort("No samples to load after filtering.")
  dir.create(destdir, showWarnings = FALSE, recursive = TRUE)

  mats <- vector("list", nrow(st))
  for (i in seq_len(nrow(st))) {
    cli::cli_inform("Loading {st$title[i]} ({st$gsm[i]}, {i}/{nrow(st)})...")
    mats[[i]] <- .cdx_fetch_one(st$gsm[i], st$suppl_url[i], destdir, prefix = st$title[i])
  }
  # Align on shared gene symbols (the references match, so this is a no-op in
  # practice, but guards against a sample built on a different reference).
  common <- Reduce(intersect, lapply(mats, rownames))
  if (length(common) < 1000) {
    cli::cli_abort("Only {length(common)} gene symbols shared across samples; check references.")
  }
  mats <- lapply(mats, function(m) m[common, , drop = FALSE])
  counts <- do.call(cbind, mats)

  cell_meta <- do.call(rbind, lapply(seq_len(nrow(st)), function(i) {
    data.frame(
      cell = colnames(mats[[i]]),
      gsm = st$gsm[i], library = st$title[i],
      model = st$model[i], condition = as.character(st$condition[i]),
      stringsAsFactors = FALSE
    )
  }))
  cli::cli_inform(c("v" = "Loaded {ncol(counts)} cells x {nrow(counts)} genes across {nrow(st)} samples."))
  list(counts = counts, cell_meta = cell_meta)
}

#' Write the WS4 inputs that run_resistance.R expects
#'
#' @param cdx Output of [load_cdx_counts()].
#' @param dir Destination. Default "Data/cdx_scrnaseq".
#' @return (Invisibly) the paths written.
#' @export
write_cdx_inputs <- function(cdx, dir = "Data/cdx_scrnaseq") {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  counts_path <- file.path(dir, "cdx_counts.rds")
  meta_path <- file.path(dir, "cell_metadata.tsv")
  saveRDS(cdx$counts, counts_path)
  utils::write.table(cdx$cell_meta, meta_path, sep = "\t", row.names = FALSE, quote = FALSE)
  cli::cli_inform(c("v" = "Wrote {.path {counts_path}} and {.path {meta_path}}."))
  invisible(c(counts = counts_path, meta = meta_path))
}

#' One-call convenience: GEO -> written WS4 inputs
#'
#' @param gse,destdir,dir,models Passed through to the steps above.
#' @return (Invisibly) the written paths.
#' @export
fetch_cdx_data <- function(gse = "GSE138267",
                           destdir = "Data/cdx_scrnaseq/download",
                           dir = "Data/cdx_scrnaseq", models = NULL) {
  tab <- cdx_sample_table(gse, destdir = destdir)
  cdx <- load_cdx_counts(tab, destdir = destdir, models = models)
  write_cdx_inputs(cdx, dir = dir)
}
