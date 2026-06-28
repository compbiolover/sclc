# ============================================================================
# fetch_emt_signatures.R
#
# Reproducibly materialize the canonical EMT gene sets vendored under
# Data/emt_signatures/. Run ONCE in your R environment (it needs network +
# {readxl} + {msigdbr}); the resulting TSVs are committed to the repo so all
# downstream scoring is reproducible offline.
#
#   source("R/fetch_emt_signatures.R")
#   fetch_emt_signatures()                 # writes Data/emt_signatures/*.tsv
#
# Sources (see Data/emt_signatures/PROVENANCE.md for full citations):
#   76GS / KS / MLR gene lists  -> Cancer-Systems-Biology-Lab (Jolly lab),
#       github.com/Cancer-Systems-Biology-Lab/EMT_Scoring_scRNA (branch master)
#       which packages the original Byers 2013 / Tan 2014 / George 2017 sets.
#   Hallmark EMT                 -> MSigDB via {msigdbr} (version-pinned).
#
# IMPORTANT: gene lists are downloaded from authoritative sources, never
# hand-typed. If a download fails, the corresponding TSV is simply not written
# and the scorer in R/emt_scoring.R errors with an actionable message.
# ============================================================================

#' Download and vendor the canonical EMT gene sets
#'
#' @param out_dir Destination directory. Default "Data/emt_signatures".
#' @param ks_variant "tumor" (default; for patient cohorts) or "cellLine"
#'   (for cell-line panels). Only the selected variant is downloaded and written
#'   to `tan_ks_*.tsv` (the files [load_emt_signatures()] reads); re-run with the
#'   other value to switch.
#' @param branch Git branch of the source repo. Default "master".
#' @return (Invisibly) a named list of the gene sets written.
#' @export
fetch_emt_signatures <- function(out_dir = "Data/emt_signatures",
                                 ks_variant = c("tumor", "cellLine"),
                                 branch = "master") {
  ks_variant <- match.arg(ks_variant)
  for (p in c("readxl", "msigdbr")) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required. Install it first (install.packages('%s')).",
                   p, p), call. = FALSE)
    }
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  base <- sprintf("https://raw.githubusercontent.com/Cancer-Systems-Biology-Lab/EMT_Scoring_scRNA/%s/Gene_signatures",
                  branch)
  tmp <- tempfile("emt_sig_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  # Join URL parts with "/" -- NEVER file.path(), whose separator is OS-specific
  # (backslashes on Windows would corrupt the URL).
  url_join <- function(...) paste(c(...), collapse = "/")

  dl <- function(url, dest) {
    status <- tryCatch(
      utils::download.file(url, dest, mode = "wb", quiet = TRUE),
      error = function(e) stop(sprintf("Download failed for %s: %s",
                                       url, conditionMessage(e)), call. = FALSE)
    )
    if (!identical(as.integer(status), 0L) || !file.exists(dest) ||
        file.size(dest) == 0) {
      stop(sprintf("Download failed (status %s, empty/missing file) for %s",
                   status, url), call. = FALSE)
    }
    dest
  }
  # Correct known OCR / Excel-date corruptions present in the upstream source
  # spreadsheets, restoring the intended valid HGNC symbols. Keys are matched
  # against whitespace-stripped tokens. These are EXACT-TOKEN corrections, never
  # heuristic O->0 / l->1 substitution, which would wrongly mangle legitimate
  # symbols such as ELMO3, MYO6, NQO1, ENO2. Each entry is documented in
  # Data/emt_signatures/PROVENANCE.md.
  symbol_fixes <- c(
    "Fl1R"     = "F11R",      # OCR l/1+space: F11R (junctional adhesion mol. A)
    "Clorf116" = "C1orf116",  # OCR letter l -> digit 1
    "Clorf172" = "C1orf172",  # OCR letter l + space -> digit 1
    "TMEM3OB"  = "TMEM30B",   # OCR letter O -> digit 0
    "TACSTDI"  = "TACSTD1",   # OCR letter I -> digit 1
    "38961"    = "SEPT1"      # Excel date corruption: serial 38961 = 1-Sep
  )
  clean_symbols <- function(genes) {
    g <- gsub("\\s+", "", trimws(as.character(genes)))   # HGNC symbols have no spaces
    g <- stats::na.omit(g)
    idx <- g %in% names(symbol_fixes)
    g[idx] <- unname(symbol_fixes[g[idx]])
    bad <- g == "" | grepl("^[0-9]+$", g)   # empty or uncorrected numeric artifacts
    if (any(bad)) {
      warning(sprintf("Dropping %d non-symbol token(s): %s", sum(bad),
                      paste(unique(g[bad]), collapse = ", ")), call. = FALSE)
    }
    unique(g[!bad])
  }
  write_genes <- function(genes, file) {
    genes <- clean_symbols(genes)
    utils::write.table(data.frame(gene = genes), file.path(out_dir, file),
                       sep = "\t", row.names = FALSE, quote = FALSE)
    genes
  }
  # Run one source resiliently: a failure warns and returns the fallback so the
  # rest still vendor (matches the "failed download -> TSV not written" note).
  safely <- function(expr, label, fallback = NULL) {
    tryCatch(expr, error = function(e) {
      warning(sprintf("%s skipped: %s", label, conditionMessage(e)), call. = FALSE)
      fallback
    })
  }

  result <- list()

  # ---- 76GS (Byers 2013): genes in column 2, with header ------------------
  result$gs_76 <- safely({
    f76 <- dl(url_join(base, "76GS", "EMT_signature_76GS.xlsx"),
              file.path(tmp, "76gs.xlsx"))
    write_genes(as.data.frame(readxl::read_excel(f76, col_names = TRUE))[[2]],
                "byers_76gs.tsv")
  }, "76GS")

  # ---- KS (Tan 2014): col1 gene, col2 category (Epi/Mes), no header -------
  ks <- safely({
    ks_file <- sprintf("EM_gene_signature_%s_KS.xlsx", ks_variant)
    fks <- dl(url_join(base, "KS", ks_file), file.path(tmp, "ks.xlsx"))
    xks <- as.data.frame(readxl::read_excel(fks, col_names = FALSE))
    cat_col <- tolower(trimws(as.character(xks[[2]])))
    list(epi = write_genes(xks[[1]][grepl("^epi", cat_col)], "tan_ks_epithelial.tsv"),
         mes = write_genes(xks[[1]][grepl("^mes", cat_col)], "tan_ks_mesenchymal.tsv"))
  }, "KS", fallback = list(epi = NULL, mes = NULL))
  result$ks_epithelial  <- ks$epi
  result$ks_mesenchymal <- ks$mes

  # ---- MLR (George 2017): predictor gene list (coefficients NOT included) -
  result$mlr <- safely({
    fmlr <- dl(url_join(base, "MLR", "genes_for_EMT_score.txt"),
               file.path(tmp, "mlr.txt"))
    mlr_genes <- clean_symbols(scan(fmlr, sep = "\n", what = character(), quiet = TRUE))
    utils::write.table(data.frame(gene = mlr_genes, role = "predictor"),
                       file.path(out_dir, "george_mlr.tsv"),
                       sep = "\t", row.names = FALSE, quote = FALSE)
    message("Note: george_mlr.tsv holds the MLR predictor genes only. score_mlr() ",
            "additionally needs the published class coefficients (coef_* columns); ",
            "until added it stays disabled. See PROVENANCE.md.")
    mlr_genes
  }, "MLR")

  # ---- Hallmark EMT (MSigDB via msigdbr) ----------------------------------
  result$hallmark <- safely({
    hm <- tryCatch(
      msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
      error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "H")
    )
    write_genes(hm$gene_symbol[hm$gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"],
                "hallmark_emt.tsv")
  }, "Hallmark")

  message(sprintf("Wrote EMT signatures to %s/ : 76GS=%d, KS(epi=%d, mes=%d), MLR=%d, Hallmark=%d",
                  out_dir, length(result$gs_76), length(result$ks_epithelial),
                  length(result$ks_mesenchymal), length(result$mlr),
                  length(result$hallmark)))
  invisible(result)
}
