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
#'   (for cell-line panels). Both are downloaded; this selects which becomes
#'   the default `tan_ks_*.tsv` read by [load_emt_signatures()].
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

  dl <- function(url, dest) {
    utils::download.file(url, dest, mode = "wb", quiet = TRUE)
    dest
  }
  write_genes <- function(genes, file) {
    genes <- unique(stats::na.omit(trimws(as.character(genes))))
    genes <- genes[genes != ""]
    utils::write.table(data.frame(gene = genes), file.path(out_dir, file),
                       sep = "\t", row.names = FALSE, quote = FALSE)
    genes
  }

  result <- list()

  # ---- 76GS (Byers 2013): genes in column 2, with header ------------------
  f76 <- dl(file.path(base, "76GS", "EMT_signature_76GS.xlsx"),
            file.path(tmp, "76gs.xlsx"))
  x76 <- as.data.frame(readxl::read_excel(f76, col_names = TRUE))
  genes76 <- x76[[2]]
  result$gs_76 <- write_genes(genes76, "byers_76gs.tsv")

  # ---- KS (Tan 2014): col1 gene, col2 category (Epi/Mes), no header -------
  ks_file <- sprintf("EM_gene_signature_%s_KS.xlsx", ks_variant)
  fks <- dl(file.path(base, "KS", ks_file), file.path(tmp, "ks.xlsx"))
  xks <- as.data.frame(readxl::read_excel(fks, col_names = FALSE))
  cat_col <- tolower(trimws(as.character(xks[[2]])))
  epi <- xks[[1]][grepl("^epi", cat_col)]
  mes <- xks[[1]][grepl("^mes", cat_col)]
  result$ks_epithelial  <- write_genes(epi, "tan_ks_epithelial.tsv")
  result$ks_mesenchymal <- write_genes(mes, "tan_ks_mesenchymal.tsv")

  # ---- MLR (George 2017): predictor gene list (coefficients NOT included) -
  fmlr <- dl(file.path(base, "MLR", "genes_for_EMT_score.txt"),
             file.path(tmp, "mlr.txt"))
  mlr_genes <- scan(fmlr, sep = "\n", what = "vector", quiet = TRUE)
  mlr_genes <- unique(trimws(mlr_genes)); mlr_genes <- mlr_genes[mlr_genes != ""]
  utils::write.table(
    data.frame(gene = mlr_genes, role = "predictor"),
    file.path(out_dir, "george_mlr.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE)
  result$mlr <- mlr_genes
  message("Note: george_mlr.tsv holds the MLR predictor genes only. score_mlr() ",
          "additionally needs the published class coefficients (coef_* columns); ",
          "until added it stays disabled. See PROVENANCE.md.")

  # ---- Hallmark EMT (MSigDB via msigdbr) ----------------------------------
  hm <- tryCatch(
    msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
    error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  )
  hm_emt <- unique(hm$gene_symbol[hm$gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"])
  result$hallmark <- write_genes(hm_emt, "hallmark_emt.tsv")

  message(sprintf("Wrote EMT signatures to %s/ : 76GS=%d, KS(epi=%d, mes=%d), MLR=%d, Hallmark=%d",
                  out_dir, length(result$gs_76), length(result$ks_epithelial),
                  length(result$ks_mesenchymal), length(result$mlr),
                  length(result$hallmark)))
  invisible(result)
}
