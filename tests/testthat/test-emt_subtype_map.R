# tests/testthat/test-emt_subtype_map.R
#
# Tests for the SCLC subtype caller, NE score, and EMT->subtype mapping. Uses
# synthetic data with a known structure; no external gene lists required.

if (!exists("call_sclc_subtype", mode = "function")) {
  source(testthat::test_path("..", "..", "R", "emt_subtype_map.R"))
}

# A/N/P expression matrix where POU2F3 has a high CONSTANT baseline, so a naive
# raw-argmax caller would wrongly label everything "P"; z-scoring fixes it.
make_subtype_matrix <- function() {
  set.seed(5)
  genes <- c("ascl1", "neurod1", "pou2f3", "yap1")   # lowercase (scl_common style)
  samp <- c(paste0("A", 1:3), paste0("N", 1:3), paste0("P", 1:3), "UT1")
  m <- matrix(5 + stats::rnorm(length(genes) * length(samp), 0, 0.1),
              nrow = length(genes), dimnames = list(genes, samp))
  m["pou2f3", ] <- m["pou2f3", ] + 20          # high constant baseline (raw-argmax trap)
  m["ascl1",   1:3] <- m["ascl1",   1:3] + 10  # A samples
  m["neurod1", 4:6] <- m["neurod1", 4:6] + 10  # N samples
  m["pou2f3",  7:9] <- m["pou2f3",  7:9] + 10  # P samples
  m["yap1",    10]  <- m["yap1",    10]  + 20  # SMARCA4-UT-like: high YAP1, low A/N/P
  m
}

test_that("call_sclc_subtype uses z-scores, not raw magnitude", {
  m <- make_subtype_matrix()
  res <- suppressWarnings(call_sclc_subtype(m))        # UT1 triggers the SMARCA4-UT note
  expect_equal(res$subtype[res$sample == "A1"], "A")  # despite POU2F3 baseline
  expect_equal(res$subtype[res$sample == "N2"], "N")
  expect_equal(res$subtype[res$sample == "P3"], "P")
  # a raw-argmax caller would have called A1 "P" (POU2F3 ~25 > ASCL1 ~15)
  expect_gt(m["pou2f3", "A1"], m["ascl1", "A1"])
})

test_that("call_sclc_subtype flags SMARCA4-UT (high YAP1, low A/N/P)", {
  m <- make_subtype_matrix()
  res <- suppressWarnings(call_sclc_subtype(m))
  expect_true(res$smarca4_ut_flag[res$sample == "UT1"])
  expect_false(any(res$smarca4_ut_flag[res$sample %in% c("A1", "N2", "P3")]))
})

test_that("call_sclc_subtype validates inputs", {
  m <- make_subtype_matrix()
  expect_error(call_sclc_subtype(m[, 1:2]), ">= 3 samples")
  expect_error(call_sclc_subtype(m[c("ascl1", "yap1"), ]), "marker")  # NEUROD1/POU2F3 missing
})

test_that("call_sclc_subtype requires markers named exactly A/N/P", {
  m <- make_subtype_matrix()
  expect_error(call_sclc_subtype(m, markers = list(A = "ascl1", N = "neurod1")), "A, N, P")
  expect_error(call_sclc_subtype(m, markers = list(X = "ascl1", Y = "neurod1", Z = "pou2f3")), "A, N, P")
})

test_that("call_sclc_subtype accepts a sparse single-cell matrix (densifies markers only)", {
  skip_if_not_installed("Matrix")
  m <- make_subtype_matrix()
  sp <- methods::as(Matrix::Matrix(m, sparse = TRUE), "CsparseMatrix")
  res <- suppressWarnings(call_sclc_subtype(sp))            # same answers as the dense call
  expect_equal(res$subtype[res$sample == "A1"], "A")
  expect_equal(res$subtype[res$sample == "N2"], "N")
  expect_equal(res$subtype[res$sample == "P3"], "P")
  expect_true(res$smarca4_ut_flag[res$sample == "UT1"])
})

test_that("map_emt_to_subtype binds cleanly when NE values are sparse", {
  m <- make_subtype_matrix()
  subt <- suppressWarnings(call_sclc_subtype(m))
  emt <- data.frame(sample = subt$sample,
                    consensus = seq(0, 1, length.out = nrow(subt)),
                    stringsAsFactors = FALSE)
  ne <- stats::setNames(rep(NA_real_, nrow(subt)), subt$sample)
  ne[1:2] <- c(0.1, 0.2)                              # only 2 non-NA -> every subtype <= 2
  res <- map_emt_to_subtype(emt, subt, ne = ne)       # must not error on rbind
  expect_true("cor_emt_ne" %in% names(res$by_subtype))
  expect_true(all(is.na(res$by_subtype$cor_emt_ne)))
})

test_that("ne_score requires a template and is oriented NE-high", {
  set.seed(1)
  ne_genes <- paste0("NEG", 1:10); non_genes <- paste0("NONG", 1:10)
  genes <- c(ne_genes, non_genes)
  tmpl <- data.frame(gene = genes,
                     ne_ref = c(rep(10, 10), rep(1, 10)),
                     nonne_ref = c(rep(1, 10), rep(10, 10)))
  # one NE-like sample, one non-NE-like sample
  expr <- cbind(
    NE_like  = c(rep(10, 10), rep(1, 10)) + stats::rnorm(20, 0, 0.1),
    NON_like = c(rep(1, 10), rep(10, 10)) + stats::rnorm(20, 0, 0.1)
  )
  rownames(expr) <- genes
  expect_error(ne_score(expr, ne_template = NULL), "NE template")
  s <- ne_score(expr, ne_template = tmpl)
  expect_gt(s[["NE_like"]], s[["NON_like"]])
  expect_gt(s[["NE_like"]], 0)
})

test_that("vendored Zhang NE template loads and scores correctly", {
  ne_path <- testthat::test_path("..", "..", "Data", "sclc_signatures", "zhang_ne_50.tsv")
  skip_if(!file.exists(ne_path), "zhang_ne_50.tsv not vendored")
  tmpl <- load_ne_template(ne_path)
  expect_equal(nrow(tmpl), 50)
  expect_true(all(c("gene", "ne_ref", "nonne_ref") %in% names(tmpl)))
  # a sample matching the NE reference scores higher than one matching non-NE
  set.seed(9)
  expr <- cbind(NE = tmpl$ne_ref + stats::rnorm(50, 0, 0.2),
                NON = tmpl$nonne_ref + stats::rnorm(50, 0, 0.2))
  rownames(expr) <- tmpl$gene
  s <- ne_score(expr, ne_template = tmpl)
  expect_gt(s[["NE"]], s[["NON"]])
  expect_gt(s[["NE"]], 0); expect_lt(s[["NON"]], 0)
})

test_that("ne_score_singlecell scores NE-like cells higher (per-cell, UCell)", {
  skip_if_not_installed("UCell")
  set.seed(3)
  ne_g <- paste0("NEG", 1:10); non_g <- paste0("NONG", 1:10)
  tmpl <- data.frame(gene = c(ne_g, non_g),
                     class = c(rep("NE", 10), rep("non_NE", 10)), stringsAsFactors = FALSE)
  genes <- c(ne_g, non_g, paste0("BG", 1:30))
  m <- matrix(stats::rpois(length(genes) * 40, 1), nrow = length(genes),
              dimnames = list(genes, paste0("c", 1:40)))
  m[ne_g, 1:20]    <- m[ne_g, 1:20]    + 8     # cells 1-20 are NE-like
  m[non_g, 21:40]  <- m[non_g, 21:40]  + 8     # cells 21-40 are non-NE-like
  s <- ne_score_singlecell(m, ne_template = tmpl, method = "UCell")
  expect_length(s, 40)
  expect_gt(mean(s[1:20]), mean(s[21:40]))     # NE-like cells score higher
  expect_error(ne_score_singlecell(m, ne_template = NULL), "NE template")
})

test_that("map_emt_to_subtype joins and summarizes by subtype", {
  m <- make_subtype_matrix()
  subt <- suppressWarnings(call_sclc_subtype(m))
  emt <- data.frame(sample = subt$sample,
                    consensus = seq(0, 1, length.out = nrow(subt)),
                    emt_state = "hybrid", stringsAsFactors = FALSE)
  ne_vec <- stats::setNames(seq(1, -1, length.out = nrow(subt)), subt$sample)
  res <- map_emt_to_subtype(emt, subt, ne = ne_vec)
  expect_true(all(c("per_sample", "by_subtype") %in% names(res)))
  expect_true(all(c("A", "N", "P") %in% res$by_subtype$subtype))
  expect_true("cor_emt_ne" %in% names(res$by_subtype))
  expect_equal(nrow(res$per_sample), nrow(subt))
})

test_that("map_emt_to_subtype validates the ne argument", {
  m <- make_subtype_matrix()
  subt <- suppressWarnings(call_sclc_subtype(m))
  emt <- data.frame(sample = subt$sample,
                    consensus = seq(0, 1, length.out = nrow(subt)),
                    stringsAsFactors = FALSE)
  expect_error(map_emt_to_subtype(emt, subt, ne = seq_len(nrow(subt))), "named")  # unnamed vector
  expect_error(
    map_emt_to_subtype(emt, subt, ne = data.frame(s = subt$sample, x = 1)), "sample"  # wrong columns
  )
})

test_that("aggregate_pseudobulk sums/means cells per group (dense + sparse)", {
  m <- matrix(c(1,5, 2,6, 3,7, 4,8), nrow = 2,
              dimnames = list(c("g1","g2"), paste0("c", 1:4)))   # 2 genes x 4 cells
  groups <- c("A","A","B","B")
  pb <- aggregate_pseudobulk(m, groups, fun = "sum")
  expect_equal(colnames(pb), c("A","B"))
  expect_equal(pb["g1","A"], 3); expect_equal(pb["g1","B"], 7)
  expect_equal(pb["g2","A"], 11); expect_equal(pb["g2","B"], 15)
  pm <- aggregate_pseudobulk(m, groups, fun = "mean")
  expect_equal(pm["g1","A"], 1.5)                       # (1+2)/2
  # sparse input gives the same answer; NA group cells are dropped
  skip_if_not_installed("Matrix")
  sp <- methods::as(Matrix::Matrix(m, sparse = TRUE), "CsparseMatrix")
  expect_equal(as.numeric(aggregate_pseudobulk(sp, groups)["g1", ]), c(3, 7))
  gna <- c("A","A","B",NA)
  expect_equal(as.numeric(aggregate_pseudobulk(m, gna)["g1", ]), c(3, 3))  # c4 dropped -> B=3
  expect_error(aggregate_pseudobulk(m, c("A","B","C")), "must equal ncol")
})

test_that("emt_ne_coupling reports per-group correlation + dispersion", {
  set.seed(8)
  # group HET: strong negative EMT~NE with real spread; group HOM: tiny spread
  het <- data.frame(grp = "HET", emt = stats::rnorm(200), stringsAsFactors = FALSE)
  het$ne <- -0.8 * het$emt + stats::rnorm(200, 0, 0.3)
  hom <- data.frame(grp = "HOM", emt = stats::rnorm(200, 0, 0.02),
                    ne = stats::rnorm(200, 0, 0.02), stringsAsFactors = FALSE)
  pc <- rbind(het, hom)
  res <- emt_ne_coupling(pc, by = "grp", method = "pearson")
  expect_setequal(res$group, c("HET","HOM"))
  expect_lt(res$cor_emt_ne[res$group == "HET"], -0.5)  # strong negative where heterogeneous
  expect_gt(res$emt_sd[res$group == "HET"], res$emt_sd[res$group == "HOM"])  # HET more spread
  expect_true(all(c("n","cor_emt_ne","emt_sd","ne_sd") %in% names(res)))
  expect_error(emt_ne_coupling(pc, by = "missing"), "Missing")
})

test_that("call_sclc_subtype errors actionably on sparse + genes_are_rows=FALSE", {
  skip_if_not_installed("Matrix")
  m <- make_subtype_matrix()
  sp <- methods::as(Matrix::Matrix(m, sparse = TRUE), "CsparseMatrix")
  expect_error(call_sclc_subtype(sp, genes_are_rows = FALSE), "genes in rows")
})
