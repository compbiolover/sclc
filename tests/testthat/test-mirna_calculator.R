# tests/testthat/test-mirna_calculator.R

source("../../R/mirna_calculator.R")

# =============================================================================
# Tests for validate_mirna_inputs
# =============================================================================

test_that("validate_mirna_inputs rejects invalid mirmap_path", {
  expect_error(
    validate_mirna_inputs(
      mirmap_path = 123,
      dbdemc_path = tempfile(),
      mirdb_path = NULL,
      mirtarbase_path = NULL,
      cancer_type = "lung cancer",
      status = "up",
      ts_context_threshold = -0.2,
      mirdb_score_threshold = 80,
      min_databases = 2
    ),
    "mirmap_path must be a single character string"
  )
})

test_that("validate_mirna_inputs rejects non-existent mirmap file", {
  expect_error(
    validate_mirna_inputs(
      mirmap_path = "/nonexistent/path.csv",
      dbdemc_path = tempfile(),
      mirdb_path = NULL,
      mirtarbase_path = NULL,
      cancer_type = "lung cancer",
      status = "up",
      ts_context_threshold = -0.2,
      mirdb_score_threshold = 80,
      min_databases = 2
    ),
    "miRmap file not found"
  )
})

test_that("validate_mirna_inputs rejects invalid status", {
  mirmap_f <- tempfile(fileext = ".csv")
  dbdemc_f <- tempfile(fileext = ".txt")
  writeLines("dummy", mirmap_f)
  writeLines("dummy", dbdemc_f)

  expect_error(
    validate_mirna_inputs(
      mirmap_path = mirmap_f,
      dbdemc_path = dbdemc_f,
      mirdb_path = NULL,
      mirtarbase_path = NULL,
      cancer_type = "lung cancer",
      status = "sideways",
      ts_context_threshold = -0.2,
      mirdb_score_threshold = 80,
      min_databases = 2
    ),
    "status must be 'up' or 'down'"
  )

  unlink(c(mirmap_f, dbdemc_f))
})

test_that("validate_mirna_inputs rejects invalid mirdb_score_threshold", {
  mirmap_f <- tempfile(fileext = ".csv")
  dbdemc_f <- tempfile(fileext = ".txt")
  writeLines("dummy", mirmap_f)
  writeLines("dummy", dbdemc_f)

  expect_error(
    validate_mirna_inputs(
      mirmap_path = mirmap_f,
      dbdemc_path = dbdemc_f,
      mirdb_path = NULL,
      mirtarbase_path = NULL,
      cancer_type = "lung cancer",
      status = "up",
      ts_context_threshold = -0.2,
      mirdb_score_threshold = 150,
      min_databases = 2
    ),
    "mirdb_score_threshold must be a single numeric value between 0 and 100"
  )

  unlink(c(mirmap_f, dbdemc_f))
})

test_that("validate_mirna_inputs rejects invalid min_databases", {
  mirmap_f <- tempfile(fileext = ".csv")
  dbdemc_f <- tempfile(fileext = ".txt")
  writeLines("dummy", mirmap_f)
  writeLines("dummy", dbdemc_f)

  expect_error(
    validate_mirna_inputs(
      mirmap_path = mirmap_f,
      dbdemc_path = dbdemc_f,
      mirdb_path = NULL,
      mirtarbase_path = NULL,
      cancer_type = "lung cancer",
      status = "up",
      ts_context_threshold = -0.2,
      mirdb_score_threshold = 80,
      min_databases = 5
    ),
    "min_databases must be 1, 2, or 3"
  )

  unlink(c(mirmap_f, dbdemc_f))
})

test_that("validate_mirna_inputs accepts valid inputs", {
  mirmap_f <- tempfile(fileext = ".csv")
  dbdemc_f <- tempfile(fileext = ".txt")
  writeLines("dummy", mirmap_f)
  writeLines("dummy", dbdemc_f)

  expect_invisible(
    validate_mirna_inputs(
      mirmap_path = mirmap_f,
      dbdemc_path = dbdemc_f,
      mirdb_path = NULL,
      mirtarbase_path = NULL,
      cancer_type = "lung cancer",
      status = "up",
      ts_context_threshold = -0.2,
      mirdb_score_threshold = 80,
      min_databases = 2
    )
  )

  unlink(c(mirmap_f, dbdemc_f))
})

test_that("validate_mirna_inputs rejects non-existent mirdb file when provided", {
  mirmap_f <- tempfile(fileext = ".csv")
  dbdemc_f <- tempfile(fileext = ".txt")
  writeLines("dummy", mirmap_f)
  writeLines("dummy", dbdemc_f)

  expect_error(
    validate_mirna_inputs(
      mirmap_path = mirmap_f,
      dbdemc_path = dbdemc_f,
      mirdb_path = "/nonexistent/mirdb.txt",
      mirtarbase_path = NULL,
      cancer_type = "lung cancer",
      status = "up",
      ts_context_threshold = -0.2,
      mirdb_score_threshold = 80,
      min_databases = 2
    ),
    "miRDB file not found"
  )

  unlink(c(mirmap_f, dbdemc_f))
})


# =============================================================================
# Tests for build_consensus
# =============================================================================

test_that("build_consensus requires at least one database", {
  expect_error(
    build_consensus(
      targetscan_df = data.frame(),
      mirdb_df = NULL,
      mirtarbase_df = NULL,
      verbose = FALSE
    ),
    "No interactions found in any database"
  )
})

test_that("build_consensus warns when only 1 database but min_databases > 1", {
  ts_df <- data.frame(
    mirna = c("hsa-miR-1", "hsa-miR-1"),
    gene = c("TP53", "BRCA1"),
    context_score = c(-0.3, -0.5),
    source = "targetscan",
    stringsAsFactors = FALSE
  )

  expect_warning(
    build_consensus(ts_df, mirdb_df = NULL, mirtarbase_df = NULL,
                    min_databases = 2, verbose = FALSE),
    "Only 1 database provided"
  )
})

test_that("build_consensus correctly counts database support", {
  ts_df <- data.frame(
    mirna = c("hsa-miR-1", "hsa-miR-1", "hsa-miR-2"),
    gene = c("TP53", "BRCA1", "TP53"),
    context_score = c(-0.3, -0.5, -0.4),
    source = "targetscan",
    stringsAsFactors = FALSE
  )

  mirdb_df <- data.frame(
    mirna = c("hsa-miR-1", "hsa-miR-2"),
    gene = c("TP53", "MYC"),
    mirdb_score = c(90, 85),
    source = "mirdb",
    stringsAsFactors = FALSE
  )

  result <- build_consensus(ts_df, mirdb_df, mirtarbase_df = NULL,
                            min_databases = 2, verbose = FALSE)

  # hsa-miR-1:TP53 should be in both databases
  tp53_row <- result[result$mirna == "hsa-miR-1" & result$gene == "TP53", ]
  expect_equal(tp53_row$n_databases, 2)

  # hsa-miR-1:BRCA1 only in TargetScan -> should be filtered out
  expect_false("BRCA1" %in% result$gene[result$mirna == "hsa-miR-1"])

  # hsa-miR-2:MYC only in miRDB -> should be filtered out
  expect_false("MYC" %in% result$gene)
})

test_that("build_consensus passes experimentally validated interactions regardless of count", {
  ts_df <- data.frame(
    mirna = "hsa-miR-1",
    gene = "TP53",
    context_score = -0.3,
    source = "targetscan",
    stringsAsFactors = FALSE
  )

  mirtarbase_df <- data.frame(
    mirna = c("hsa-miR-1", "hsa-miR-99"),
    gene = c("TP53", "EGFR"),
    evidence_type = c("Luciferase reporter assay", "Western blot"),
    source = "mirtarbase",
    stringsAsFactors = FALSE
  )

  result <- build_consensus(ts_df, mirdb_df = NULL, mirtarbase_df,
                            min_databases = 2, verbose = FALSE)

  # EGFR only in miRTarBase but experimentally validated -> should pass
  expect_true("EGFR" %in% result$gene)
  egfr_row <- result[result$gene == "EGFR", ]
  expect_true(egfr_row$has_experimental_validation)
})

test_that("build_consensus confidence_score increases with more database support", {
  ts_df <- data.frame(
    mirna = c("hsa-miR-1", "hsa-miR-1"),
    gene = c("TP53", "BRCA1"),
    context_score = c(-0.3, -0.3),
    source = "targetscan",
    stringsAsFactors = FALSE
  )

  mirdb_df <- data.frame(
    mirna = "hsa-miR-1",
    gene = "TP53",
    mirdb_score = 90,
    source = "mirdb",
    stringsAsFactors = FALSE
  )

  # Use min_databases = 1 so both pass
  result <- build_consensus(ts_df, mirdb_df, mirtarbase_df = NULL,
                            min_databases = 1, verbose = FALSE)

  tp53_score <- result$confidence_score[result$gene == "TP53"]
  brca1_score <- result$confidence_score[result$gene == "BRCA1"]

  # TP53 supported by 2 databases should score higher than BRCA1 with 1
  expect_gt(tp53_score, brca1_score)
})


# =============================================================================
# Tests for compute_gene_ranking
# =============================================================================

test_that("compute_gene_ranking returns normalized scores summing to 1", {
  consensus_df <- data.frame(
    mirna = c("hsa-miR-1", "hsa-miR-1", "hsa-miR-2", "hsa-miR-2"),
    gene = c("TP53", "BRCA1", "TP53", "MYC"),
    n_databases = c(2, 2, 3, 1),
    databases = c("ts,mdb", "ts,mdb", "ts,mdb,mtb", "ts"),
    has_experimental_validation = c(FALSE, FALSE, TRUE, FALSE),
    confidence_score = c(0.5, 0.4, 0.8, 0.3),
    stringsAsFactors = FALSE
  )

  ranking <- compute_gene_ranking(consensus_df, verbose = FALSE)

  expect_type(ranking, "double")
  expect_true(abs(sum(ranking) - 1) < 1e-10)
  expect_true(all(ranking > 0))
  expect_true(!is.unsorted(rev(ranking)))  # sorted descending
})

test_that("compute_gene_ranking ranks genes with more/stronger interactions higher", {
  consensus_df <- data.frame(
    mirna = c("hsa-miR-1", "hsa-miR-2", "hsa-miR-3", "hsa-miR-1"),
    gene = c("TP53", "TP53", "TP53", "BRCA1"),
    n_databases = c(3, 3, 2, 1),
    databases = c("ts,mdb,mtb", "ts,mdb,mtb", "ts,mdb", "ts"),
    has_experimental_validation = c(TRUE, TRUE, FALSE, FALSE),
    confidence_score = c(0.9, 0.8, 0.6, 0.3),
    stringsAsFactors = FALSE
  )

  ranking <- compute_gene_ranking(consensus_df, verbose = FALSE)

  # TP53 targeted by 3 miRNAs with high confidence should rank above BRCA1
  expect_gt(ranking["TP53"], ranking["BRCA1"])
})

test_that("compute_gene_ranking returns named vector", {
  consensus_df <- data.frame(
    mirna = c("hsa-miR-1"),
    gene = c("TP53"),
    n_databases = 2,
    databases = "ts,mdb",
    has_experimental_validation = FALSE,
    confidence_score = 0.5,
    stringsAsFactors = FALSE
  )

  ranking <- compute_gene_ranking(consensus_df, verbose = FALSE)
  expect_true(!is.null(names(ranking)))
  expect_equal(names(ranking), "TP53")
  expect_equal(unname(ranking), 1.0)
})


# =============================================================================
# Tests for load_mirdb
# =============================================================================

test_that("load_mirdb filters by score threshold", {
  # Create temp miRDB file
  mirdb_f <- tempfile(fileext = ".txt")
  writeLines(c(
    "hsa-miR-1\tTP53\t95",
    "hsa-miR-1\tBRCA1\t60",
    "hsa-miR-2\tMYC\t85",
    "hsa-miR-3\tEGFR\t45"
  ), mirdb_f)

  result <- load_mirdb(mirdb_f, mirnas = c("hsa-miR-1", "hsa-miR-2", "hsa-miR-3"),
                       score_threshold = 80, verbose = FALSE)

  expect_equal(nrow(result), 2)
  expect_true("TP53" %in% result$gene)
  expect_true("MYC" %in% result$gene)
  expect_false("BRCA1" %in% result$gene)
  expect_false("EGFR" %in% result$gene)

  unlink(mirdb_f)
})

test_that("load_mirdb filters by miRNA list", {
  mirdb_f <- tempfile(fileext = ".txt")
  writeLines(c(
    "hsa-miR-1\tTP53\t95",
    "hsa-miR-99\tBRCA1\t90"
  ), mirdb_f)

  result <- load_mirdb(mirdb_f, mirnas = c("hsa-miR-1"),
                       score_threshold = 80, verbose = FALSE)

  expect_equal(nrow(result), 1)
  expect_equal(result$gene, "TP53")

  unlink(mirdb_f)
})

test_that("load_mirdb rejects file with fewer than 3 columns", {
  mirdb_f <- tempfile(fileext = ".txt")
  writeLines(c("hsa-miR-1\tTP53"), mirdb_f)

  expect_error(
    load_mirdb(mirdb_f, mirnas = c("hsa-miR-1"), verbose = FALSE),
    "at least 3 columns"
  )

  unlink(mirdb_f)
})


# =============================================================================
# Tests for load_mirtarbase
# =============================================================================

test_that("load_mirtarbase filters strong evidence", {
  mtb_f <- tempfile(fileext = ".csv")
  mtb_data <- data.frame(
    miRNA = c("hsa-miR-1", "hsa-miR-1", "hsa-miR-2"),
    Target.Gene = c("TP53", "BRCA1", "MYC"),
    Experiments = c("Luciferase reporter assay", "Microarray", "Western blot"),
    stringsAsFactors = FALSE
  )
  write.csv(mtb_data, mtb_f, row.names = FALSE)

  result <- load_mirtarbase(mtb_f, mirnas = c("hsa-miR-1", "hsa-miR-2"),
                            strong_only = TRUE, verbose = FALSE)

  expect_equal(nrow(result), 2)
  expect_true("TP53" %in% result$gene)
  expect_true("MYC" %in% result$gene)
  expect_false("BRCA1" %in% result$gene)  # Microarray = weak evidence

  unlink(mtb_f)
})

test_that("load_mirtarbase includes all evidence when strong_only = FALSE", {
  mtb_f <- tempfile(fileext = ".csv")
  mtb_data <- data.frame(
    miRNA = c("hsa-miR-1", "hsa-miR-1"),
    Target.Gene = c("TP53", "BRCA1"),
    Experiments = c("Luciferase reporter assay", "Microarray"),
    stringsAsFactors = FALSE
  )
  write.csv(mtb_data, mtb_f, row.names = FALSE)

  result <- load_mirtarbase(mtb_f, mirnas = c("hsa-miR-1"),
                            strong_only = FALSE, verbose = FALSE)

  expect_equal(nrow(result), 2)

  unlink(mtb_f)
})

test_that("load_mirtarbase deduplicates mirna-gene pairs", {
  mtb_f <- tempfile(fileext = ".csv")
  mtb_data <- data.frame(
    miRNA = c("hsa-miR-1", "hsa-miR-1"),
    Target.Gene = c("TP53", "TP53"),
    Experiments = c("Luciferase reporter assay", "Western blot"),
    stringsAsFactors = FALSE
  )
  write.csv(mtb_data, mtb_f, row.names = FALSE)

  result <- load_mirtarbase(mtb_f, mirnas = c("hsa-miR-1"),
                            strong_only = TRUE, verbose = FALSE)

  expect_equal(nrow(result), 1)

  unlink(mtb_f)
})


# =============================================================================
# Tests for load_targetscan_bulk
# =============================================================================

test_that("load_targetscan_bulk loads predictions with context++ score filtering", {
  # Create temp TargetScan bulk file with context++ scores
  ts_f <- tempfile(fileext = ".txt")
  ts_data <- data.frame(
    `miR Family` = c("hsa-miR-1/206", "hsa-miR-1/206", "hsa-miR-1/206"),
    `Gene Symbol` = c("TP53", "BRCA1", "MYC"),
    `Cumulative weighted context++ score` = c(-0.45, -0.10, -0.30),
    `Species ID` = c(9606, 9606, 9606),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  write.table(ts_data, ts_f, sep = "\t", row.names = FALSE, quote = FALSE)

  # Create miR family mapping file
  fam_f <- tempfile(fileext = ".txt")
  fam_data <- data.frame(
    `miR family` = c("hsa-miR-1/206", "hsa-miR-1/206"),
    `MiRBase ID` = c("hsa-miR-1-3p", "hsa-miR-206"),
    `Species ID` = c(9606, 9606),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  write.table(fam_data, fam_f, sep = "\t", row.names = FALSE, quote = FALSE)

  result <- load_targetscan_bulk(
    targetscan_path = ts_f,
    ts_family_path = fam_f,
    mirnas = c("hsa-miR-1-3p"),
    context_threshold = -0.2,
    verbose = FALSE
  )

  # TP53 (-0.45) and MYC (-0.30) pass threshold; BRCA1 (-0.10) does not
  expect_true("TP53" %in% result$gene)
  expect_true("MYC" %in% result$gene)
  expect_false("BRCA1" %in% result$gene)
  expect_equal(unique(result$source), "targetscan")

  unlink(c(ts_f, fam_f))
})

test_that("load_targetscan_bulk returns empty data frame when no miRNAs match", {
  ts_f <- tempfile(fileext = ".txt")
  ts_data <- data.frame(
    `miR Family` = c("hsa-miR-99"),
    `Gene Symbol` = c("TP53"),
    `Species ID` = c(9606),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  write.table(ts_data, ts_f, sep = "\t", row.names = FALSE, quote = FALSE)

  result <- load_targetscan_bulk(
    targetscan_path = ts_f,
    ts_family_path = NULL,
    mirnas = c("hsa-miR-1"),
    verbose = FALSE
  )

  expect_equal(nrow(result), 0)

  unlink(ts_f)
})

test_that("load_targetscan_bulk works without family mapping (direct miRNA match)", {
  ts_f <- tempfile(fileext = ".txt")
  ts_data <- data.frame(
    `miR Family` = c("hsa-miR-1-3p", "hsa-miR-1-3p"),
    `Gene Symbol` = c("TP53", "BRCA1"),
    `Species ID` = c(9606, 9606),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  write.table(ts_data, ts_f, sep = "\t", row.names = FALSE, quote = FALSE)

  result <- load_targetscan_bulk(
    targetscan_path = ts_f,
    ts_family_path = NULL,
    mirnas = c("hsa-miR-1-3p"),
    verbose = FALSE
  )

  # No context++ score column -> all predictions kept
  expect_equal(nrow(result), 2)
  expect_true(all(c("TP53", "BRCA1") %in% result$gene))

  unlink(ts_f)
})

# =============================================================================
# Tests for filter_mirnas_from_geo (pre-computed DE results mode)
# =============================================================================

test_that("filter_mirnas_from_geo loads pre-computed DE results and filters correctly", {
  de_f <- tempfile(fileext = ".csv")
  de_data <- data.frame(
    mirna = c("hsa-miR-1", "hsa-miR-21", "hsa-miR-200c", "hsa-miR-99"),
    log2fc = c(2.5, -1.8, 0.3, -3.0),
    adj_p_value = c(0.001, 0.01, 0.2, 0.001),
    stringsAsFactors = FALSE
  )
  write.csv(de_data, de_f, row.names = FALSE)

  result <- filter_mirnas_from_geo(
    de_results_path = de_f,
    fdr_threshold = 0.05,
    min_log2fc = 1.0,
    verbose = FALSE
  )

  # miR-1 (log2FC=2.5, FDR=0.001) -> passes (up)
  # miR-21 (log2FC=-1.8, FDR=0.01) -> passes (down)
  # miR-200c (log2FC=0.3) -> fails min_log2fc
  # miR-99 (log2FC=-3.0, FDR=0.001) -> passes (down)
  expect_equal(nrow(result), 3)
  expect_true("hsa-miR-1" %in% result$mirna)
  expect_true("hsa-miR-21" %in% result$mirna)
  expect_true("hsa-miR-99" %in% result$mirna)
  expect_false("hsa-miR-200c" %in% result$mirna)

  # Check direction assignment
  expect_equal(result$direction[result$mirna == "hsa-miR-1"], "up")
  expect_equal(result$direction[result$mirna == "hsa-miR-21"], "down")

  # Check abs_log2fc
  expect_equal(result$abs_log2fc[result$mirna == "hsa-miR-1"], 2.5)

  unlink(de_f)
})

test_that("filter_mirnas_from_geo respects mirna_remove parameter", {
  de_f <- tempfile(fileext = ".csv")
  de_data <- data.frame(
    mirna = c("hsa-miR-1", "hsa-miR-21"),
    log2fc = c(2.5, -1.8),
    adj_p_value = c(0.001, 0.01),
    stringsAsFactors = FALSE
  )
  write.csv(de_data, de_f, row.names = FALSE)

  result <- filter_mirnas_from_geo(
    de_results_path = de_f,
    fdr_threshold = 0.05,
    min_log2fc = 1.0,
    mirna_remove = "hsa-miR-1",
    verbose = FALSE
  )

  expect_equal(nrow(result), 1)
  expect_equal(result$mirna, "hsa-miR-21")

  unlink(de_f)
})

test_that("filter_mirnas_from_geo returns both directions", {
  de_f <- tempfile(fileext = ".csv")
  de_data <- data.frame(
    mirna = c("hsa-miR-1", "hsa-miR-21", "hsa-miR-99"),
    log2fc = c(3.0, -2.5, -1.5),
    adj_p_value = c(0.001, 0.001, 0.01),
    stringsAsFactors = FALSE
  )
  write.csv(de_data, de_f, row.names = FALSE)

  result <- filter_mirnas_from_geo(
    de_results_path = de_f,
    fdr_threshold = 0.05,
    min_log2fc = 1.0,
    verbose = FALSE
  )

  # All three pass: both up and down are included
  expect_equal(nrow(result), 3)
  expect_true("up" %in% result$direction)
  expect_true("down" %in% result$direction)

  unlink(de_f)
})

test_that("filter_mirnas_from_geo rejects malformed DE results file", {
  de_f <- tempfile(fileext = ".csv")
  bad_data <- data.frame(
    gene = c("TP53"),
    score = c(0.5),
    stringsAsFactors = FALSE
  )
  write.csv(bad_data, de_f, row.names = FALSE)

  expect_error(
    filter_mirnas_from_geo(de_results_path = de_f, verbose = FALSE),
    "must have columns"
  )

  unlink(de_f)
})

test_that("filter_mirnas_from_geo respects max_mirnas limit", {
  de_f <- tempfile(fileext = ".csv")
  de_data <- data.frame(
    mirna = paste0("hsa-miR-", 1:50),
    log2fc = seq(5, 1, length.out = 50),
    adj_p_value = seq(0.001, 0.04, length.out = 50),
    stringsAsFactors = FALSE
  )
  write.csv(de_data, de_f, row.names = FALSE)

  result <- filter_mirnas_from_geo(
    de_results_path = de_f,
    fdr_threshold = 0.05,
    min_log2fc = 1.0,
    max_mirnas = 10,
    verbose = FALSE
  )

  expect_lte(nrow(result), 10)

  unlink(de_f)
})


test_that("validate_mirna_inputs rejects non-existent targetscan_path", {
  mirmap_f <- tempfile(fileext = ".csv")
  dbdemc_f <- tempfile(fileext = ".txt")
  writeLines("dummy", mirmap_f)
  writeLines("dummy", dbdemc_f)

  expect_error(
    validate_mirna_inputs(
      mirmap_path = mirmap_f,
      dbdemc_path = dbdemc_f,
      mirdb_path = NULL,
      mirtarbase_path = NULL,
      targetscan_path = "/nonexistent/ts.txt",
      cancer_type = "lung cancer",
      status = "up",
      ts_context_threshold = -0.2,
      mirdb_score_threshold = 80,
      min_databases = 2
    ),
    "TargetScan bulk file not found"
  )

  unlink(c(mirmap_f, dbdemc_f))
})
