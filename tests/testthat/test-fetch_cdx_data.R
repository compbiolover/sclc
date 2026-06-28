# tests/testthat/test-fetch_cdx_data.R
#
# Tests for the WS4 CDX loader's PARSING + 10x reading (no network). The sample
# fixture is the real GSE138267 (title, treatment) table, so the base-model /
# condition / pairing logic is checked against the actual cohort structure.

if (!exists("cdx_classify_samples", mode = "function")) {
  source(testthat::test_path("..", "..", "R", "fetch_cdx_data.R"))
}

# The genuine 25-sample GSE138267 sheet (title + treatment), as observed on GEO.
real_cdx_samples <- function() {
  data.frame(
    title = c("SC16.LB17028","SC16.LB17029","SC39.LB17018","SC39.LB17019",
              "SC39.LB17020","SC49.LB17023","SC49.LB17024","SC49.LB17025",
              "SC4.LB17009","SC4.LB17013","SC4.LB17014","SC55-1.LB17001",
              "SC55-2.LB17002","SC55-3.LB17026","SC4_LY2606368.LB17003",
              "SC4_LY2606368.LB17016","SC4_Talazoparib.LB17004",
              "SC4_Talazoparib.LB17011","SC4_Talazoparib.LB17015","SC55.LB19001",
              "SC68_cr.LB19003","SC68.LB19002","SC75.LB19004","SC53cis","SC53"),
    treatment = c(rep("untreated", 12), "cisplatin", "cisplatin, relapsed",
                  rep("prexasertib, relapsed", 2), rep("talazoparib, relapsed", 3),
                  "untreated", "cisplatin, relapsed", "vehicle-treated",
                  "untreated", "cisplatin, relapsed", "vehicle-treated"),
    stringsAsFactors = FALSE
  )
}

# ---- base model extraction (the critical SC-stem join) ---------------------

test_that("cdx_base_model takes the SC-number stem, not the dotted prefix", {
  expect_equal(
    cdx_base_model(c("SC4.LB17009", "SC4_Talazoparib.LB17004", "SC4_LY2606368.LB17003")),
    c("SC4", "SC4", "SC4")                       # all one model despite drug tags
  )
  expect_equal(cdx_base_model(c("SC49.LB17023", "SC4.LB17009")), c("SC49", "SC4"))  # SC49 != SC4
  expect_equal(cdx_base_model(c("SC53cis", "SC53", "SC55-3.LB17026", "SC68_cr.LB19003")),
               c("SC53", "SC53", "SC55", "SC68"))
  expect_error(cdx_base_model(c("SC4.LB1", "weird")), "SC<number> stem")
})

# ---- condition mapping -----------------------------------------------------

test_that("cdx_condition maps untreated/vehicle->sensitive, relapsed->resistant, else NA", {
  cc <- cdx_condition(c("untreated", "vehicle-treated", "cisplatin, relapsed",
                        "talazoparib, relapsed", "cisplatin"))
  expect_equal(as.character(cc),
               c("sensitive", "sensitive", "resistant", "resistant", NA))
  expect_identical(levels(cc), c("sensitive", "resistant"))
})

# ---- full classification against the real cohort ---------------------------

test_that("cdx_classify_samples reproduces the real 4-paired / 4-unpaired structure", {
  cls <- cdx_classify_samples(real_cdx_samples())
  expect_equal(sort(unique(cls$model)),
               c("SC16","SC39","SC4","SC49","SC53","SC55","SC68","SC75"))
  # exactly four models carry both a sensitive and a resistant sample
  expect_equal(sort(unique(cls$model[cls$paired])), c("SC4","SC53","SC55","SC68"))
  # the lone on-treatment (not-relapsed) cisplatin sample is excluded
  excl <- cls[cls$excluded, ]
  expect_equal(excl$title, "SC55-2.LB17002")
  expect_true(is.na(cls$condition[cls$title == "SC55-2.LB17002"]))
  # SC4 has 3 sensitive + 5 resistant (untreated vs the two relapsed drug arms)
  sc4 <- cls[cls$model == "SC4", ]
  expect_equal(sum(sc4$condition == "sensitive", na.rm = TRUE), 3)
  expect_equal(sum(sc4$condition == "resistant", na.rm = TRUE), 5)
  # sensitive-only models are not paired
  expect_false(any(cls$paired[cls$model %in% c("SC16","SC39","SC49","SC75")]))
})

test_that("cdx_classify_samples validates columns", {
  expect_error(cdx_classify_samples(data.frame(x = 1)), "Missing column")
})

# ---- 10x triplet reader (synthetic, no network) ----------------------------

test_that("read_10x_triplet reads MM, collapses duplicate symbols, prefixes cells", {
  skip_if_not_installed("Matrix")
  d <- file.path(tempdir(), "triplet_test"); dir.create(d, showWarnings = FALSE)
  # 3 gene rows (two share symbol GENEA), 2 cells.
  writeLines(c("ENSG1\tGENEA", "ENSG2\tGENEB", "ENSG3\tGENEA"),
             file.path(d, "genes.tsv"))
  writeLines(c("BC1-1", "BC2-1"), file.path(d, "barcodes.tsv"))
  writeLines(c("%%MatrixMarket matrix coordinate integer general", "%",
               "3 2 4",
               "1 1 5",    # GENEA(row1), cell1
               "3 1 2",    # GENEA(row3), cell1  -> collapses with row1
               "1 2 1",    # GENEA(row1), cell2
               "2 2 3"),   # GENEB,       cell2
             file.path(d, "matrix.mtx"))
  m <- read_10x_triplet(d, cell_prefix = "SAMP")
  expect_equal(rownames(m), c("GENEA", "GENEB"))            # collapsed + ordered
  expect_equal(colnames(m), c("SAMP|BC1-1", "SAMP|BC2-1"))  # prefixed, unique
  expect_equal(as.numeric(m["GENEA", "SAMP|BC1-1"]), 7)     # 5 + 2 summed
  expect_equal(as.numeric(m["GENEA", "SAMP|BC2-1"]), 1)
  expect_equal(as.numeric(m["GENEB", "SAMP|BC2-1"]), 3)
})

test_that("read_10x_triplet also reads the dense .txt layout (LB19 batch)", {
  skip_if_not_installed("Matrix")
  d <- file.path(tempdir(), "dense_test"); dir.create(d, showWarnings = FALSE)
  writeLines(c("ENSG1\tGENEA", "ENSG2\tGENEB", "ENSG3\tGENEA"),
             file.path(d, "genes.txt"))
  writeLines(c("BC1-1", "BC2-1"), file.path(d, "barcodes.txt"))
  # dense genes x cells table, no header: same data as the MM case above.
  writeLines(c("5\t1", "0\t3", "2\t0"), file.path(d, "matrix.txt"))
  m <- read_10x_triplet(d, cell_prefix = "SAMP")
  expect_equal(rownames(m), c("GENEA", "GENEB"))
  expect_equal(colnames(m), c("SAMP|BC1-1", "SAMP|BC2-1"))
  expect_equal(as.numeric(m["GENEA", "SAMP|BC1-1"]), 7)     # 5 + 2 (dup collapse)
  expect_equal(as.numeric(m["GENEB", "SAMP|BC2-1"]), 3)
})

test_that("read_10x_triplet finds sample-prefixed filenames (SC53 batch)", {
  skip_if_not_installed("Matrix")
  d <- file.path(tempdir(), "prefixed_test"); dir.create(d, showWarnings = FALSE)
  # Files prefixed with the sample name, as in the GSM4851xxx batch.
  writeLines(c("ENSG1\tGENEA", "ENSG2\tGENEB"), file.path(d, "SC53cis.genes.tsv"))
  writeLines(c("BC1-1", "BC2-1"), file.path(d, "SC53cis.barcodes.tsv"))
  writeLines(c("%%MatrixMarket matrix coordinate integer general", "%",
               "2 2 2", "1 1 4", "2 2 6"), file.path(d, "SC53cis.matrix.mtx"))
  m <- read_10x_triplet(d, cell_prefix = "SC53cis")
  expect_equal(rownames(m), c("GENEA", "GENEB"))
  expect_equal(as.numeric(m["GENEA", "SC53cis|BC1-1"]), 4)
  expect_equal(as.numeric(m["GENEB", "SC53cis|BC2-1"]), 6)
})

test_that("read_10x_triplet reads a gzipped MatrixMarket matrix", {
  skip_if_not_installed("Matrix")
  d <- file.path(tempdir(), "gz_test"); unlink(d, recursive = TRUE)
  dir.create(d, showWarnings = FALSE)
  writeLines(c("ENSG1\tGENEA", "ENSG2\tGENEB"), file.path(d, "genes.tsv"))
  writeLines(c("BC1", "BC2"), file.path(d, "barcodes.tsv"))
  con <- gzfile(file.path(d, "matrix.mtx.gz"), "w")    # only a gzipped matrix present
  writeLines(c("%%MatrixMarket matrix coordinate integer general", "%",
               "2 2 2", "1 1 9", "2 2 4"), con); close(con)
  m <- read_10x_triplet(d)
  expect_equal(as.numeric(m["GENEA", "BC1"]), 9)
  expect_equal(as.numeric(m["GENEB", "BC2"]), 4)
})
