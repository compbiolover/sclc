# tests/testthat/test-fetch_chan_data.R
#
# Tests for the Chan-atlas (external validation) loader: the pure treatment->group
# mapping, and the CSR .h5ad reader against a tiny synthetic file (no network).

if (!exists("chan_treatment_group", mode = "function")) {
  source(testthat::test_path("..", "..", "R", "fetch_chan_data.R"))
}

test_that("chan_treatment_group maps Naive vs everything-else, keeping NA", {
  g <- chan_treatment_group(c("Naive", "Platinum Doublet",
                              "Platinum Doublet,Immunotherapy", NA, "Naive"))
  expect_equal(as.character(g), c("naive", "treated", "treated", NA, "naive"))
  expect_identical(levels(g), c("naive", "treated"))
})

test_that("read_h5ad_counts reads CSR, maps symbols, collapses duplicates", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("Matrix")
  f <- tempfile(fileext = ".h5ad"); on.exit(unlink(f), add = TRUE)
  rhdf5::h5createFile(f)
  rhdf5::h5createGroup(f, "raw"); rhdf5::h5createGroup(f, "raw/X")
  # 2 cells x 3 genes, CSR: cell1 -> {gene0:5, gene2:2}; cell2 -> {gene1:3}
  rhdf5::h5write(c(5, 2, 3),              f, "raw/X/data")
  rhdf5::h5write(as.integer(c(0, 2, 1)), f, "raw/X/indices")
  rhdf5::h5write(as.integer(c(0, 2, 3)), f, "raw/X/indptr")
  fid <- rhdf5::H5Fopen(f); gid <- rhdf5::H5Gopen(fid, "raw/X")
  rhdf5::h5writeAttribute("csr_matrix", gid, "encoding-type")
  rhdf5::h5writeAttribute(as.integer(c(2, 3)), gid, "shape")
  rhdf5::H5Gclose(gid); rhdf5::H5Fclose(fid)
  rhdf5::h5createGroup(f, "var")
  rhdf5::h5write(c("GENEA", "GENEB", "GENEA"), f, "var/feature_name")  # dup GENEA
  rhdf5::h5createGroup(f, "obs")
  rhdf5::h5write(c("c1", "c2"), f, "obs/_index")

  m <- read_h5ad_counts(f, layer = "raw")
  expect_equal(rownames(m), c("GENEA", "GENEB"))      # collapsed, in first-seen order
  expect_equal(colnames(m), c("c1", "c2"))
  expect_equal(as.numeric(m["GENEA", "c1"]), 7)       # 5 (gene0) + 2 (gene2)
  expect_equal(as.numeric(m["GENEB", "c1"]), 0)
  expect_equal(as.numeric(m["GENEB", "c2"]), 3)
  # subsetting to one cell works
  m2 <- read_h5ad_counts(f, keep_cells = 2L, layer = "raw")
  expect_equal(colnames(m2), "c2")
  expect_equal(as.numeric(m2["GENEB", "c2"]), 3)
})

test_that("read_h5ad_counts resolves the CELLxGENE obs `_index` attribute", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("Matrix")
  f <- tempfile(fileext = ".h5ad"); on.exit(unlink(f), add = TRUE)
  rhdf5::h5createFile(f)
  rhdf5::h5createGroup(f, "raw"); rhdf5::h5createGroup(f, "raw/X")
  rhdf5::h5write(c(5, 3), f, "raw/X/data")
  rhdf5::h5write(as.integer(c(0, 1)), f, "raw/X/indices")
  rhdf5::h5write(as.integer(c(0, 1, 2)), f, "raw/X/indptr")
  fid <- rhdf5::H5Fopen(f); gid <- rhdf5::H5Gopen(fid, "raw/X")
  rhdf5::h5writeAttribute("csr_matrix", gid, "encoding-type")
  rhdf5::h5writeAttribute(as.integer(c(2, 2)), gid, "shape")
  rhdf5::H5Gclose(gid); rhdf5::H5Fclose(fid)
  rhdf5::h5createGroup(f, "var"); rhdf5::h5write(c("GENEA", "GENEB"), f, "var/feature_name")
  # cell ids live in an obs column named by the /obs `_index` attribute (here "Cell")
  rhdf5::h5createGroup(f, "obs"); rhdf5::h5write(c("cellA", "cellB"), f, "obs/Cell")
  fid <- rhdf5::H5Fopen(f); gid <- rhdf5::H5Gopen(fid, "obs")
  rhdf5::h5writeAttribute("Cell", gid, "_index")
  rhdf5::H5Gclose(gid); rhdf5::H5Fclose(fid)

  m <- read_h5ad_counts(f, layer = "raw")
  expect_equal(colnames(m), c("cellA", "cellB"))     # used the attribute, not a literal _index
})
