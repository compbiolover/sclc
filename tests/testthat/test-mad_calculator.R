# tests/testthat/test-mad_calculator.R

test_that("mad_calculator handles valid input correctly", {
  # Setup
  test_matrix <- matrix(
    rnorm(100),
    nrow = 10,
    dimnames = list(paste0("gene_", 1:10))
  )

  # Test basic functionality
  result <- mad_calculator(test_matrix)
  expect_type(result, "double")
  expect_length(result, 10)
  expect_true(all(result >= 0))
  expect_true(abs(sum(result) - 1) < 1e-10)

  # Test without normalization
  unnorm_result <- mad_calculator(test_matrix, normalize = FALSE)
  expect_true(any(unnorm_result != result))
})

test_that("mad_calculator validates input correctly", {
  expect_error(
    mad_calculator("not a matrix"),
    "must be a matrix"
  )

  # Create character matrix without warning
  char_matrix <- matrix(
    rep(letters[1:4], 1),
    nrow = 2,
    ncol = 2
  )
  expect_error(
    mad_calculator(char_matrix),
    "must contain numeric values"
  )

  expect_error(
    mad_calculator(matrix(1:4, 2, 2), normalize = "yes"),
    "must be TRUE or FALSE"
  )

  expect_error(
    mad_calculator(matrix(1:4, 2, 2), genes_are_rows = "yes"),
    "must be TRUE or FALSE"
  )
})

test_that("mad_calculator handles min_expression filtering correctly", {
  # Create test matrix with known values
  test_matrix <- matrix(
    c(1, 1, 1, 1,    # MAD = 0
      1, 2, 3, 4,    # MAD ≈ 1.48
      5, 5, 5, 5),   # MAD = 0
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("gene_", 1:3), paste0("cell_", 1:4))
  )

  # Debug: Print raw MAD values
  raw_mads <- apply(test_matrix, 1, mad)
  print("Raw MAD values:")
  print(raw_mads)

  # Test with min_expression = 0 (should keep all values)
  result_0 <- mad_calculator(test_matrix, min_expression = 0)
  print("Result with min_expression = 0:")
  print(result_0)
  expect_length(result_0, 3)  # Should keep all genes

  # Test with min_expression = 0.5 (should keep only gene 2)
  result_1 <- mad_calculator(test_matrix, min_expression = 0.5)
  print("Result with min_expression = 0.5:")
  print(result_1)
  expect_length(result_1, 1)
  expect_equal(names(result_1), "gene_2")

  # Test with min_expression = 2 (should keep no genes)
  result_2 <- mad_calculator(test_matrix, min_expression = 2)
  print("Result with min_expression = 2:")
  print(result_2)
  expect_length(result_2, 0)

  # Test that transposed matrix gives same results
  result_t <- mad_calculator(t(test_matrix),
                             min_expression = 0.5,
                             genes_are_rows = FALSE)
  expect_equal(result_1, result_t)
})

test_that("mad_calculator handles different gene orientations correctly", {
  # Create test data with more variation
  set.seed(123)
  genes <- 5
  cells <- 10
  row_matrix <- matrix(
    rnorm(genes * cells),
    nrow = genes,
    dimnames = list(paste0("gene_", 1:genes), paste0("cell_", 1:cells))
  )

  # Calculate both ways
  result_rows <- mad_calculator(row_matrix, genes_are_rows = TRUE)
  result_cols <- mad_calculator(t(row_matrix), genes_are_rows = FALSE)

  # Results should be identical
  expect_equal(result_rows, result_cols)
})