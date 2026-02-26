# test-mock_contracts.R
# Verifies that mock constructors stay in sync with real function output shapes.
# No database required. Reads .__real_names attributes from mock helpers and
# compares against the actual column names of the returned objects.

library(testthat)

# ── Test 1 ── enrichment result x-frame columns match .__real_names ─────────
test_that("create_mock_enrichment_result()$x has exactly the columns in .__real_names", {
  result <- create_mock_enrichment_result(n = 5)
  expected_cols <- attr(result, ".__real_names")
  expect_equal(colnames(result$x), expected_cols)
})

# ── Test 2 ── gene info full columns match .__real_names ────────────────────
test_that("create_mock_gene_info_full() has exactly the columns in .__real_names", {
  df <- create_mock_gene_info_full(n = 10)
  expected_cols <- attr(df, ".__real_names")
  actual_cols   <- colnames(df)
  expect_equal(actual_cols, expected_cols)
})

# ── Test 3 ── conversion table data columns match .__real_names ─────────────
test_that("create_mock_conversion_table_data() has exactly the columns in .__real_names", {
  df <- create_mock_conversion_table_data(n = 5)
  expected_cols <- attr(df, ".__real_names")
  actual_cols   <- colnames(df)
  expect_equal(actual_cols, expected_cols)
})

# ── Test 4 ── enrichment groupings column count ─────────────────────────────
test_that("create_mock_enrichment_result()$groupings has 2 columns (Group, Pathway)", {
  result <- create_mock_enrichment_result(n = 5)
  # The real FindOverlap() groupings slot has at minimum Pathway/description column
  # and a count column. Our mock has Group + Pathway = 2 columns.
  expect_equal(ncol(result$groupings), 2L)
  expect_true(all(c("Group", "Pathway") %in% colnames(result$groupings)))
})

# ── Test 5 ── mock constructors return objects of the expected class ─────────
test_that("all mock constructors return objects of expected class", {
  expect_true(is.list(create_mock_enrichment_result()))
  expect_true(is.list(create_mock_sidebar_values()))
  expect_true(is.list(create_mock_enrichment_values()))
  expect_true(is.data.frame(create_mock_gene_info_full()))
  expect_true(is.data.frame(create_mock_conversion_table_data()))
  expect_true(is.list(create_mock_gmt()))
  expect_true(is.data.frame(create_mock_promoter_db()))
})

# ── Test 6 ── constructors are deterministic (same structure both times) ─────
test_that("mock constructors produce structurally identical objects on repeated calls", {
  r1 <- create_mock_enrichment_result(n = 5)
  r2 <- create_mock_enrichment_result(n = 5)
  expect_equal(colnames(r1$x), colnames(r2$x))
  expect_equal(nrow(r1$x),    nrow(r2$x))

  g1 <- create_mock_gene_info_full(n = 8)
  g2 <- create_mock_gene_info_full(n = 8)
  expect_equal(colnames(g1), colnames(g2))
  expect_equal(nrow(g1),     nrow(g2))

  c1 <- create_mock_conversion_table_data(n = 4)
  c2 <- create_mock_conversion_table_data(n = 4)
  expect_equal(colnames(c1), colnames(c2))
  expect_equal(nrow(c1),     nrow(c2))
})
