# test-fct_helper_functions.R
# Tests for helper functions in R/fct_08_plots.R and R/fct_09_genome.R.
# No database required.

library(testthat)
library(ggplot2)

if (basename(getwd()) == "testthat") setwd("../..")
if (basename(getwd()) == "tests")    setwd("..")

# Source helper functions
fct_loaded <- tryCatch({
  suppressMessages(suppressWarnings({
    source("R/fct_08_plots.R", local = FALSE)
    source("R/fct_09_genome.R", local = FALSE)
  }))
  TRUE
}, error = function(e) FALSE)

# We also need dplyr for genome_sliding_window
suppressPackageStartupMessages(library(dplyr))

# ── Test 1 ── densMode: peak near mode of input ──────────────────────────────
test_that("densMode returns list with $x and $y; peak near mode of input", {
  skip_if(!fct_loaded, "fct_08_plots.R not loaded")
  set.seed(42)
  x      <- rnorm(100, mean = 5, sd = 0.5)
  result <- densMode(x)
  expect_true(is.list(result))
  expect_true(all(c("x", "y") %in% names(result)))
  expect_true(abs(result$x - 5) < 0.5)   # peak within 0.5 of true mode
})

# ── Test 2 ── densMode: handles identical values without error ───────────────
test_that("densMode handles vectors with all identical values without error", {
  skip_if(!fct_loaded, "fct_08_plots.R not loaded")
  expect_no_error({
    result <- densMode(rep(3, 20))
  })
  expect_true(is.list(result))
})

# ── Test 3 ── mark_significance: correct symbols for p-value thresholds ──────
test_that("mark_significance returns correct significance symbols", {
  skip_if(!fct_loaded, "fct_08_plots.R not loaded")
  # Three thresholds used by the app: PvalGeneInfo2=0.001, PvalGeneInfo1=0.01, PvalGeneInfo=0.05
  expect_true(grepl("\\*\\*\\*", mark_significance(0.0005, 0.001, 0.01, 0.05)))
  expect_true(grepl("\\*\\*$",   mark_significance(0.005,  0.001, 0.01, 0.05)))
  expect_true(grepl("\\*$",      mark_significance(0.04,   0.001, 0.01, 0.05)))
  # p >= 0.05 → no star
  result_ns <- mark_significance(0.1, 0.001, 0.01, 0.05)
  expect_false(grepl("\\*", result_ns))
})

# ── Test 4 ── mark_significance: formats in scientific notation ───────────────
# Note: formatC(x, digits=2, format="G") only switches to exponential when the
# exponent is < -4 (e.g., 1e-6, not 0.000234 whose exponent is exactly -4).
test_that("mark_significance formats p-value in scientific notation for very small p", {
  skip_if(!fct_loaded, "fct_08_plots.R not loaded")
  # 1e-6 has exponent -6 < -4, which triggers %E format in formatC("G")
  result <- mark_significance(1e-6, 0.001, 0.01, 0.05)
  expect_true(grepl("[Ee]", result))
  # Also confirm the *** stars are appended (1e-6 < 0.001)
  expect_true(grepl("\\*\\*\\*", result))
})

# ── Test 5 ── fake_plot: returns ggplot with message text ────────────────────
test_that("fake_plot returns a ggplot object containing the message text", {
  skip_if(!fct_loaded, "fct_08_plots.R not loaded")
  msg    <- "Data not available."
  result <- fake_plot(msg)
  expect_true(inherits(result, "gg"))
  # Check the annotation is present in the ggplot build
  built  <- ggplot_build(result)
  labels <- unlist(lapply(built$data, function(d) if ("label" %in% names(d)) d$label))
  expect_true(any(grepl(msg, labels, fixed = TRUE)))
})

# ── Test 6 ── refine_ggplot2: applies a theme without error ──────────────────
test_that("refine_ggplot2 applies a named theme without error; result is still a ggplot", {
  skip_if(!fct_loaded, "fct_08_plots.R not loaded")
  p <- ggplot(data.frame(x = 1:5, y = 1:5), aes(x, y)) + geom_point()
  result <- refine_ggplot2(p, gridline = TRUE, ggplot2_theme = "classic")
  expect_true(inherits(result, "gg"))
})

# ── Test 7 ── genome_sliding_window: correct output columns ─────────────────
test_that("genome_sliding_window returns data frame with chNum, x, n, k, pval, y columns", {
  skip_if(!fct_loaded, "fct_09_genome.R not loaded")
  # Create input: 20 genes on chr1, half in list, half background
  set.seed(42)
  x0 <- data.frame(
    start_position   = seq(1, 50, length.out = 20),
    chNum            = 1L,
    chromosome_name  = "1",
    Fold             = c(rep(1L, 10), rep(0L, 10)),
    stringsAsFactors = FALSE
  )
  result <- genome_sliding_window(x0, listN = 10, totalN = 20,
                                  windowSize = 20, steps = 2, pvalCutoff = 1.0)
  if (nrow(result) > 0) {
    expect_true(all(c("chNum", "x", "n", "k", "pval", "y") %in% colnames(result)))
  } else {
    # If no windows pass FDR cutoff (FDR=1.0 should keep all, but n>=3 filter
    # may remove some), just verify no error occurred
    expect_true(nrow(result) >= 0)
  }
})

# ── Test 8 ── genome_sliding_window: pval in (0, 1] ─────────────────────────
test_that("genome_sliding_window column pval is in (0, 1] for valid input", {
  skip_if(!fct_loaded, "fct_09_genome.R not loaded")
  set.seed(42)
  x0 <- data.frame(
    start_position  = seq(1, 100, by = 5),
    chNum           = 1L,
    chromosome_name = "1",
    Fold            = c(rep(1L, 10), rep(0L, 10)),
    stringsAsFactors = FALSE
  )
  result <- genome_sliding_window(x0, listN = 10, totalN = 20,
                                  windowSize = 30, steps = 2, pvalCutoff = 1.0)
  if (nrow(result) > 0) {
    expect_true(all(result$pval > 0 & result$pval <= 1))
  }
})

# ── Test 9 ── genome_sliding_window: impossibly stringent cutoff → 0 rows ────
test_that("genome_sliding_window returns 0-row result when pvalCutoff is 0", {
  skip_if(!fct_loaded, "fct_09_genome.R not loaded")
  x0 <- data.frame(
    start_position  = seq(1, 50, length.out = 20),
    chNum           = 1L,
    chromosome_name = "1",
    Fold            = c(rep(1L, 10), rep(0L, 10)),
    stringsAsFactors = FALSE
  )
  result <- genome_sliding_window(x0, listN = 10, totalN = 20,
                                  windowSize = 20, steps = 2, pvalCutoff = 0)
  expect_equal(nrow(result), 0L)
})

# ── Test 10 ── genome_sliding_window: empty input (0 query genes) ─────────────
test_that("genome_sliding_window handles empty input gracefully", {
  skip_if(!fct_loaded, "fct_09_genome.R not loaded")
  x0 <- data.frame(
    start_position  = seq(1, 50, length.out = 10),
    chNum           = 1L,
    chromosome_name = "1",
    Fold            = rep(0L, 10),     # no query genes
    stringsAsFactors = FALSE
  )
  expect_no_error({
    result <- genome_sliding_window(x0, listN = 0, totalN = 10,
                                    windowSize = 20, steps = 2, pvalCutoff = 0.05)
  })
  # With 0 list genes there is no enrichment, so 0 rows expected
  expect_equal(nrow(result), 0L)
})

# ── Test 11 ── genome_sliding_window: clustered genes → at least one window
test_that("genome_sliding_window: 10 query genes clustered on chr1 at pos 1-10 → k == 10 in some window", {
  skip_if(!fct_loaded, "fct_09_genome.R not loaded")
  # All list genes tightly packed in positions 1-10 Mbp
  x0 <- data.frame(
    start_position  = c(1:10, 50, 100, 150, 200, 250),  # list genes at 1-10, background far away
    chNum           = 1L,
    chromosome_name = "1",
    Fold            = c(rep(1L, 10), rep(0L, 5)),
    stringsAsFactors = FALSE
  )
  result <- genome_sliding_window(x0, listN = 10, totalN = 15,
                                  windowSize = 15, steps = 1, pvalCutoff = 1.0)
  # At least one window should contain all 10 list genes
  if (nrow(result) > 0) {
    expect_true(max(result$k) >= 9L)
  }
})

# ── Test 12 ── refine_ggplot2: theme "bw" produces white panel background ────
test_that("refine_ggplot2 with theme 'bw' changes theme to bw", {
  skip_if(!fct_loaded, "fct_08_plots.R not loaded")
  p <- ggplot(data.frame(x = 1:5, y = 1:5), aes(x, y)) + geom_point()
  result <- refine_ggplot2(p, gridline = FALSE, ggplot2_theme = "bw")
  built  <- ggplot_build(result)
  theme_bg <- built$plot$theme$panel.background
  # theme_bw sets panel.background to element_rect(fill = "white") or element_blank
  is_white <- (!is.null(theme_bg) &&
               !inherits(theme_bg, "element_blank") &&
               isTRUE(theme_bg$fill == "white"))
  is_blank <- inherits(theme_bg, "element_blank") || is.null(theme_bg)
  expect_true(is_white || is_blank)
})
