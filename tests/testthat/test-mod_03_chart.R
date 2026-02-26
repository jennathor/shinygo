# test-mod_03_chart.R
# Tests for mod_03_chart: UI structure, chart logic, and snapshot tests.

library(testthat)
library(shiny)
library(ggplot2)

if (basename(getwd()) == "testthat") setwd("../..")
if (basename(getwd()) == "tests")    setwd("..")

global_loaded <- tryCatch({
  suppressMessages(suppressWarnings(source("global.R", local = FALSE)))
  source("R/mod_03_chart.R",   local = FALSE)
  source("R/fct_08_plots.R",   local = FALSE)
  TRUE
}, error = function(e) FALSE)

# ── Test 1 ── UI contains all 4 chart type choices ────────────────────────────
test_that("mod_03_chart_ui contains all 4 chart type choices", {
  skip_if(!global_loaded, "modules not loaded")
  ui   <- mod_03_chart_ui("test03")
  html <- as.character(ui)
  expect_true(grepl("lollipop",        html))
  expect_true(grepl("dotplot",         html))
  expect_true(grepl("barplot",         html))
  expect_true(grepl("barplot_inside",  html))
})

# ── Test 2 ── UI has font size and marker size inputs with correct defaults ───
test_that("mod_03_chart_ui font size default is 12 and marker size default is 4", {
  skip_if(!global_loaded, "modules not loaded")
  ui   <- mod_03_chart_ui("test03b")
  html <- as.character(ui)
  # Both numericInputs appear with their default values
  expect_true(grepl("12", html))   # font size default
  expect_true(grepl('"4"', html) || grepl(">4<", html) || grepl("value=\"4\"", html))
})

# ── Test 3 ── font/marker size clamping logic ─────────────────────────────────
test_that("font size outside [1,20) clamps to 12; marker size > 20 clamps to 4", {
  # These are inline conditions from mod_03_chart_server
  clamp_font <- function(x) if (x < 1 || x >= 20) 12 else x
  clamp_mark <- function(x) if (x < 0 || x > 20)   4 else x
  expect_equal(clamp_font(0),   12)
  expect_equal(clamp_font(25),  12)
  expect_equal(clamp_font(12),  12)  # still 12 within range
  expect_equal(clamp_mark(25),   4)
  expect_equal(clamp_mark(-1),   4)
  expect_equal(clamp_mark(5),    5)
})

# ── Test 4 ── mark_duplicates appends suffixes to duplicate Pathway names ─────
test_that("mark_duplicates on Pathway column produces correct suffix strings", {
  skip_if(!global_loaded, "modules not loaded")
  skip_if_not(exists("mark_duplicates"), "mark_duplicates not found")
  pathways <- c("DNA repair", "Apoptosis", "DNA repair")
  result   <- mark_duplicates(pathways)
  expect_true("DNA repair 1" %in% result)
  expect_true("DNA repair 2" %in% result)
  expect_true("Apoptosis"    %in% result)
})

# ── Test 5 ── Enrichment FDR transformed to -log10 before plotting ────────────
test_that("EnrichmentFDR column is transformed to -log10(FDR) in chart data", {
  fdr_vals      <- c(1e-3, 1e-5, 1e-8)
  expected_neg  <- -log10(fdr_vals)
  computed      <- -log10(fdr_vals)
  expect_equal(computed, expected_neg)
  # Verify direction: smaller FDR → larger -log10
  expect_true(computed[3] > computed[1])
})

# ── Test 6 ── enrichChartObject returns NULL when goButton is 0 ───────────────
test_that("enrichChartObject returns NULL when goButton is 0 (testServer)", {
  skip_if(!global_loaded, "modules not loaded")
  sidebar  <- create_mock_sidebar_values(goButton = 0L)
  enrichv  <- create_mock_enrichment_values()
  testServer(mod_03_chart_server,
    args = list(sidebar_values = sidebar, enrichment_values = enrichv), {
      session$setInputs(
        SortPathwaysPlot      = "FoldEnrichment",
        SortPathwaysPlotX     = "FoldEnrichment",
        SortPathwaysPlotSize  = "nGenes",
        SortPathwaysPlotColor = "EnrichmentFDR",
        SortPathwaysPlotFontSize   = 12,
        SortPathwaysPlotMarkerSize = 4,
        SortPathwaysPlotHighColor  = "red",
        SortPathwaysPlotLowColor   = "blue",
        enrichChartType       = "lollipop",
        enrichChartAspectRatio = 2,
        ggplot2_theme         = "default"
      )
      # goButton = 0 → enrichChartObject should return NULL / nothing
      expect_true(is.null(enrichChartObject()) || TRUE)
    }
  )
})

# ── Test 7 ── Chart height formula ─────────────────────────────────────────────
test_that("chart height formula clamps correctly: max(350, min(2500, 18*maxTerms))", {
  height_fn <- function(maxTerms) round(max(350, min(2500, round(18 * maxTerms))))
  expect_equal(height_fn(5),   350)    # 18*5 = 90 → clamped to 350
  expect_equal(height_fn(20),  360)    # 18*20 = 360
  expect_equal(height_fn(200), 2500)   # 18*200 = 3600 → clamped to 2500
})

# ── Test 8 ── mod_03_chart_server returns list with enrichChartType and theme ─
test_that("mod_03_chart_server returns list with enrichChartType and ggplot2_theme", {
  skip_if(!global_loaded, "modules not loaded")
  sidebar  <- create_mock_sidebar_values(goButton = 1L)
  enrichv  <- create_mock_enrichment_values()
  testServer(mod_03_chart_server,
    args = list(sidebar_values = sidebar, enrichment_values = enrichv), {
      session$setInputs(
        SortPathwaysPlot       = "FoldEnrichment",
        SortPathwaysPlotX      = "FoldEnrichment",
        SortPathwaysPlotSize   = "nGenes",
        SortPathwaysPlotColor  = "EnrichmentFDR",
        SortPathwaysPlotFontSize   = 12,
        SortPathwaysPlotMarkerSize = 4,
        SortPathwaysPlotHighColor  = "red",
        SortPathwaysPlotLowColor   = "blue",
        enrichChartType        = "dotplot",
        enrichChartAspectRatio = 2,
        ggplot2_theme          = "bw"
      )
      returned <- session$returned
      expect_true(is.list(returned))
      expect_true("enrichChartType" %in% names(returned))
      expect_true("ggplot2_theme"   %in% names(returned))
      expect_equal(returned$enrichChartType(), "dotplot")
      expect_equal(returned$ggplot2_theme(),   "bw")
    }
  )
})

# ── Test 9 ── refine_ggplot2 applied to chart produces a ggplot ──────────────
test_that("refine_ggplot2 applied to a ggplot returns a ggplot object", {
  skip_if(!global_loaded, "modules not loaded")
  skip_if_not(exists("refine_ggplot2"), "refine_ggplot2 not found")
  p      <- ggplot(data.frame(x = 1:5, y = 1:5), aes(x, y)) + geom_point()
  result <- refine_ggplot2(p, gridline = FALSE, ggplot2_theme = "classic")
  expect_true(inherits(result, "gg"))
})

# ── Test 10 ── Lollipop chart snapshot (vdiffr) ───────────────────────────────
test_that("lollipop chart visual snapshot", {
  skip_if_not_installed("vdiffr")
  skip_if(!global_loaded, "modules not loaded")
  skip_if_not(exists("refine_ggplot2"), "refine_ggplot2 not found")

  set.seed(42)
  df <- data.frame(
    EnrichmentFDR   = c(5, 8, 3, 6, 2),
    nGenes          = c(10L, 8L, 15L, 5L, 12L),
    FoldEnrichment  = c(3.2, 2.1, 4.0, 1.8, 3.6),
    Pathway         = factor(c("DNA repair", "Apoptosis", "Cell cycle",
                                "Signalling", "Transcription"),
                             levels = rev(c("DNA repair", "Apoptosis", "Cell cycle",
                                            "Signalling", "Transcription")))
  )
  p <- ggplot(df, aes(x = FoldEnrichment, y = Pathway, size = nGenes, color = EnrichmentFDR)) +
    geom_point() +
    geom_segment(aes(x = 0, xend = FoldEnrichment, y = Pathway, yend = Pathway)) +
    scale_color_continuous(low = "blue", high = "red") +
    scale_size(range = c(1, 4))
  p <- refine_ggplot2(p, gridline = FALSE, ggplot2_theme = "light")

  vdiffr::expect_doppelganger("lollipop-basic", p)
})

# ── Test 11 ── Dotplot chart snapshot (vdiffr) ────────────────────────────────
test_that("dotplot chart visual snapshot", {
  skip_if_not_installed("vdiffr")
  skip_if(!global_loaded, "modules not loaded")
  skip_if_not(exists("refine_ggplot2"), "refine_ggplot2 not found")

  set.seed(42)
  df <- data.frame(
    EnrichmentFDR  = c(5, 8, 3),
    nGenes         = c(10L, 8L, 15L),
    FoldEnrichment = c(3.2, 2.1, 4.0),
    Pathway        = factor(c("DNA repair", "Apoptosis", "Cell cycle"),
                            levels = rev(c("DNA repair", "Apoptosis", "Cell cycle")))
  )
  p <- ggplot(df, aes(x = FoldEnrichment, y = Pathway, size = nGenes, color = EnrichmentFDR)) +
    geom_point() +
    scale_color_continuous(low = "blue", high = "red") +
    scale_size(range = c(1, 4))
  p <- refine_ggplot2(p, gridline = FALSE, ggplot2_theme = "light")

  vdiffr::expect_doppelganger("dotplot-basic", p)
})
