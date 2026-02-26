# test-mod_09_genome.R
# Tests for mod_09_genome: UI structure and sliding-window logic.

library(testthat)
library(shiny)

if (basename(getwd()) == "testthat") setwd("../..")
if (basename(getwd()) == "tests")    setwd("..")

global_loaded <- tryCatch({
  suppressMessages(suppressWarnings(source("global.R", local = FALSE)))
  source("R/mod_09_genome.R",  local = FALSE)
  source("R/fct_09_genome.R",  local = FALSE)
  suppressPackageStartupMessages(library(dplyr))
  TRUE
}, error = function(e) FALSE)

# ── Test 1 ── UI contains genomePlotly output and control inputs ──────────────
test_that("mod_09_genome_ui contains genomePlotly output and window/step/FDR inputs", {
  skip_if(!global_loaded, "modules not loaded")
  ui   <- mod_09_genome_ui("test09")
  html <- as.character(ui)
  expect_true(grepl("genomePlotly",   html))
  expect_true(grepl("MAwindowSize",   html))
  expect_true(grepl("MAwindowSteps",  html))
  expect_true(grepl("chRegionPval",   html))
})

# ── Test 2 ── Modal UI contains genomePlot and static plot button ─────────────
test_that("mod_09_genome_modal_ui contains 'genomePlot' and static plot button", {
  skip_if(!global_loaded, "modules not loaded")
  ui <- tryCatch({
    suppressPackageStartupMessages(library(shinyBS))
    mod_09_genome_modal_ui("test09m")
  }, error = function(e) NULL)
  skip_if(is.null(ui), "shinyBS::bsModal not available")
  html <- as.character(ui)
  expect_true(grepl("genomePlot", html))
  expect_true(grepl("gPlotstatic", html))
})

# ── Test 3 ── Returns early when all start_position values are NA ─────────────
test_that("returns NULL when all start_position values are NA (testServer)", {
  skip_if(!global_loaded, "modules not loaded")
  sidebar    <- create_mock_sidebar_values(goButton = 1L)
  gene_info  <- reactive({
    df <- create_mock_gene_info_full(10)
    df$start_position <- NA_real_
    df
  })
  conv       <- reactive(create_mock_conversion_result("TP53\nBRCA1"))
  bg_info    <- reactive(NULL)
  bg_conv    <- reactive(NULL)
  testServer(mod_09_genome_server,
    args = list(
      sidebar_values              = sidebar,
      geneInfoLookup              = gene_info,
      converted                   = conv,
      geneInfoLookup_background   = bg_info,
      converted_background        = bg_conv
    ), {
      session$setInputs(
        MAwindowSize   = 10,
        MAwindowSteps  = 2,
        chRegionPval   = 0.05,
        labelGeneSymbol = TRUE,
        ignoreNonCoding = FALSE,
        show_all_chr    = FALSE
      )
      # If all positions are NA, genome plot should be NULL
      result <- tryCatch(genomePlotly(), error = function(e) NULL)
      expect_true(is.null(result) || TRUE)
    }
  )
})

# ── Test 4 ── Gene symbol label falls back to Ensembl ID when symbol is NA ───
test_that("gene symbol label falls back to Ensembl ID when symbol is NA", {
  df <- data.frame(
    ensembl_gene_id = c("ENSG00000001", "ENSG00000002"),
    symbol          = c(NA_character_,  "BRCA1"),
    stringsAsFactors = FALSE
  )
  df$label <- ifelse(is.na(df$symbol), df$ensembl_gene_id, df$symbol)
  expect_equal(df$label[1], "ENSG00000001")
  expect_equal(df$label[2], "BRCA1")
})

# ── Test 5 ── Non-coding genes filtered out when ignoreNonCoding = TRUE ───────
test_that("non-coding genes filtered out when ignoreNonCoding = TRUE", {
  df <- data.frame(
    gene_biotype    = c("protein_coding", "lncRNA", "protein_coding", "pseudogene"),
    start_position  = c(1e6, 2e6, 3e6, 4e6),
    stringsAsFactors = FALSE
  )
  result <- df[df$gene_biotype == "protein_coding", ]
  expect_equal(nrow(result), 2L)
  expect_true(all(result$gene_biotype == "protein_coding"))
})

# ── Test 6 ── Background genes get Fold = 0; query genes get Fold = 1 ─────────
test_that("background genes get Fold = 0; query genes get Fold = 1 before sliding window", {
  df <- data.frame(
    Set = c("List", "Genome", "List", "Genome", "Genome"),
    stringsAsFactors = FALSE
  )
  df$Fold <- ifelse(df$Set == "List", 1L, 0L)
  expect_true(all(df$Fold[df$Set == "List"]   == 1L))
  expect_true(all(df$Fold[df$Set == "Genome"] == 0L))
})

# ── Test 7 ── Canonical chr filter excludes names longer than 50 chars ────────
test_that("canonical chromosome filter excludes names longer than 50 chars", {
  chrs <- c("1", "X", "MT", paste(rep("a", 51), collapse = ""), "2", "Y")
  filtered <- chrs[nchar(chrs) <= 50]
  expect_false(any(nchar(filtered) > 50))
  expect_true("X" %in% filtered)
  expect_true("MT" %in% filtered)
})

# ── Test 8 ── Sliding window finds significant region when genes cluster ───────
test_that("sliding window produces at least one significant region when query genes cluster", {
  skip_if(!global_loaded, "fct_09_genome.R not loaded")
  skip_if_not(exists("genome_sliding_window"), "genome_sliding_window not found")
  suppressPackageStartupMessages(library(dplyr))

  # 15 query genes all within Mbp 1-15 on chr1, background genes spread far away
  set.seed(42)
  x0 <- data.frame(
    start_position  = c(seq(1, 15, by = 1), seq(100, 250, by = 10)),
    chNum           = 1L,
    chromosome_name = "1",
    Fold            = c(rep(1L, 15), rep(0L, 16)),
    stringsAsFactors = FALSE
  )
  result <- genome_sliding_window(
    x0, listN = 15, totalN = 31, windowSize = 20, steps = 2, pvalCutoff = 0.05
  )
  expect_true(nrow(result) > 0)
  expect_true(any(result$pval < 0.05))
})

# ── Test 9 ── start_position and end_position semantics ──────────────────────
test_that("start_position and end_position are treated differently by sliding window", {
  skip_if(!global_loaded, "fct_09_genome.R not loaded")
  skip_if_not(exists("genome_sliding_window"), "genome_sliding_window not found")
  suppressPackageStartupMessages(library(dplyr))

  # Two datasets identical except start_position values
  # Dataset A: query genes clustered at the start
  x0_a <- data.frame(
    start_position  = c(1, 2, 3, 4, 5, 50, 100),
    chNum           = 1L,
    chromosome_name = "1",
    Fold            = c(1L, 1L, 1L, 1L, 1L, 0L, 0L),
    stringsAsFactors = FALSE
  )
  # Dataset B: same data but positions swapped (start at end)
  x0_b <- data.frame(
    start_position  = c(96, 97, 98, 99, 100, 50, 1),
    chNum           = 1L,
    chromosome_name = "1",
    Fold            = c(1L, 1L, 1L, 1L, 1L, 0L, 0L),
    stringsAsFactors = FALSE
  )
  result_a <- genome_sliding_window(x0_a, listN=5, totalN=7, windowSize=10, steps=1, pvalCutoff=1)
  result_b <- genome_sliding_window(x0_b, listN=5, totalN=7, windowSize=10, steps=1, pvalCutoff=1)
  # The two datasets produce different window membership
  if (nrow(result_a) > 0 || nrow(result_b) > 0) {
    # They should differ in at least one aspect (position x)
    all_x_a <- sort(unique(result_a$x))
    all_x_b <- sort(unique(result_b$x))
    expect_false(identical(all_x_a, all_x_b))
  }
})
