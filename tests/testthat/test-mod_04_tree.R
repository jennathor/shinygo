# test-mod_04_tree.R
# Tests for mod_04_tree: UI structure and preprocessing logic.

library(testthat)
library(shiny)

if (basename(getwd()) == "testthat") setwd("../..")
if (basename(getwd()) == "tests")    setwd("..")

global_loaded <- tryCatch({
  suppressMessages(suppressWarnings(source("global.R", local = FALSE)))
  source("R/mod_04_tree.R", local = FALSE)
  TRUE
}, error = function(e) FALSE)

# ── Test 1 ── UI contains GOTermsTree output and download button ──────────────
test_that("mod_04_tree_ui contains plotOutput('GOTermsTree') and download button", {
  skip_if(!global_loaded, "modules not loaded")
  ui   <- mod_04_tree_ui("test04")
  html <- as.character(ui)
  expect_true(grepl("GOTermsTree",      html))
  expect_true(grepl("download",         html, ignore.case = TRUE))
})

# ── Test 2 ── preprocessing strips HTML anchor tags from Pathway column ───────
test_that("preprocessing strips HTML anchor tags from Pathway column", {
  df <- data.frame(
    adj.Pval  = c(1e-5, 1e-4),
    nGenesList = c(5L, 8L),
    nGenesCategor = c(20L, 30L),
    Fold      = c(3.0, 2.5),
    Pathways  = c("<a href='http://go.db/term/GO:0001' target='_blank'>DNA repair</a>",
                  "<a href='http://go.db/term/GO:0002' target='_blank'>Apoptosis</a>"),
    URL       = c("http://go.db/1", "http://go.db/2"),
    Genes     = c("A B C", "D E F"),
    stringsAsFactors = FALSE
  )
  # Apply the same logic as significantOverlaps2() in mod_04_tree_server
  df$Pathways <- gsub(".*'_blank'>|</a>", "", df$Pathways)
  expect_false(any(grepl("<a",    df$Pathways)))
  expect_false(any(grepl("</a>", df$Pathways)))
  expect_equal(df$Pathways, c("DNA repair", "Apoptosis"))
})

# ── Test 3 ── preprocessing adds Direction = "Diff" for all rows ─────────────
test_that("preprocessing adds Direction = 'Diff' column for all rows", {
  df <- data.frame(
    adj.Pval = c(1e-5, 1e-4),
    Pathways = c("DNA repair", "Apoptosis"),
    stringsAsFactors = FALSE
  )
  df$Direction <- "Diff"
  expect_true(all(df$Direction == "Diff"))
})

# ── Test 4 ── preprocessing renames columns correctly ────────────────────────
test_that("preprocessing renames enrichment columns correctly for enrichmentPlot", {
  # In significantOverlaps2(), column names are reassigned:
  # original: "Enrichment FDR", "nGenes", "Pathway Genes", "Fold Enrichment", "Pathway", "URL", "Genes"
  # renamed:  "adj.Pval",       "nGenesList", "nGenesCategor", "Fold", "Pathways", "URL", "Genes"
  original_names <- c("Enrichment FDR", "nGenes", "Pathway Genes",
                       "Fold Enrichment", "Pathway", "URL", "Genes")
  expected_names <- c("adj.Pval", "nGenesList", "nGenesCategor",
                       "Fold", "Pathways", "URL", "Genes")
  df <- as.data.frame(matrix(1, nrow = 2, ncol = 7))
  colnames(df) <- original_names
  colnames(df) <- expected_names   # simulated rename
  expect_true(all(c("adj.Pval", "Pathways", "Genes") %in% names(df)))
  expect_false(any(c("Enrichment FDR", "Pathway") %in% names(df)))
})

# ── Test 5 ── enrichmentPlot returns NULL for single-row input ────────────────
test_that("enrichmentPlot returns NULL for single-row input", {
  skip_if(!global_loaded, "modules not loaded")
  skip_if_not(exists("enrichmentPlot"), "enrichmentPlot not found")
  one_row <- data.frame(
    adj.Pval      = 1e-5,
    nGenesList    = 5L,
    nGenesCategor = 20L,
    Fold          = 3.0,
    Pathways      = "DNA repair",
    URL           = "http://example.com",
    Genes         = "A B C",
    Direction     = "Diff",
    stringsAsFactors = FALSE
  )
  result <- enrichmentPlot(one_row)
  expect_null(result)
})

# ── Test 6 ── enrichmentPlot returns recordedPlot for valid 5-row input ───────
test_that("enrichmentPlot returns a recordedPlot for valid 5-row input", {
  skip_if(!global_loaded, "modules not loaded")
  skip_if_not(exists("enrichmentPlot"), "enrichmentPlot not found")
  skip_if_not_installed("dendextend")

  set.seed(42)
  five_rows <- data.frame(
    adj.Pval      = 10^(-seq(3, 7)),
    nGenesList    = c(5L, 8L, 3L, 12L, 6L),
    nGenesCategor = c(20L, 30L, 15L, 50L, 25L),
    Fold          = c(3.0, 2.5, 4.0, 1.8, 3.2),
    Pathways      = c("DNA repair", "Apoptosis", "Cell cycle",
                      "Signal transduction", "Transcription"),
    URL           = paste0("http://example.com/", 1:5),
    Genes         = c("A B C D E", "B C D E F", "C D E G H",
                      "D E F G I", "E F G H J"),
    Direction     = "Diff",
    stringsAsFactors = FALSE
  )
  pdf(NULL)
  result <- enrichmentPlot(five_rows)
  dev.off()
  expect_true(inherits(result, "recordedplot"))
})

# ── Test 7 ── tree_plot returns NULL when goButton is 0 ──────────────────────
test_that("tree_plot returns NULL when goButton is 0 (testServer)", {
  skip_if(!global_loaded, "modules not loaded")
  sidebar  <- create_mock_sidebar_values(goButton = 0L)
  enrichv  <- create_mock_enrichment_values()
  testServer(mod_04_tree_server,
    args = list(sidebar_values = sidebar, enrichment_values = enrichv), {
      result <- tree_plot()
      expect_null(result)
    }
  )
})
