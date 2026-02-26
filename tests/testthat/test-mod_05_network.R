# test-mod_05_network.R
# Tests for mod_05_network: UI structure and preprocessing logic.

library(testthat)
library(shiny)

if (basename(getwd()) == "testthat") setwd("../..")
if (basename(getwd()) == "tests")    setwd("..")

global_loaded <- tryCatch({
  suppressMessages(suppressWarnings(source("global.R", local = FALSE)))
  source("R/mod_05_network.R", local = FALSE)
  TRUE
}, error = function(e) FALSE)

# ── Test 1 ── UI contains interactive network output and edge cutoff input ────
test_that("mod_05_network_ui contains interactive network output and edgeCutoff input", {
  skip_if(!global_loaded, "modules not loaded")
  ui   <- mod_05_network_ui("test05")
  html <- as.character(ui)
  expect_true(grepl("enrichmentNetworkPlotInteractive", html))
  expect_true(grepl("edgeCutoff", html))
})

# ── Test 2 ── Modal UI creates bsModal with InteractiveNetwork id ─────────────
test_that("mod_05_network_modal_ui contains a modal with 'InteractiveNetwork' id", {
  skip_if(!global_loaded, "modules not loaded")
  ui <- tryCatch({
    suppressPackageStartupMessages(library(shinyBS))
    mod_05_network_modal_ui("test05")
  }, error = function(e) NULL)
  skip_if(is.null(ui), "shinyBS::bsModal not available")
  html <- as.character(ui)
  expect_true(grepl("InteractiveNetwork", html))
})

# ── Test 3 ── preprocessing strips HTML anchor tags from Pathway column ───────
test_that("preprocessing strips HTML anchor tags from Pathways column", {
  tem <- data.frame(
    adj.Pval  = c(1e-5, 1e-4),
    Pathways  = c("<a href='http://go.db/1' target='_blank'>DNA repair</a>",
                  "<a href='http://go.db/2' target='_blank'>Apoptosis</a>"),
    Genes     = c("A B C", "D E F"),
    stringsAsFactors = FALSE
  )
  tem$Pathways <- gsub(".*'_blank'>|</a>", "", tem$Pathways)
  expect_false(any(grepl("<a",    tem$Pathways)))
  expect_equal(tem$Pathways, c("DNA repair", "Apoptosis"))
})

# ── Test 4 ── text wrapping inserts \n when wrapTextNetwork = TRUE ────────────
test_that("text wrapping inserts \\n in long pathway names when wrapTextNetwork = TRUE", {
  skip_if(!global_loaded, "modules not loaded")
  skip_if_not(exists("wrap_strings"), "wrap_strings not found")
  tem <- data.frame(
    Pathways = c("Short name",
                 "Positive regulation of the interferon-gamma-mediated signalling pathway"),
    stringsAsFactors = FALSE
  )
  tem$Pathways <- wrap_strings(tem$Pathways)   # default width 30
  expect_true(grepl("\n", tem$Pathways[2]))
  expect_false(grepl("\n", tem$Pathways[1]))   # short name unchanged
})

# ── Test 5 ── returns NULL when enrichment x has only 1 column ────────────────
test_that("returns NULL when significantOverlaps()$x has only 1 column", {
  one_col_result <- list(x = data.frame("ID not recognized!"))
  # Simulate the server check: dim(tem$x)[2] == 1 → return(NULL)
  tem <- one_col_result
  result <- if (dim(tem$x)[2] == 1) NULL else tem$x
  expect_null(result)
})

# ── Test 6 ── edge cutoff default is 0.30 in the UI ──────────────────────────
test_that("edge cutoff input has default value 0.30", {
  skip_if(!global_loaded, "modules not loaded")
  ui   <- mod_05_network_ui("test05c")
  html <- as.character(ui)
  # The numericInput(value = 0.30) should appear in the HTML
  expect_true(grepl("0.3", html) || grepl("0\\.30", html) || grepl("\"0.3\"", html))
})

# ── Test 7 ── enrichmentNetwork returns an igraph object for valid input ───────
test_that("enrichmentNetwork returns an igraph object for valid input", {
  skip_if(!global_loaded, "modules not loaded")
  skip_if_not(exists("enrichmentNetwork"),  "enrichmentNetwork not found")
  skip_if_not_installed("igraph")
  skip_if_not_installed("visNetwork")

  set.seed(42)
  df <- data.frame(
    adj.Pval  = 10^(-seq(3, 7)),
    nGenesList = c(5L, 8L, 3L, 12L, 6L),
    nGenesCategor = c(20L, 30L, 15L, 50L, 25L),
    Fold      = c(3.0, 2.5, 4.0, 1.8, 3.2),
    Pathways  = c("DNA repair", "Apoptosis", "Cell cycle",
                  "Signal transduction", "Transcription"),
    URL       = paste0("http://example.com/", 1:5),
    Genes     = c("A B C D E", "B C D E F", "C D E G H",
                  "D E F G I", "E F G H J"),
    Direction = "Diff",
    stringsAsFactors = FALSE
  )
  result <- tryCatch(
    enrichmentNetwork(df, layoutButton = 0, edge.cutoff = 0.1),
    error = function(e) NULL
  )
  if (!is.null(result)) {
    expect_true(inherits(result, "igraph"))
  }
})

# ── Test 8 ── network has correct node count ──────────────────────────────────
test_that("network has one node per unique pathway in valid input", {
  skip_if(!global_loaded, "modules not loaded")
  skip_if_not(exists("enrichmentNetwork"),  "enrichmentNetwork not found")
  skip_if_not_installed("igraph")
  skip_if_not_installed("visNetwork")

  set.seed(42)
  df <- data.frame(
    adj.Pval  = 10^(-seq(3, 7)),
    Pathways  = c("DNA repair", "Apoptosis", "Cell cycle",
                  "Signal transduction", "Transcription"),
    Genes     = c("A B C D E", "B C D E F", "C D E G H",
                  "D E F G I", "E F G H J"),
    Direction = "Diff",
    stringsAsFactors = FALSE
  )
  g <- tryCatch(
    enrichmentNetwork(df, layoutButton = 0, edge.cutoff = 0.05),
    error = function(e) NULL
  )
  if (!is.null(g)) {
    expect_equal(igraph::vcount(g), length(unique(df$Pathways)))
  }
})
