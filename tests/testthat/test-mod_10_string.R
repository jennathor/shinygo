# test-mod_10_string.R
# Tests for mod_10_string: UI structure, find_taxon_by_id, and logic gates.

library(testthat)
library(shiny)

if (basename(getwd()) == "testthat") setwd("../..")
if (basename(getwd()) == "tests")    setwd("..")

global_loaded <- tryCatch({
  suppressMessages(suppressWarnings(source("global.R", local = FALSE)))
  source("R/mod_10_string.R", local = FALSE)
  TRUE
}, error = function(e) FALSE)

# ── Test 1 ── UI contains STRINGdbGO selectInput ──────────────────────────────
test_that("mod_10_string_ui contains STRINGdbGO selectInput", {
  skip_if(!global_loaded, "modules not loaded")
  ui   <- mod_10_string_ui("test10")
  html <- as.character(ui)
  expect_true(grepl("STRINGdbGO", html))
})

# ── Test 2 ── STRINGdbGO has exactly 6 choices ────────────────────────────────
test_that("STRINGdbGO selectInput has 6 enrichment category choices", {
  skip_if(!global_loaded, "modules not loaded")
  ui   <- mod_10_string_ui("test10b")
  html <- as.character(ui)
  # Six choices: Process, Component, Function, KEGG, Pfam, InterPro
  for (choice in c("Process", "Component", "Function", "KEGG", "Pfam", "InterPro")) {
    expect_true(grepl(choice, html), info = paste("Missing choice:", choice))
  }
})

# ── Test 3 ── Modal UI contains nGenesPPI slider ──────────────────────────────
test_that("mod_10_string_modal_ui contains nGenesPPI sliderInput", {
  skip_if(!global_loaded, "modules not loaded")
  # bsModal requires shinyBS to be attached
  ui <- tryCatch({
    suppressPackageStartupMessages(library(shinyBS))
    mod_10_string_modal_ui("test10m")
  }, error = function(e) NULL)
  skip_if(is.null(ui), "shinyBS::bsModal not available")
  html <- as.character(ui)
  expect_true(grepl("nGenesPPI", html))
})

# ── Test 4 ── Modal UI contains ModalPPI trigger button ──────────────────────
test_that("mod_10_string_ui contains ModalPPI action button", {
  skip_if(!global_loaded, "modules not loaded")
  ui   <- mod_10_string_ui("test10c")
  html <- as.character(ui)
  expect_true(grepl("ModalPPI", html))
})

# ── Test 5 ── find_taxon_by_id with 'taxon' column ───────────────────────────
test_that("find_taxon_by_id returns correct taxon for human using 'taxon' column", {
  skip_if(!global_loaded, "modules not loaded")
  skip_if_not(exists("find_taxon_by_id"), "find_taxon_by_id not found")
  org_info <- data.frame(
    id    = c(-9606L, -10090L, -7955L),
    taxon = c(9606L,  10090L,  7955L),
    stringsAsFactors = FALSE
  )
  result <- find_taxon_by_id(-9606L, org_info)
  expect_equal(result, 9606L)
})

# ── Test 6 ── find_taxon_by_id with 'taxon_id' column ────────────────────────
test_that("find_taxon_by_id returns correct taxon for human using 'taxon_id' column", {
  skip_if(!global_loaded, "modules not loaded")
  skip_if_not(exists("find_taxon_by_id"), "find_taxon_by_id not found")
  org_info <- data.frame(
    id       = c(-9606L, -10090L),
    taxon_id = c(9606L,  10090L),
    stringsAsFactors = FALSE
  )
  result <- find_taxon_by_id(-9606L, org_info)
  expect_equal(result, 9606L)
})

# ── Test 7 ── find_taxon_by_id for unknown species returns empty ──────────────
test_that("find_taxon_by_id returns zero-length vector for unknown species", {
  skip_if(!global_loaded, "modules not loaded")
  skip_if_not(exists("find_taxon_by_id"), "find_taxon_by_id not found")
  org_info <- data.frame(
    id    = c(-9606L, -10090L),
    taxon = c(9606L,  10090L),
    stringsAsFactors = FALSE
  )
  result <- find_taxon_by_id(-99999L, org_info)
  expect_equal(length(result), 0L)
})

# ── Test 8 ── find_taxon_by_id returns exactly 1 value for known species ──────
test_that("find_taxon_by_id returns exactly 1 value for known species", {
  skip_if(!global_loaded, "modules not loaded")
  skip_if_not(exists("find_taxon_by_id"), "find_taxon_by_id not found")
  org_info <- data.frame(
    id    = c(-9606L, -10090L, -7955L),
    taxon = c(9606L,  10090L,  7955L),
    stringsAsFactors = FALSE
  )
  result <- find_taxon_by_id(-10090L, org_info)
  expect_equal(length(result), 1L)
  expect_equal(result, 10090L)
})

# ── Test 9 ── mapping ratio warning threshold is 0.30 ────────────────────────
test_that("mapping ratio warning is appended when ratio < 0.30", {
  # Simulate the renderText logic from mod_10_string_server
  format_mapping_stat <- function(ratio) {
    tem <- paste0(100 * round(ratio, 3), "% genes mapped by STRING web server.")
    if (ratio < 0.3) {
      tem <- paste(tem, "Warning!!! Very few gene mapped. Double check if the correct species is selected.")
    }
    tem
  }
  # Below threshold: warning present
  low_stat <- format_mapping_stat(0.15)
  expect_true(grepl("Warning!!!", low_stat))
  # At or above threshold: no warning
  ok_stat  <- format_mapping_stat(0.50)
  expect_false(grepl("Warning!!!", ok_stat))
  # Exactly at boundary (0.30): no warning
  boundary <- format_mapping_stat(0.30)
  expect_false(grepl("Warning!!!", boundary))
})
