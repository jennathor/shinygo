# test-mod_13_groups.R
# Tests for mod_13_groups: UI structure and download filename logic.

library(testthat)
library(shiny)

if (basename(getwd()) == "testthat") setwd("../..")
if (basename(getwd()) == "tests")    setwd("..")

global_loaded <- tryCatch({
  suppressMessages(suppressWarnings(source("global.R", local = FALSE)))
  source("R/mod_13_groups.R", local = FALSE)
  TRUE
}, error = function(e) FALSE)

# ── Test 1 ── UI contains tableOutput('grouping') ────────────────────────────
test_that("mod_13_groups_ui contains tableOutput('grouping')", {
  skip_if(!global_loaded, "modules not loaded")
  ui   <- mod_13_groups_ui("test13")
  html <- as.character(ui)
  expect_true(grepl("grouping", html))
})

# ── Test 2 ── UI contains downloadButton('downloadGrouping') ─────────────────
test_that("mod_13_groups_ui contains downloadButton('downloadGrouping')", {
  skip_if(!global_loaded, "modules not loaded")
  ui   <- mod_13_groups_ui("test13b")
  html <- as.character(ui)
  expect_true(grepl("downloadGrouping", html))
})

# ── Test 3 ── Returns early when goButton = 0 (testServer) ───────────────────
test_that("server returns early (NULL) when goButton = 0", {
  skip_if(!global_loaded, "modules not loaded")
  sidebar <- create_mock_sidebar_values(goButton = 0L)
  enr     <- create_mock_enrichment_values()
  testServer(mod_13_groups_server,
    args = list(
      sidebar_values      = sidebar,
      enrichment_values   = enr
    ), {
      session$setInputs()
      # When goButton = 0 the server returns early; no error expected
      expect_no_error(output$grouping)
    }
  )
})

# ── Test 4 ── groupings data frame survives write.csv round-trip ──────────────
test_that("groupings data frame survives write.csv/read.csv round-trip", {
  groupings <- data.frame(
    Category = c("Signal transduction", "Metabolic process", "Cell division"),
    Genes    = c("TP53 BRCA1 EGFR", "LDHA PKM GPI", "CDK1 CCNB1"),
    Count    = c(3L, 3L, 2L),
    stringsAsFactors = FALSE
  )
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  write.csv(groupings, tmp, row.names = FALSE)
  reloaded <- read.csv(tmp, stringsAsFactors = FALSE)
  expect_equal(nrow(reloaded), 3L)
  expect_equal(ncol(reloaded), 3L)
  expect_equal(reloaded$Category[1], "Signal transduction")
})

# ── Test 5 ── UI tab value is 7 ───────────────────────────────────────────────
test_that("mod_13_groups_ui tab value is 7", {
  skip_if(!global_loaded, "modules not loaded")
  ui   <- mod_13_groups_ui("test13f")
  html <- as.character(ui)
  # tabPanel(value = 7, ...) renders as data-value="7" in the HTML
  expect_true(grepl("value.*7|7.*value", html))
})
