# test-mod_11_about.R
# Tests for mod_11_about: static UI content.

library(testthat)
library(shiny)

if (basename(getwd()) == "testthat") setwd("../..")
if (basename(getwd()) == "tests")    setwd("..")

global_loaded <- tryCatch({
  suppressMessages(suppressWarnings(source("global.R", local = FALSE)))
  source("R/mod_11_about.R", local = FALSE)
  TRUE
}, error = function(e) FALSE)

# ── Test 1 ── UI contains citation Bioinformatics 36:2628–2629, 2020 ─────────
test_that("mod_11_about_ui contains the Bioinformatics citation", {
  skip_if(!global_loaded, "modules not loaded")
  ui   <- mod_11_about_ui("test11")
  html <- as.character(ui)
  # Ge SX, Jung D & Yao R, Bioinformatics 36:2628–2629, 2020
  expect_true(grepl("Bioinformatics", html))
  expect_true(grepl("2020", html))
  expect_true(grepl("Ge SX", html))
})

# ── Test 2 ── UI contains GitHub repository link ──────────────────────────────
test_that("mod_11_about_ui contains GitHub repository link", {
  skip_if(!global_loaded, "modules not loaded")
  ui   <- mod_11_about_ui("test11b")
  html <- as.character(ui)
  expect_true(grepl("github.com", html, ignore.case = TRUE))
  expect_true(grepl("shinygo", html, ignore.case = TRUE))
})
