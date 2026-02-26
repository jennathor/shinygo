#!/usr/bin/env Rscript
# Convenient test runner for ShinyGO
# Usage: Rscript run_tests.R [test_file]
#
# Examples:
#   Rscript run_tests.R                           # Run all tests
#   Rscript run_tests.R test-mod_01_sidebar.R     # Run specific test file

library(testthat)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

cat("╔══════════════════════════════════════╗\n")
cat("║     ShinyGO Test Suite Runner        ║\n")
cat("╚══════════════════════════════════════╝\n\n")

# Check if testthat is installed
if (!requireNamespace("testthat", quietly = TRUE)) {
  cat("❌ testthat package not installed.\n")
  cat("Installing testthat...\n")
  install.packages("testthat", repos = "https://cloud.r-project.org")
}

# Ensure we're in the app root directory
if (!file.exists("global.R") || !file.exists("ui.R") || !file.exists("server.R")) {
  cat("❌ Error: Must run from ShinyGO root directory\n")
  quit(status = 1)
}

# Check if test directory exists
if (!dir.exists("tests/testthat")) {
  cat("❌ Error: tests/testthat directory not found\n")
  quit(status = 1)
}

# Run tests
if (length(args) == 0) {
  cat("Running all tests...\n\n")

  results <- test_dir(
    "tests/testthat",
    reporter = "progress",
    stop_on_failure = FALSE
  )

  cat("\n")
  cat("══════════════════════════════════════\n")
  cat("Test Summary:\n")
  cat("══════════════════════════════════════\n")

} else {
  test_file <- args[1]
  test_path <- file.path("tests/testthat", test_file)

  if (!file.exists(test_path)) {
    cat("❌ Error: Test file not found:", test_path, "\n")
    quit(status = 1)
  }

  cat("Running tests from:", test_file, "\n\n")

  results <- test_file(
    test_path,
    reporter = "progress",
    stop_on_failure = FALSE
  )
}

cat("\n✅ Testing complete!\n")
