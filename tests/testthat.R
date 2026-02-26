# Main test runner for ShinyGO tests
# Run this file to execute all tests

library(testthat)

# Set working directory to app root
if (basename(getwd()) == "tests") {
  setwd("..")
}

cat("Starting ShinyGO test suite...\n")
cat("Working directory:", getwd(), "\n\n")

# Run all tests in testthat directory
test_dir(
  "tests/testthat",
  reporter = "progress",
  stop_on_failure = FALSE
)

cat("\n✅ Test suite completed!\n")
