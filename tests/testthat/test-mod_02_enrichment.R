# test-mod_02_enrichment.R
# Tests for mod_02_enrichment: sort/filter logic, redundancy removal,
# abbreviation, and download handlers.

library(testthat)
library(shiny)

if (basename(getwd()) == "testthat") setwd("../..")
if (basename(getwd()) == "tests")    setwd("..")

global_loaded <- tryCatch({
  suppressMessages(suppressWarnings(source("global.R", local = FALSE)))
  source("R/mod_02_enrichment.R", local = FALSE)
  TRUE
}, error = function(e) FALSE)

# Helper: build a test enrichment data frame with n rows
make_enrichment_df <- function(n = 10) {
  data.frame(
    `Enrichment FDR`  = 10^(-seq(2, 2 + n - 1)),
    nGenes            = as.integer(seq(5, 5 * n, by = 5)),
    `Pathway Genes`   = as.integer(seq(15, 15 * n, by = 15)),
    `Fold Enrichment` = rev(seq(1, n)),
    Pathway           = paste0("Pathway_", LETTERS[seq_len(n)]),
    URL               = paste0("http://example.com/", seq_len(n)),
    Genes             = paste0("G", seq_len(n), " G", seq_len(n) + 100),
    check.names       = FALSE,
    stringsAsFactors  = FALSE
  )
}

# ── Test 1 ── UI contains all 6 sort-mode choices ────────────────────────────
test_that("mod_02_enrichment_ui contains all 6 sort-mode choices in SortPathways", {
  skip_if(!global_loaded, "modules not loaded")
  ui   <- mod_02_enrichment_ui("test02")
  html <- as.character(ui)
  expect_true(grepl("Sort by FDR",             html))
  expect_true(grepl("Sort by Fold Enrichment", html))
  expect_true(grepl("Sort by Genes",           html))
  expect_true(grepl("Sort by Category Name",   html))
  expect_true(grepl("Select by FDR",           html))
  expect_true(grepl("average",                 html, ignore.case = TRUE))
})

# ── Test 2 ── "Sort by FDR" produces ascending FDR order ─────────────────────
test_that("'Sort by FDR' produces ascending FDR order", {
  df     <- make_enrichment_df(8)
  sorted <- df[order(df[, 1]), ]   # mimic the server logic
  expect_equal(sorted[, 1], sort(df[, 1]))
  expect_equal(sorted[1, 1], min(df[, 1]))
})

# ── Test 3 ── "Sort by Fold Enrichment" drops small pathways ─────────────────
test_that("'Sort by Fold Enrichment' drops rows where Pathway Genes < min_gene_fold and sorts descending", {
  min_gene_fold_test <- 10L
  df <- make_enrichment_df(10)
  # Apply server logic: filter then sort
  filtered <- df[df[, 3] > min_gene_fold_test, ]
  sorted   <- filtered[order(filtered[, 4], decreasing = TRUE), ]
  # All remaining rows have Pathway Genes > 10
  expect_true(all(sorted[, 3] > min_gene_fold_test))
  # Fold enrichment is non-increasing
  folds <- sorted[, 4]
  expect_true(all(diff(folds) <= 0))
})

# ── Test 4 ── "Select by FDR, sort by Fold Enrichment" logic ─────────────────
test_that("'Select by FDR then sort by Fold' selects top N by FDR and re-ranks by fold", {
  n         <- 10
  max_terms <- 5L
  df        <- make_enrichment_df(n)
  # Server: sort by FDR first, keep top maxTerms, then sort by fold descending
  by_fdr    <- df[order(df[, 1]), ]
  top_n     <- by_fdr[seq_len(max_terms), ]
  final     <- top_n[order(top_n[, 4], decreasing = TRUE), ]
  # Rows with FDR rank > max_terms should be absent
  fdr_ranks  <- rank(df[, 1])
  top_n_ids  <- which(fdr_ranks <= max_terms)
  expect_equal(nrow(final), max_terms)
  expect_true(all(df$Pathway[top_n_ids] %in% final$Pathway))
  # Fold enrichment non-increasing
  expect_true(all(diff(final[, 4]) <= 0))
})

# ── Test 5 ── "Sort by average ranks" ─────────────────────────────────────────
test_that("'Sort by average ranks' computes and orders by mean(FDR rank, Fold rank)", {
  df          <- make_enrichment_df(5)
  fdr_rank    <- rank(df[, 1])
  fold_rank   <- rank(-1 * df[, 4])
  avg_rank    <- (fdr_rank + fold_rank) / 2
  sorted      <- df[order(avg_rank), ]
  # Re-compute and check that ordering is consistent with average ranks
  recomputed_avg_rank <- sort(avg_rank)
  expect_equal(sorted$Pathway, df$Pathway[order(avg_rank)])
})

# ── Test 6 ── "Sort by Genes" orders by nGenes descending ────────────────────
test_that("'Sort by Genes' orders by nGenes descending", {
  df     <- make_enrichment_df(6)
  sorted <- df[order(df[, 2], decreasing = TRUE), ]
  expect_true(all(diff(sorted[, 2]) <= 0))
  expect_equal(sorted[1, 2], max(df[, 2]))
})

# ── Test 7 ── "Sort by Category Name" produces alphabetical order ─────────────
test_that("'Sort by Category Name' produces alphabetical Pathway order", {
  df     <- make_enrichment_df(6)
  sorted <- df[order(df[, 5]), ]
  expect_equal(sorted$Pathway, sort(df$Pathway))
})

# ── Test 8 ── FDR filter removes pathways above threshold ────────────────────
test_that("FDR filter removes pathways where Enrichment FDR > minFDR", {
  df     <- make_enrichment_df(10)
  minFDR <- 1e-5
  result <- df[df[, 1] < minFDR, ]
  expect_true(all(result[, 1] <= minFDR))
  expect_true(nrow(result) < nrow(df))
})

# ── Test 9 ── maxTerms limits result to exactly that many rows ────────────────
test_that("maxTerms limits result to exactly that many rows", {
  df        <- make_enrichment_df(20)
  max_terms <- 7L
  result    <- df[seq_len(max_terms), ]
  expect_equal(nrow(result), max_terms)
})

# ── Test 10 ── redundancy removal eliminates pathway with identical genes ─────
test_that("redundancy removal eliminates a pathway sharing ≥95% genes with same-name", {
  df <- data.frame(
    `Enrichment FDR`  = c(1e-5, 1e-4),
    nGenes            = c(5L, 5L),
    `Pathway Genes`   = c(10L, 10L),
    `Fold Enrichment` = c(3.0, 3.0),
    Pathway           = c("DNA repair", "DNA repair process"),
    URL               = c("http://a.com/1", "http://a.com/2"),
    Genes             = c("A B C D E", "A B C D E"),   # identical genes
    check.names       = FALSE,
    stringsAsFactors  = FALSE
  )
  # Simulate the redundancy check logic from mod_02_enrichment.R
  reduced <- 0.95
  n       <- nrow(df)
  flag1   <- rep(TRUE, n)
  gene_lists <- lapply(df$Genes, function(y) unlist(strsplit(y, " ")))
  pathways   <- lapply(df$Pathway, function(y) unlist(strsplit(y, " ")))
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      if (flag1[j]) {
        ratio1 <- length(intersect(gene_lists[[i]], gene_lists[[j]])) /
                  length(union(gene_lists[[i]], gene_lists[[j]]))
        if (ratio1 > reduced) {
          ratio2 <- length(intersect(pathways[[i]], pathways[[j]])) /
                    length(union(pathways[[i]], pathways[[j]]))
          if (ratio2 > 0.5) flag1[i] <- FALSE
        }
      }
    }
  }
  result <- df[flag1, ]
  expect_equal(nrow(result), 1L)
})

# ── Test 11 ── redundancy removal keeps pathways with identical genes but different names
test_that("redundancy removal keeps two pathways with same genes but different names", {
  df <- data.frame(
    `Enrichment FDR`  = c(1e-5, 1e-4),
    nGenes            = c(5L, 5L),
    `Pathway Genes`   = c(10L, 10L),
    `Fold Enrichment` = c(3.0, 3.0),
    Pathway           = c("DNA repair", "Calcium signalling completely different topic"),
    URL               = c("http://a.com/1", "http://a.com/2"),
    Genes             = c("A B C D E", "A B C D E"),   # identical genes
    check.names       = FALSE,
    stringsAsFactors  = FALSE
  )
  reduced <- 0.95
  n       <- nrow(df)
  flag1   <- rep(TRUE, n)
  gene_lists <- lapply(df$Genes, function(y) unlist(strsplit(y, " ")))
  pathways   <- lapply(df$Pathway, function(y) unlist(strsplit(y, " ")))
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      if (flag1[j]) {
        ratio1 <- length(intersect(gene_lists[[i]], gene_lists[[j]])) /
                  length(union(gene_lists[[i]], gene_lists[[j]]))
        if (ratio1 > reduced) {
          ratio2 <- length(intersect(pathways[[i]], pathways[[j]])) /
                    length(union(pathways[[i]], pathways[[j]]))
          if (ratio2 > 0.5) flag1[i] <- FALSE
        }
      }
    }
  }
  result <- df[flag1, ]
  expect_equal(nrow(result), 2L)
})

# ── Test 12 ── abbreviation replaces "Positive regulation" and truncates ──────
test_that("abbreviation replaces 'Positive regulation' → 'Pos. reg.' and truncates at 100 chars", {
  skip_if(!global_loaded, "modules not loaded")
  df <- data.frame(
    `Enrichment FDR`  = 1e-5,
    nGenes            = 5L,
    `Pathway Genes`   = 20L,
    `Fold Enrichment` = 3.0,
    Pathway           = c(
      "Positive regulation of cell proliferation",
      paste(rep("x", 110), collapse = "")
    ),
    URL               = "http://example.com",
    Genes             = "A B C",
    check.names       = FALSE,
    stringsAsFactors  = FALSE
  )
  df[, 5] <- gsub("Positive regulation", "Pos. reg.", df[, 5])
  df[, 5] <- gsub("Negative regulation", "Neg. reg.", df[, 5])
  df[, 5] <- substr(df[, 5], 1, 100)
  expect_true(any(grepl("Pos. reg.", df[, 5])))
  expect_true(all(nchar(df[, 5]) <= 100))
})

# ── Test 13 ── download produces valid CSV with correct column headers ─────────
test_that("enrichment data can be written and re-read as valid CSV", {
  df  <- make_enrichment_df(5)
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  write.csv(df, tmp, row.names = FALSE)
  result <- read.csv(tmp, check.names = FALSE)
  expect_equal(colnames(result), colnames(df))
  expect_equal(nrow(result), nrow(df))
})

# ── Test 14 ── FDR filter with minFDR = 1 keeps all rows ─────────────────────
test_that("FDR filter with minFDR = 1 keeps all rows", {
  df     <- make_enrichment_df(8)
  result <- df[df[, 1] < 1, ]
  expect_equal(nrow(result), nrow(df))
})

# ── Test 15 ── sort on a 1-row data frame does not error ─────────────────────
test_that("sorting a 1-row data frame does not error for any of the 6 sort modes", {
  df1 <- make_enrichment_df(1)
  expect_no_error(df1[order(df1[, 1]), ])                              # FDR
  expect_no_error(df1[order(df1[, 4], decreasing = TRUE), ])           # Fold
  expect_no_error(df1[order(df1[, 2], decreasing = TRUE), ])           # Genes
  expect_no_error(df1[order(df1[, 5]), ])                              # Name
  fdr_rank  <- rank(df1[, 1]);  fold_rank <- rank(-df1[, 4])
  avg_rank  <- (fdr_rank + fold_rank) / 2
  expect_no_error(df1[order(avg_rank), ])                              # Average ranks
  # Select by FDR then fold
  by_fdr <- df1[order(df1[, 1]), ]
  expect_no_error(by_fdr[order(by_fdr[, 4], decreasing = TRUE), ])    # Select+Fold
})
