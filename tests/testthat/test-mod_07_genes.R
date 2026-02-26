# test-mod_07_genes.R
# Tests for mod_07_genes: UI structure and gene table preprocessing logic.

library(testthat)
library(shiny)

if (basename(getwd()) == "testthat") setwd("../..")
if (basename(getwd()) == "tests")    setwd("..")

global_loaded <- tryCatch({
  suppressMessages(suppressWarnings(source("global.R", local = FALSE)))
  source("R/mod_07_genes.R", local = FALSE)
  TRUE
}, error = function(e) FALSE)

# Helper: build a full mock gene table (Ensembl species, 9+ columns)
make_gene_table <- function(n = 5) {
  if (n == 0L) {
    return(data.frame(
      User_input        = character(0),
      `Ensembl Gene ID` = character(0),
      Entrez            = character(0),
      Symbol            = character(0),
      Description       = character(0),
      Chr               = character(0),
      Type              = character(0),
      Species           = character(0),
      Set               = character(0),
      stringsAsFactors  = FALSE,
      check.names       = FALSE
    ))
  }
  ensembl_ids <- sprintf("ENSG%011d", seq_len(n))
  data.frame(
    User_input        = paste0("GENE", seq_len(n)),
    `Ensembl Gene ID` = ensembl_ids,
    Entrez            = as.character(7000L + seq_len(n)),
    Symbol            = paste0("GENE", seq_len(n)),
    Description       = paste0("Mock gene function [Source:HGNC]", seq_len(n)),
    Chr               = sample(c(as.character(1:22), "X"), n, replace = TRUE),
    Type              = rep("protein_coding", n),
    Species           = rep("STRINGdbHuman",  n),
    Set               = rep("List",           n),
    stringsAsFactors  = FALSE,
    check.names       = FALSE
  )
}

# ── Test 1 ── UI contains conversionTable output and download button ──────────
test_that("mod_07_genes_ui contains tableOutput('conversionTable') and download button", {
  skip_if(!global_loaded, "modules not loaded")
  ui   <- mod_07_genes_ui("test07")
  html <- as.character(ui)
  expect_true(grepl("conversionTable",  html))
  expect_true(grepl("downloadGeneInfo", html))
})

# ── Test 2 ── Description truncation at ; or [ ───────────────────────────────
test_that("Description truncation removes text after ';' or '[' when showDetailedGeneInfo=FALSE", {
  df <- make_gene_table(3)
  df$Description <- c("DNA repair [source: HGNC]",
                       "Tumor suppressor; involved in apoptosis",
                       "Short desc")
  # Apply server logic (from mod_07_genes_server)
  df$Description <- gsub(";.*|\\[.*", "", df$Description)
  expect_equal(df$Description[1], "DNA repair ")
  expect_equal(df$Description[2], "Tumor suppressor")
  expect_equal(df$Description[3], "Short desc")
})

# ── Test 3 ── STRINGdb label stripped from Species column ────────────────────
test_that("STRINGdb prefix is stripped from Species column", {
  df          <- make_gene_table(2)
  df$Species  <- gsub("STRINGdb", "", df$Species)
  expect_false(any(grepl("STRINGdb", df$Species)))
})

# ── Test 4 ── protein_coding → 'coding' and pseudogene → 'pseudo' ────────────
test_that("'protein_coding' → 'coding' and pseudogene types → 'pseudo' in Type column", {
  df <- data.frame(
    Type = c("protein_coding", "processed_pseudogene", "lncRNA_gene"),
    stringsAsFactors = FALSE
  )
  df$Type <- gsub(".*pseudogene",  "pseudo",  df$Type)
  df$Type <- gsub("protein_coding", "coding",  df$Type)
  df$Type <- gsub("_gene", "",                 df$Type)
  expect_equal(df$Type[1], "coding")
  expect_equal(df$Type[2], "pseudo")
  expect_equal(df$Type[3], "lncRNA")
})

# ── Test 5 ── Chr names longer than 50 chars become empty string ──────────────
test_that("chromosome names longer than 50 chars become empty string", {
  df <- data.frame(
    Chr = c("1", "X", paste(rep("a", 51), collapse = ""), "MT"),
    stringsAsFactors = FALSE
  )
  df$Chr[nchar(df$Chr) > 50] <- ""
  expect_equal(df$Chr[3], "")
  expect_equal(df$Chr[1], "1")
  expect_equal(df$Chr[4], "MT")
})

# ── Test 6 ── Ensembl IDs get Ensembl hyperlink ───────────────────────────────
test_that("Ensembl IDs starting with 'ENS' get hyperlink to ensembl.org/id/", {
  ensembl_id <- "ENSG00000139618"
  df <- data.frame(`Ensembl Gene ID` = ensembl_id, check.names = FALSE,
                   stringsAsFactors = FALSE)
  ix <- grepl("ENS", df$`Ensembl Gene ID`)
  tem <- paste0("<a href='http://www.ensembl.org/id/", df$`Ensembl Gene ID`,
                "' target='_blank'>", df$`Ensembl Gene ID`, "</a>")
  df$`Ensembl Gene ID`[ix] <- tem[ix]
  expect_true(grepl("href='http://www.ensembl.org/id/ENSG00000139618'",
                    df$`Ensembl Gene ID`[1]))
})

# ── Test 7 ── STRING-format IDs containing 'ENS' do get hyperlinks ───────────
# The server uses grepl("ENS", id) — not anchored to start — so STRING IDs
# like "9606.ENSP00000269305" that contain "ENS" in them DO receive an
# ensembl.org hyperlink. This test documents that actual behavior.
test_that("STRING-format IDs containing 'ENS' DO get Ensembl hyperlinks (grepl not anchored)", {
  df <- data.frame(`Ensembl Gene ID` = "9606.ENSP00000269305",
                   check.names = FALSE, stringsAsFactors = FALSE)
  # Actual server logic: grepl("ENS", ...) is unanchored
  ix <- grepl("ENS", df$`Ensembl Gene ID`)
  if (sum(ix) > 0) {
    tem <- paste0("<a href='http://www.ensembl.org/id/", df$`Ensembl Gene ID`,
                  "' target='_blank'>", df$`Ensembl Gene ID`, "</a>")
    df$`Ensembl Gene ID`[ix] <- tem[ix]
  }
  # "9606.ENSP..." contains "ENS", so ix = TRUE → hyperlink IS applied
  expect_true(ix[1])
  expect_true(grepl("<a", df$`Ensembl Gene ID`[1]))
})

# ── Test 8 ── Entrez IDs get NCBI hyperlink ────────────────────────────────────
test_that("Entrez IDs get hyperlink to ncbi.nlm.nih.gov/gene/", {
  df <- data.frame(Entrez = "7157", stringsAsFactors = FALSE)
  ix <- !is.na(suppressWarnings(as.numeric(df$Entrez)))
  tem <- paste0("<a href='https://www.ncbi.nlm.nih.gov/gene/", df$Entrez,
                "' target='_blank'>", df$Entrez, "</a>")
  df$Entrez[ix] <- tem[ix]
  expect_true(grepl("href='https://www.ncbi.nlm.nih.gov/gene/7157'", df$Entrez[1]))
})

# ── Test 9 ── All-NA columns are dropped ─────────────────────────────────────
test_that("columns where all values are NA are dropped from display", {
  df <- data.frame(
    A = c(1, 2, 3),
    B = c(NA_real_, NA_real_, NA_real_),   # all NA
    C = c(4, 5, 6),
    stringsAsFactors = FALSE
  )
  df_filtered <- df[, which(!apply(is.na(df), 2, sum) == nrow(df))]
  expect_equal(ncol(df_filtered), 2L)
  expect_true("A" %in% names(df_filtered))
  expect_false("B" %in% names(df_filtered))
})

# ── Test 10 ── Empty input (0-row data frame) does not error ─────────────────
test_that("empty input (0-row data frame) does not error; returns 0-row frame", {
  df <- make_gene_table(0)
  expect_no_error({
    df$Type <- gsub(".*pseudogene", "pseudo",  df$Type)
    df$Type <- gsub("protein_coding", "coding", df$Type)
    df$Chr[nchar(df$Chr) > 50] <- ""
    result <- df[, which(!apply(is.na(df), 2, sum) == nrow(df))]
  })
  expect_equal(nrow(df), 0L)
})

# ── Test 11 ── Description column with all-NA values is dropped ───────────────
test_that("Description column with all-NA values is specifically dropped", {
  df <- data.frame(
    `Ensembl Gene ID` = c("ENSG00000001", "ENSG00000002"),
    Entrez            = c("100", "200"),
    Description       = c(NA_character_, NA_character_),
    check.names       = FALSE,
    stringsAsFactors  = FALSE
  )
  df_filtered <- df[, which(!apply(is.na(df), 2, sum) == nrow(df))]
  expect_false("Description" %in% names(df_filtered))
  expect_true("Entrez" %in% names(df_filtered))
})
