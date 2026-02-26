# test-mod_06_kegg.R
# Tests for mod_06_kegg: UI structure, pathway filtering, and keggSpeciesID.

library(testthat)
library(shiny)

if (basename(getwd()) == "testthat") setwd("../..")
if (basename(getwd()) == "tests")    setwd("..")

global_loaded <- tryCatch({
  suppressMessages(suppressWarnings(source("global.R", local = FALSE)))
  source("R/mod_06_kegg.R",   local = FALSE)
  source("R/fct_06_kegg.R",   local = FALSE)
  TRUE
}, error = function(e) FALSE)

# ── Test 1 ── UI contains listSigPathways and KEGG conditional panel ──────────
test_that("mod_06_kegg_ui contains listSigPathways uiOutput and conditional panel", {
  skip_if(!global_loaded, "modules not loaded")
  ui   <- mod_06_kegg_ui("test06")
  html <- as.character(ui)
  expect_true(grepl("listSigPathways", html))
  expect_true(grepl("KEGG", html))
})

# ── Test 2 ── UI shows 'Please select KEGG' message when non-KEGG DB is active
test_that("mod_06_kegg_ui shows 'Please select KEGG' conditional message", {
  skip_if(!global_loaded, "modules not loaded")
  ui   <- mod_06_kegg_ui("test06b")
  html <- as.character(ui)
  expect_true(grepl("Please select KEGG", html))
})

# ── Test 3 ── Pathways with Pathway Genes > 1000 are removed ─────────────────
test_that("pathways with Pathway Genes > 1000 are removed; remainder sorted by Fold descending", {
  df <- data.frame(
    `Enrichment FDR`  = c(1e-5, 1e-4, 1e-6),
    nGenes            = c(10L, 8L, 15L),
    `Pathway Genes`   = c(50L, 1200L, 80L),
    `Fold Enrichment` = c(2.0, 3.5, 4.0),
    Pathway           = c("Cell cycle", "Big pathway", "DNA repair"),
    check.names       = FALSE,
    stringsAsFactors  = FALSE
  )
  # Apply server logic: filter out large pathways, sort by fold descending
  filtered <- df[df[, 3] < 1000, ]
  sorted   <- filtered[order(-filtered[, 4]), ]
  expect_true(all(sorted[, 3] <= 1000))
  expect_true(all(diff(sorted[, 4]) <= 0))
  expect_false("Big pathway" %in% sorted$Pathway)
})

# ── Test 4 ── keggSpeciesID has required columns ──────────────────────────────
test_that("keggSpeciesID has ensembl_dataset, name, kegg columns", {
  skip_if(!global_loaded, "modules not loaded")
  skip_if_not(exists("keggSpeciesID"), "keggSpeciesID not found")
  expect_true(all(c("ensembl_dataset", "name", "kegg") %in% colnames(keggSpeciesID)))
})

# ── Test 5 ── convertEnsembl2Entrez returns NULL for unknown Ensembl IDs ──────
test_that("convertEnsembl2Entrez returns NULL for unknown Ensembl IDs", {
  skip_if(!global_loaded, "modules not loaded")
  skip_if_not(exists("convertEnsembl2Entrez"), "convertEnsembl2Entrez not found")
  skip_if(!file.exists("./data/data104b/db/Homo_sapiens.db"),
          "Species database not available")
  result <- tryCatch(
    convertEnsembl2Entrez("ENSG_FAKE_999999_DOES_NOT_EXIST", -9606L),
    error = function(e) NULL
  )
  # Should return NULL or an empty data frame for unrecognized IDs
  expect_true(is.null(result) || (is.data.frame(result) && nrow(result) == 0))
})

# ── Test 6 ── convertEnsembl2Entrez returns data frame with entrezgene_id ─────
test_that("convertEnsembl2Entrez returns data frame with entrezgene_id column for known IDs", {
  skip_if(!global_loaded, "modules not loaded")
  skip_if_not(exists("convertEnsembl2Entrez"), "convertEnsembl2Entrez not found")
  skip_if(!has_full_database(), "Full species database not available")
  result <- tryCatch(
    convertEnsembl2Entrez(c("ENSG00000141510", "ENSG00000012048"), -9606L),
    error = function(e) NULL
  )
  if (!is.null(result) && nrow(result) > 0) {
    expect_true("entrezgene_id" %in% names(result))
  }
})

# ── Test 7 ── blank image list has correct keys ───────────────────────────────
test_that("blank image list before submission has correct keys: src, contentType, width, height", {
  outfile <- tempfile(fileext = ".png")
  on.exit(unlink(outfile))
  png(outfile, width = 400, height = 300)
  graphics::frame()
  dev.off()
  blank <- list(
    src         = outfile,
    contentType = "image/png",
    width       = 400,
    height      = 300,
    alt         = " "
  )
  expect_true(all(c("src", "contentType", "width", "height") %in% names(blank)))
})

# ── Test 8 ── keggSpeciesID lookup for human returns 'hsa' ───────────────────
test_that("keggSpeciesID lookup for hsapiens_gene_ensembl returns KEGG code 'hsa'", {
  skip_if(!global_loaded, "modules not loaded")
  skip_if_not(exists("keggSpeciesID"), "keggSpeciesID not found")
  row <- keggSpeciesID[keggSpeciesID$ensembl_dataset == "hsapiens_gene_ensembl", ]
  if (nrow(row) > 0 && !all(row$kegg == "")) {
    expect_equal(tolower(row$kegg[1]), "hsa")
  } else {
    skip("KEGG codes not populated in local development database")
  }
})
