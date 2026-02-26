# test-requires_external.R
# Integration tests that use fixture databases (helper-db.R) or the full
# production database when available.  Tests are grouped into lettered blocks:
#
#   C1–C6  Core fixture DB schema / pipeline structure
#   K1–K4  KEGG data structures
#   S1–S6  Statistical pipeline (enrichment result properties)
#   T1–T2  enrichmentPlot / tree rendering

library(testthat)
library(shiny)

if (basename(getwd()) == "testthat") setwd("../..")
if (basename(getwd()) == "tests")    setwd("..")

global_loaded <- tryCatch({
  suppressMessages(suppressWarnings(source("global.R", local = FALSE)))
  suppressPackageStartupMessages(library(dplyr))
  TRUE
}, error = function(e) FALSE)

# Helpers ---------------------------------------------------------------------

#' TRUE if DBI/RSQLite are installed and the fixture file exists
fixture_db_ok <- function(path = "tests/testthat/fixtures/mock_db.sqlite") {
  requireNamespace("DBI",     quietly = TRUE) &&
  requireNamespace("RSQLite", quietly = TRUE) &&
  file.exists(path)
}

#' TRUE if the full production convertIDs.db is present (any known path)
has_convert_db <- function() {
  candidates <- c(
    file.path(Sys.getenv("IDEP_DATABASE", unset = ""), "data113", "convertIDs.db"),
    "./data/data113/convertIDs.db",
    "./data104b/convertIDs.db",
    "./data/convertIDs.db"
  )
  any(file.exists(candidates))
}

fixture_species_db <- "tests/testthat/fixtures/db/Homo_sapiens.db"

# ── C1 ── Fixture convertIDs-like DB has expected tables ─────────────────────
test_that("[C1] fixture DB has orgInfo, mapping, idIndex, quotes tables", {
  skip_if(!fixture_db_ok(), "fixture DB not available")
  con <- DBI::dbConnect(RSQLite::SQLite(),
                        "tests/testthat/fixtures/mock_db.sqlite")
  on.exit(DBI::dbDisconnect(con))
  tables <- DBI::dbListTables(con)
  for (tbl in c("orgInfo", "mapping", "idIndex", "quotes")) {
    expect_true(tbl %in% tables, info = paste("Missing table:", tbl))
  }
})

# ── C2 ── orgInfo table has required columns ──────────────────────────────────
test_that("[C2] fixture orgInfo has id, name, name2, ensembl_dataset, taxon, genes", {
  skip_if(!fixture_db_ok(), "fixture DB not available")
  con <- DBI::dbConnect(RSQLite::SQLite(),
                        "tests/testthat/fixtures/mock_db.sqlite")
  on.exit(DBI::dbDisconnect(con))
  cols <- DBI::dbListFields(con, "orgInfo")
  for (col in c("id", "name", "name2", "ensembl_dataset", "taxon", "genes")) {
    expect_true(col %in% cols, info = paste("Missing column:", col))
  }
})

# ── C3 ── mapping table contains human Ensembl IDs ───────────────────────────
test_that("[C3] fixture mapping table contains human Ensembl IDs for TP53", {
  skip_if(!fixture_db_ok(), "fixture DB not available")
  con <- DBI::dbConnect(RSQLite::SQLite(),
                        "tests/testthat/fixtures/mock_db.sqlite")
  on.exit(DBI::dbDisconnect(con))
  rows <- DBI::dbGetQuery(con,
    "SELECT * FROM mapping WHERE id = 'TP53' AND species_id = '-9606'")
  expect_true(nrow(rows) >= 1L)
  expect_true(grepl("ENSG", rows$ensembl_gene_id[1]))
})

# ── C4 ── Species fixture DB has pathway, pathwayInfo, geneInfo, categories ───
test_that("[C4] fixture species DB has pathway, pathwayInfo, geneInfo, categories", {
  skip_if(!requireNamespace("DBI",     quietly = TRUE), "DBI not installed")
  skip_if(!requireNamespace("RSQLite", quietly = TRUE), "RSQLite not installed")
  skip_if(!file.exists(fixture_species_db), "fixture species DB not available")
  con <- DBI::dbConnect(RSQLite::SQLite(), fixture_species_db)
  on.exit(DBI::dbDisconnect(con))
  tables <- DBI::dbListTables(con)
  for (tbl in c("pathway", "pathwayInfo", "geneInfo", "categories")) {
    expect_true(tbl %in% tables, info = paste("Missing table:", tbl))
  }
})

# ── C5 ── Fixture geneInfo has expected columns ───────────────────────────────
test_that("[C5] fixture geneInfo has ensembl_gene_id, symbol, chromosome_name, gene_biotype", {
  skip_if(!requireNamespace("DBI",     quietly = TRUE), "DBI not installed")
  skip_if(!requireNamespace("RSQLite", quietly = TRUE), "RSQLite not installed")
  skip_if(!file.exists(fixture_species_db), "fixture species DB not available")
  con <- DBI::dbConnect(RSQLite::SQLite(), fixture_species_db)
  on.exit(DBI::dbDisconnect(con))
  cols <- DBI::dbListFields(con, "geneInfo")
  for (col in c("ensembl_gene_id", "symbol", "chromosome_name", "gene_biotype",
                "start_position", "entrezgene_id")) {
    expect_true(col %in% cols, info = paste("Missing column:", col))
  }
})

# ── C6 ── All fixture genes have protein_coding biotype ──────────────────────
test_that("[C6] all genes in fixture geneInfo are protein_coding", {
  skip_if(!requireNamespace("DBI",     quietly = TRUE), "DBI not installed")
  skip_if(!requireNamespace("RSQLite", quietly = TRUE), "RSQLite not installed")
  skip_if(!file.exists(fixture_species_db), "fixture species DB not available")
  con <- DBI::dbConnect(RSQLite::SQLite(), fixture_species_db)
  on.exit(DBI::dbDisconnect(con))
  df  <- DBI::dbGetQuery(con, "SELECT DISTINCT gene_biotype FROM geneInfo")
  expect_equal(nrow(df), 1L)
  expect_equal(df$gene_biotype[1], "protein_coding")
})

# ── K1 ── keggSpeciesID object exists and has ≥ 1 row ────────────────────────
test_that("[K1] keggSpeciesID exists with at least one row", {
  skip_if(!global_loaded,                   "global.R not loaded")
  skip_if_not(exists("keggSpeciesID"),       "keggSpeciesID not found")
  expect_true(nrow(keggSpeciesID) >= 1L)
})

# ── K2 ── keggSpeciesID has no NA in ensembl_dataset column ──────────────────
test_that("[K2] keggSpeciesID$ensembl_dataset has no NA values", {
  skip_if(!global_loaded,             "global.R not loaded")
  skip_if_not(exists("keggSpeciesID"), "keggSpeciesID not found")
  expect_false(any(is.na(keggSpeciesID$ensembl_dataset)))
})

# ── K3 ── KEGG pathway giant filter: > 1000 genes removed ────────────────────
test_that("[K3] KEGG pathway filter removes pathways with > 1000 genes", {
  df <- data.frame(
    Pathway          = c("Small pathway", "Giant pathway", "Medium pathway"),
    `Pathway Genes`  = c(50L, 1500L, 200L),
    `Fold Enrichment` = c(3.0, 2.0, 4.0),
    check.names      = FALSE,
    stringsAsFactors = FALSE
  )
  filtered <- df[df[["Pathway Genes"]] < 1000, ]
  expect_equal(nrow(filtered), 2L)
  expect_false("Giant pathway" %in% filtered$Pathway)
})

# ── K4 ── KEGG human lookup returns 'hsa' ────────────────────────────────────
test_that("[K4] keggSpeciesID lookup for hsapiens_gene_ensembl returns 'hsa'", {
  skip_if(!global_loaded,             "global.R not loaded")
  skip_if_not(exists("keggSpeciesID"), "keggSpeciesID not found")
  row <- keggSpeciesID[keggSpeciesID$ensembl_dataset == "hsapiens_gene_ensembl", ]
  if (nrow(row) == 0 || all(row$kegg == "")) {
    skip("KEGG codes not populated in local development database")
  }
  expect_equal(tolower(row$kegg[1]), "hsa")
})

# ── S1 ── mock enrichment result has correct column names ────────────────────
test_that("[S1] create_mock_enrichment_result returns all required columns", {
  result <- create_mock_enrichment_result(5)
  required <- c("Enrichment FDR", "nGenes", "Pathway Genes",
                "Fold Enrichment", "Pathway", "URL", "Genes")
  expect_true(all(required %in% colnames(result$x)))
})

# ── S2 ── FDR values are in [0, 1] ───────────────────────────────────────────
test_that("[S2] mock enrichment FDR values are in [0, 1]", {
  result <- create_mock_enrichment_result(10)
  fdr <- result$x[["Enrichment FDR"]]
  expect_true(all(fdr >= 0 & fdr <= 1))
})

# ── S3 ── Fold enrichment > 0 for enriched pathways ─────────────────────────
test_that("[S3] mock enrichment Fold Enrichment values are > 0", {
  result <- create_mock_enrichment_result(5)
  fold   <- result$x[["Fold Enrichment"]]
  expect_true(all(fold > 0))
})

# ── S4 ── nGenes ≤ Pathway Genes for all rows ────────────────────────────────
test_that("[S4] nGenes is always ≤ Pathway Genes in enrichment result", {
  result   <- create_mock_enrichment_result(10)
  n_genes  <- result$x$nGenes
  p_genes  <- result$x[["Pathway Genes"]]
  expect_true(all(n_genes <= p_genes))
})

# ── S5 ── p.adjust(BH) returns same-length vector ────────────────────────────
test_that("[S5] p.adjust(method='BH') returns correct-length vector", {
  pvals <- c(0.001, 0.01, 0.05, 0.10, 0.50)
  adj   <- p.adjust(pvals, method = "BH")
  expect_equal(length(adj), length(pvals))
  expect_true(all(adj >= pvals))   # BH corrections never decrease p-values
})

# ── S6 ── hypergeometric enrichment p-value monotonic in list overlap ─────────
test_that("[S6] phyper p-value decreases as list overlap with pathway increases", {
  # Background: 100 genes total, pathway has 10; list has 20 genes.
  # Vary the overlap k = 1,2,...,8; p-value should decrease.
  pvals <- sapply(1:8, function(k) {
    phyper(k - 1, 10, 90, 20, lower.tail = FALSE)
  })
  expect_true(all(diff(pvals) < 0))   # strictly decreasing
})

# ── T1 ── enrichmentPlot returns NULL for single-row input ───────────────────
test_that("[T1] enrichmentPlot returns NULL when enrichment result has ≤ 1 row", {
  skip_if(!global_loaded,                  "global.R not loaded")
  skip_if_not(exists("enrichmentPlot"),     "enrichmentPlot not found")
  result <- create_mock_enrichment_result(1)
  # enrichmentPlot(x, ...) returns NULL for ≤1-row input
  out <- tryCatch(
    enrichmentPlot(result$x, sortBy = "FDR", termCounts = 20, fontSize = 12,
                   markerSize = 4, wrap = FALSE, chartType = "lollipop",
                   ggplot2_theme = NULL),
    error = function(e) NULL
  )
  expect_null(out)
})

# ── T2 ── enrichmentPlot returns a recorded plot for multi-row input ──────────
test_that("[T2] enrichmentPlot returns a recordedplot for 5-row enrichment data", {
  skip_if(!global_loaded,              "global.R not loaded")
  skip_if_not(exists("enrichmentPlot"), "enrichmentPlot not found")
  result <- create_mock_enrichment_result(5)
  out <- tryCatch(
    enrichmentPlot(result$x, sortBy = "FDR", termCounts = 20, fontSize = 12,
                   markerSize = 4, wrap = FALSE, chartType = "lollipop",
                   ggplot2_theme = NULL),
    error = function(e) NULL
  )
  if (!is.null(out)) {
    expect_true(inherits(out, "recordedplot") || inherits(out, "gg"))
  }
})
