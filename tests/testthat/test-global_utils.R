# test-global_utils.R
# Tests for pure utility functions in global.R.
# Most tests do not require a database connection.

library(testthat)
library(shiny)

# Set working directory to project root if needed
if (basename(getwd()) == "testthat") setwd("../..")
if (basename(getwd()) == "tests")    setwd("..")

# Source global.R; skip file if it cannot be loaded (e.g., no database)
global_loaded <- tryCatch({
  suppressMessages(suppressWarnings(source("global.R", local = FALSE)))
  TRUE
}, error = function(e) {
  message("global.R could not be sourced: ", conditionMessage(e))
  FALSE
})

# ── Test 1 ── remove_gene_version: Ensembl human with version 2 ─────────────
test_that("remove_gene_version strips .N suffix from human Ensembl IDs", {
  skip_if(!global_loaded, "global.R not loaded")
  skip_if_not(exists("remove_gene_version"), "remove_gene_version not found")
  expect_equal(remove_gene_version("ENSG00000139618.2"),  "ENSG00000139618")
  expect_equal(remove_gene_version("ENSG00000139618.10"), "ENSG00000139618")
  expect_equal(remove_gene_version("ENSMUSG00000025902.5"), "ENSMUSG00000025902")
})

# ── Test 2 ── remove_gene_version: RefSeq IDs ───────────────────────────────
test_that("remove_gene_version strips version from RefSeq NM_ and XR_ IDs", {
  skip_if(!global_loaded, "global.R not loaded")
  skip_if_not(exists("remove_gene_version"), "remove_gene_version not found")
  result_nm <- remove_gene_version("NM_000546.5")
  result_xr <- remove_gene_version("XR_007058843.1")
  expect_false(grepl("\\.\\d+$", result_nm))
  expect_false(grepl("\\.\\d+$", result_xr))
  expect_equal(result_nm, "NM_000546")
  expect_equal(result_xr, "XR_007058843")
})

# ── Test 3 ── remove_gene_version: non-matching strings pass through ─────────
test_that("remove_gene_version passes non-matching strings through unchanged", {
  skip_if(!global_loaded, "global.R not loaded")
  skip_if_not(exists("remove_gene_version"), "remove_gene_version not found")
  expect_identical(remove_gene_version("TP53"), "TP53")
  expect_identical(remove_gene_version("7157"), "7157")
})

# ── Test 4 ── cleanGeneSet: upper-case, deduplicate, remove short IDs ────────
test_that("cleanGeneSet uppercases, deduplicates, and removes <2-char IDs", {
  skip_if(!global_loaded, "global.R not loaded")
  skip_if_not(exists("cleanGeneSet"), "cleanGeneSet not found")
  result <- cleanGeneSet(c("a", "tp53", "tp53", "x"))
  expect_equal(result, "TP53")
})

# ── Test 5 ── cleanGeneSet: strips version suffixes from Ensembl IDs ─────────
# Note: cleanGeneSet() calls unique() BEFORE remove_gene_version(), so two
# differently-versioned IDs for the same gene are NOT collapsed to one entry.
# Each versioned ID is stripped individually.
test_that("cleanGeneSet strips version suffix from each Ensembl ID", {
  skip_if(!global_loaded, "global.R not loaded")
  skip_if_not(exists("cleanGeneSet"), "cleanGeneSet not found")
  result <- cleanGeneSet(c("ENSG00000139618.2", "ENSG00000139618.5"))
  # Both entries should have had their versions stripped
  expect_true(all(result == "ENSG00000139618"))
  # Each input produces one output entry (versions stripped, not deduplicated)
  expect_equal(length(result), 2L)
})

# ── Test 6 ── wrap_strings: inserts \n for long strings ─────────────────────
test_that("wrap_strings inserts \\n for long strings; short strings unchanged", {
  skip_if(!global_loaded, "global.R not loaded")
  skip_if_not(exists("wrap_strings"), "wrap_strings not found")
  long_str  <- paste(rep("word", 10), collapse = " ")  # 49 chars
  short_str <- "hello"
  expect_true(grepl("\n", wrap_strings(long_str,  width = 20)))
  expect_false(grepl("\n", wrap_strings(short_str, width = 20)))
})

# ── Test 7 ── mark_duplicates: appends indices to repeated strings ────────────
test_that("mark_duplicates appends ' 1'/' 2'/' 3' to repeated strings; unique unchanged", {
  skip_if(!global_loaded, "global.R not loaded")
  skip_if_not(exists("mark_duplicates"), "mark_duplicates not found")
  input  <- c("aa", "bb", "aa", "cc", "aa")
  result <- mark_duplicates(input)
  expect_true("aa 1" %in% result)
  expect_true("aa 2" %in% result)
  expect_true("aa 3" %in% result)
  expect_true("bb"   %in% result)
  expect_true("cc"   %in% result)
})

# ── Test 8 ── remove_pathway_id: strips prefix for GO categories ─────────────
test_that("remove_pathway_id strips first-word prefix for GOBP; passes through for other categories", {
  skip_if(!global_loaded, "global.R not loaded")
  skip_if_not(exists("remove_pathway_id"), "remove_pathway_id not found")
  pathway_name <- "GO:0006915 apoptosis"
  result_gobp  <- remove_pathway_id(pathway_name, "GOBP")
  result_react  <- remove_pathway_id(pathway_name, "Reactome")
  expect_equal(result_gobp,  "Apoptosis")   # proper() capitalises
  expect_equal(result_react, pathway_name)  # unchanged for non-GO/KEGG
})

# ── Test 9 ── proper: capitalizes first letter only ──────────────────────────
test_that("proper capitalizes first letter and leaves rest unchanged", {
  skip_if(!global_loaded, "global.R not loaded")
  skip_if_not(exists("proper"), "proper not found")
  expect_equal(proper("hello world"), "Hello world")
  expect_equal(proper("hello World"), "Hello World")
})

# ── Test 10 ── extract1: converts GMT-format names to human-readable ─────────
test_that("extract1 converts GMT-format names with >4 words to human-readable form", {
  skip_if(!global_loaded, "global.R not loaded")
  skip_if_not(exists("extract1"), "extract1 not found")
  gmt_name <- "GOBP_mmu_mgi_GO:0000183_chromatin_silencing_at_rDNA"
  result   <- extract1(gmt_name)
  expect_equal(result, "Chromatin silencing at rDNA")
})

# ── Test 11 ── hyperText: generates <a href> tags ─────────────────────────────
test_that("hyperText generates href tags for valid URLs and leaves text unchanged for mismatched lengths", {
  skip_if(!global_loaded, "global.R not loaded")
  skip_if_not(exists("hyperText"), "hyperText not found")
  text   <- c("Apoptosis", "Cell cycle")
  urls   <- c("https://example.com/1", "https://example.com/2")
  result <- hyperText(text, urls)
  expect_true(all(grepl("href=", result)))

  # Mismatched lengths → return unchanged
  result_mismatch <- hyperText(text, urls[1])
  expect_equal(result_mismatch, text)
})

# ── Test 12 ── findSpeciesByIdName / find_taxon_by_id ────────────────────────
test_that("findSpeciesByIdName returns the correct common name from mock orgInfo", {
  skip_if(!global_loaded, "global.R not loaded")
  mock_org <- data.frame(
    id    = c(-9606L, -10090L),
    name  = c("Homo_sapiens", "Mus_musculus"),
    name2 = c("Human", "Mouse"),
    taxon = c(9606L, 10090L),
    stringsAsFactors = FALSE
  )
  # Temporarily override orgInfo
  old_orgInfo <- if (exists("orgInfo")) orgInfo else NULL
  assign("orgInfo", mock_org, envir = .GlobalEnv)
  on.exit({
    if (!is.null(old_orgInfo)) assign("orgInfo", old_orgInfo, envir = .GlobalEnv)
    else if (exists("orgInfo")) rm("orgInfo", envir = .GlobalEnv)
  })

  skip_if_not(exists("findSpeciesByIdName"), "findSpeciesByIdName not found")
  result <- findSpeciesByIdName(-9606L)
  expect_equal(as.character(result), "Human")
})

test_that("find_taxon_by_id returns 9606 for human; handles both taxon and taxon_id columns", {
  skip_if_not(exists("find_taxon_by_id"), "find_taxon_by_id not found")
  org_taxon    <- data.frame(id = c(-9606L, -10090L), taxon    = c(9606L, 10090L))
  org_taxon_id <- data.frame(id = c(-9606L, -10090L), taxon_id = c(9606L, 10090L))
  expect_equal(find_taxon_by_id(-9606L, org_taxon),    9606L)
  expect_equal(find_taxon_by_id(-9606L, org_taxon_id), 9606L)
})

# ── Test 13 ── readGMT: parses a 3-pathway GMT file ──────────────────────────
test_that("readGMT parses a GMT file with 3 pathways correctly", {
  skip_if(!global_loaded, "global.R not loaded")
  skip_if_not(exists("readGMT"), "readGMT not found")
  gmt_content <- paste(
    "PATHWAY_A\tNA\tGENE1\tGENE2\tGENE3",
    "PATHWAY_B\tNA\tGENE4\tGENE5\tGENE6\tGENE7",
    "PATHWAY_C\tNA\tGENE8\tGENE9",
    sep = "\n"
  )
  tmp <- tempfile(fileext = ".gmt")
  writeLines(gmt_content, tmp)
  on.exit(unlink(tmp))
  result <- readGMT(tmp)
  expect_equal(length(result), 3L)
  expect_equal(names(result), c("PATHWAY_A", "PATHWAY_B", "PATHWAY_C"))
  expect_equal(result$PATHWAY_A, c("GENE1", "GENE2", "GENE3"))
  expect_equal(length(result$PATHWAY_B), 4L)
})

# ── Test 14 ── readGMT: handles pathway with 0 genes without error ────────────
test_that("readGMT handles an empty pathway line without error; entry is excluded", {
  skip_if(!global_loaded, "global.R not loaded")
  skip_if_not(exists("readGMT"), "readGMT not found")
  gmt_content <- paste(
    "PATHWAY_FULL\tNA\tGENE1\tGENE2\tGENE3",
    "PATHWAY_EMPTY\tNA",                          # 0 genes
    "PATHWAY_ONE\tNA\tGENE4",                     # 1 gene (also excluded: length ≤ 1)
    sep = "\n"
  )
  tmp <- tempfile(fileext = ".gmt")
  writeLines(gmt_content, tmp)
  on.exit(unlink(tmp))
  expect_no_error({
    result <- readGMT(tmp)
  })
  # Pathways with < 2 genes are excluded
  expect_equal(length(result), 1L)
  expect_equal(names(result), "PATHWAY_FULL")
})

# ── Test 15 ── showGeneIDs: returns data frame with example IDs ───────────────
test_that("showGeneIDs returns data frame with ID Type and Examples columns", {
  skip_if(!global_loaded, "global.R not loaded")
  skip_if_not(exists("showGeneIDs"), "showGeneIDs not found")
  skip_if(!file.exists("./data/data104b/db/Homo_sapiens.db"),
          "Species database not available")
  # The species DB must also have the idIndex table (present in production DBs
  # but not all local fixtures); skip gracefully if it doesn't.
  result <- tryCatch(showGeneIDs(-9606L, nGenes = 5), error = function(e) NULL)
  skip_if(is.null(result), "showGeneIDs DB query failed (idIndex missing)")
  expect_true(is.data.frame(result))
  expect_true("ID Type" %in% colnames(result))
  expect_true("Examples" %in% colnames(result))
  expect_true(nrow(result) > 0)
})

# ── Test 16 ── gmtCategory: returns standard GO categories ───────────────────
test_that("gmtCategory returns at least GOBP and KEGG for human", {
  skip_if(!global_loaded, "global.R not loaded")
  skip_if_not(exists("gmtCategory"), "gmtCategory not found")
  skip_if(!file.exists("./data/data104b/db/Homo_sapiens.db"),
          "Species database not available")
  # Need a valid converted object
  mock_conv <- create_mock_conversion_result("TP53\nBRCA1", species_id = -9606L)
  # Fake species name so grep(gmtFiles) can find it
  mock_conv$species$name <- "Homo_sapiens"
  result <- gmtCategory(mock_conv, -9606L)
  category_values <- unlist(result)
  expect_true("GOBP" %in% category_values || any(grepl("GOBP", names(result))))
})

# ── Test 17 ── gmtCategory: returns error data.frame when converted is NULL ──
test_that("gmtCategory returns error data.frame when converted is NULL", {
  skip_if(!global_loaded, "global.R not loaded")
  skip_if_not(exists("gmtCategory"), "gmtCategory not found")
  result <- gmtCategory(NULL, -9606L)
  expect_true(is.data.frame(result))
  # Should contain the "ID not recognized!" message
  expect_true(any(grepl("ID not recognized", as.character(result))))
})
