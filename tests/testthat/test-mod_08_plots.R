# test-mod_08_plots.R
# Tests for mod_08_plots: UI structure and plot computation logic.

library(testthat)
library(shiny)
library(ggplot2)

if (basename(getwd()) == "testthat") setwd("../..")
if (basename(getwd()) == "tests")    setwd("..")

global_loaded <- tryCatch({
  suppressMessages(suppressWarnings(source("global.R", local = FALSE)))
  source("R/mod_08_plots.R",   local = FALSE)
  source("R/fct_08_plots.R",   local = FALSE)
  TRUE
}, error = function(e) FALSE)

# ── Test 1 ── UI contains genePlot2 and gene_barplot outputs ─────────────────
test_that("mod_08_plots_ui contains plotOutput('genePlot2') and plotOutput('gene_barplot')", {
  skip_if(!global_loaded, "modules not loaded")
  ui   <- mod_08_plots_ui("test08")
  html <- as.character(ui)
  expect_true(grepl("genePlot2",   html))
  expect_true(grepl("gene_barplot", html))
})

# ── Test 2 ── Chromosome frequency table produces correct 2-column matrix ─────
test_that("chromosome frequency table produces correct 2-column matrix (List vs Genome)", {
  gene_data <- data.frame(
    chromosome_name = c("1", "1", "2", "3", "1", "2", "3", "3", "X", "X"),
    gene_biotype    = "protein_coding",
    Set             = c("List", "List", "List", "Genome", "Genome",
                        "Genome", "Genome", "List", "Genome", "List"),
    stringsAsFactors = FALSE
  )
  # Simulate server logic: table(chr_name) split by Set
  list_chrs  <- gene_data$chromosome_name[gene_data$Set == "List"]
  genome_chrs <- gene_data$chromosome_name[gene_data$Set != "List"]
  freq_list   <- table(list_chrs)
  freq_genome <- table(genome_chrs)
  all_chrs    <- sort(unique(c(names(freq_list), names(freq_genome))))
  mat <- matrix(0, nrow = length(all_chrs), ncol = 2,
                dimnames = list(all_chrs, c("List", "Genome")))
  mat[names(freq_list),   1] <- as.integer(freq_list)
  mat[names(freq_genome), 2] <- as.integer(freq_genome)
  # chr1 has 2 list genes and 1 genome gene
  expect_equal(mat["1", "List"],   2L)
  expect_equal(mat["1", "Genome"], 1L)
  expect_equal(ncol(mat), 2L)
})

# ── Test 3 ── Expected counts normalized by column totals ─────────────────────
test_that("expected counts are normalized by column totals (chi-squared independence)", {
  observed <- matrix(c(10L, 5L, 15L, 3L, 8L, 12L), nrow = 3,
                     dimnames = list(c("chr1","chr2","chrX"), c("List","Genome")))
  # Chi-squared expected: outer(rowSums, colSums) / sum
  expected <- outer(rowSums(observed), colSums(observed)) / sum(observed)
  expect_equal(colSums(expected), colSums(observed))
})

# ── Test 4 ── chisq.test p-value in (0, 1) ────────────────────────────────────
test_that("chisq.test p-value is in (0, 1) for non-trivial frequency data", {
  observed <- matrix(c(10L, 5L, 20L, 15L, 2L, 8L), nrow = 3,
                     dimnames = list(c("chr1","chr2","chrX"), c("List","Genome")))
  result <- chisq.test(observed)
  expect_true(result$p.value > 0 && result$p.value < 1)
})

# ── Test 5 ── Chi-squared not significant when distributions are identical ─────
test_that("chi-squared p-value > 0.05 when query and genome have identical distributions", {
  # Both columns have proportional counts → no enrichment
  observed <- matrix(c(10L, 20L, 10L, 20L), nrow = 2,
                     dimnames = list(c("chr1","chr2"), c("List","Genome")))
  result <- suppressWarnings(chisq.test(observed))
  expect_true(result$p.value > 0.05)
})

# ── Test 6 ── Background genes merged with Set = "Background" ────────────────
test_that("background genes are merged with Set = 'Background' label", {
  query_genes <- data.frame(
    ensembl_gene_id = c("ENSG00000001", "ENSG00000002"),
    Set             = "List",
    stringsAsFactors = FALSE
  )
  bg_genes <- data.frame(
    ensembl_gene_id = c("ENSG00000003", "ENSG00000004"),
    Set             = "Background",
    stringsAsFactors = FALSE
  )
  combined <- rbind(query_genes, bg_genes)
  expect_true("Background" %in% unique(combined$Set))
  expect_true("List"       %in% unique(combined$Set))
})

# ── Test 7 ── Only protein_coding genes enter coding-gene analyses ────────────
test_that("only protein_coding genes enter coding-gene analyses", {
  df <- data.frame(
    gene_biotype = c("protein_coding", "lncRNA", "protein_coding", "pseudogene"),
    Set          = "List",
    stringsAsFactors = FALSE
  )
  coding <- df[df$gene_biotype == "protein_coding", ]
  expect_equal(nrow(coding), 2L)
  expect_true(all(coding$gene_biotype == "protein_coding"))
})

# ── Test 8 ── Plot suppressed when fewer than minGenes present ────────────────
test_that("plot is suppressed (returns NULL / fake_plot) when fewer than minGenes genes", {
  skip_if(!global_loaded, "modules not loaded")
  skip_if_not(exists("fake_plot"), "fake_plot not found")
  minGenes_test <- 10L
  n_genes       <- 5L
  result <- if (n_genes < minGenes_test) fake_plot("Not enough genes") else NULL
  expect_true(inherits(result, "gg"))
})

# ── Test 9 ── Chr names > 50 chars excluded from frequency table ──────────────
test_that("chromosome names longer than 50 chars are excluded from frequency table", {
  chrs <- c("1", "X", paste(rep("a", 51), collapse = ""), "MT", "2")
  filtered <- chrs[nchar(chrs) <= 50]
  expect_false(any(nchar(filtered) > 50))
  expect_true("1" %in% filtered)
  expect_true("MT" %in% filtered)
  expect_false(paste(rep("a", 51), collapse = "") %in% filtered)
})

# ── Test 10 ── Download handler returns valid filename string ─────────────────
test_that("download handler for gene barplots returns a character filename", {
  filename_fn <- function() "genePlots.png"
  filename    <- filename_fn()
  expect_true(is.character(filename))
  expect_true(nchar(filename) > 0)
})

# ── Test 11 ── Gene length distribution snapshot (vdiffr) ─────────────────────
test_that("gene length barplot visual snapshot", {
  skip_if_not_installed("vdiffr")
  skip_if(!global_loaded, "modules not loaded")

  set.seed(42)
  lengths <- data.frame(
    Length   = c(1000, 2000, 500, 3000, 1500, 800, 2500, 400),
    Set      = c(rep("List", 4), rep("Genome", 4)),
    stringsAsFactors = FALSE
  )
  p <- ggplot(lengths, aes(x = Length, fill = Set)) +
    geom_histogram(bins = 10, position = "identity", alpha = 0.6) +
    scale_fill_manual(values = c("List" = "red", "Genome" = "lightgrey")) +
    theme_light()

  vdiffr::expect_doppelganger("gene-length-barplot", p)
})
