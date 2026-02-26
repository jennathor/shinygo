# Mock helpers for testing without full databases
# These simulate database responses for testing code structure and logic

# ── Original helpers (preserved) ─────────────────────────────────────────────

#' Create mock gene conversion result
#' Simulates what convertID() would return
create_mock_conversion_result <- function(gene_ids, species_id = -9606) {
  original_ids <- unlist(strsplit(toupper(gene_ids), "\t| |\n|\\,|;"))
  original_ids <- unique(gsub(" ", "", original_ids))
  original_ids <- original_ids[nchar(original_ids) > 1]

  if (length(original_ids) == 0) return(NULL)

  converted_ids <- if (any(grepl("^ENSG", original_ids))) {
    original_ids
  } else {
    paste0("ENSG", sprintf("%011d", seq_along(original_ids)))
  }

  conversion_table <- data.frame(
    User_input       = original_ids,
    ensembl_gene_id  = converted_ids,
    Species          = "Human",
    stringsAsFactors = FALSE
  )
  species_info <- data.frame(
    id               = species_id,
    name             = "Homo_sapiens",
    name2            = "Human",
    ensembl_dataset  = "hsapiens_gene_ensembl",
    group            = "Ensembl",
    taxon            = abs(species_id),
    genes            = 20000,
    stringsAsFactors = FALSE
  )
  species_matched <- data.frame(
    "Matched Species (%genes)" = paste0("Human (", length(original_ids), ")"),
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )
  list(
    originalIDs    = original_ids,
    IDs            = converted_ids,
    species        = species_info,
    speciesMatched = species_matched,
    conversionTable = conversion_table
  )
}

#' Mock convertID function
mock_convertID <- function(query, selectOrg) {
  if (is.null(query) || nchar(trimws(query)) == 0) return(NULL)
  create_mock_conversion_result(query, selectOrg)
}

#' Mock database query (handles missing tables gracefully)
mock_db_query <- function(conn, query) {
  query_lower <- tolower(query)
  if (grepl("from mapping", query_lower) && grepl("where id in", query_lower)) {
    return(data.frame(id = "TEST_GENE_1", ens = "ENSG00000000001",
                      idType = "1", stringsAsFactors = FALSE))
  }
  if (grepl("from idindex", query_lower)) {
    return(data.frame(id = c("1","2","3"),
                      idType = c("ensembl_gene_id","hgnc_symbol","entrezgene_id"),
                      stringsAsFactors = FALSE))
  }
  data.frame()
}

#' Check if database has a specific table
has_table <- function(conn, table_name) {
  if (is.null(conn) || inherits(conn, "try-error")) return(FALSE)
  tables <- try(DBI::dbListTables(conn), silent = TRUE)
  if (inherits(tables, "try-error")) return(FALSE)
  tolower(table_name) %in% tolower(tables)
}

#' Check if we have full database structure
has_full_database <- function() {
  test_db <- "./data/data104b/db/Homo_sapiens.db"
  if (!file.exists(test_db)) return(FALSE)
  con <- try(DBI::dbConnect(drv = RSQLite::dbDriver("SQLite"), dbname = test_db,
                            flags = RSQLite::SQLITE_RO), silent = TRUE)
  if (inherits(con, "try-error")) return(FALSE)
  has_mapping <- has_table(con, "mapping")
  DBI::dbDisconnect(con)
  has_mapping
}

#' Run test with real or mock database
with_mock_or_real <- function(test_code, mock_code) {
  if (has_full_database()) eval(test_code) else eval(mock_code)
}

#' Create mock geneInfo dataframe (minimal version)
create_mock_gene_info <- function(gene_ids) {
  n_genes <- length(gene_ids)
  data.frame(
    ensembl_gene_id  = gene_ids,
    entrezgene_id    = sample(1000:9999, n_genes),
    symbol           = paste0("GENE", seq_len(n_genes)),
    chromosome_name  = sample(c(as.character(1:22), "X", "Y"), n_genes, replace = TRUE),
    gene_biotype     = "protein_coding",
    description      = paste("Mock gene", seq_len(n_genes), "description"),
    Set              = rep("List", n_genes),
    duplicated       = rep(FALSE, n_genes),
    stringsAsFactors = FALSE
  )
}

#' Setup test environment with mocked global variables
setup_mock_globals <- function() {
  mock_orgInfo <- data.frame(
    id               = c(-9606L, -10090L, -3702L),
    name             = c("Homo_sapiens","Mus_musculus","Arabidopsis_thaliana"),
    name2            = c("Human","Mouse","Arabidopsis thaliana"),
    ensembl_dataset  = c("hsapiens_gene_ensembl","mmusculus_gene_ensembl","athaliana_gene_ensembl"),
    group            = c("Ensembl","Ensembl","Plants"),
    taxon            = c(9606L, 10090L, 3702L),
    genes            = c(20000L, 22000L, 27000L),
    stringsAsFactors = FALSE
  )
  assign("orgInfo", mock_orgInfo, envir = .GlobalEnv)
  assign("speciesChoice", setNames(as.list(mock_orgInfo$id), mock_orgInfo$name2), envir = .GlobalEnv)
  assign("datapath", "./data/data104b/", envir = .GlobalEnv)
  invisible(TRUE)
}

#' Cleanup mocked global variables
cleanup_mock_globals <- function() {
  vars_to_remove <- c("orgInfo", "speciesChoice", "datapath", "convert")
  for (var in vars_to_remove) {
    if (exists(var, envir = .GlobalEnv)) suppressWarnings(rm(list = var, envir = .GlobalEnv))
  }
  invisible(TRUE)
}

# ── New v2 helpers ────────────────────────────────────────────────────────────

#' Create a mock enrichment result list (mirrors FindOverlap() output)
#'
#' The `.__real_names` attribute records the exact column names that a real
#' FindOverlap() call produces; test-mock_contracts.R uses it to detect drift.
create_mock_enrichment_result <- function(n = 5) {
  set.seed(42)
  df <- data.frame(
    `Enrichment FDR`  = 10^(-seq(3, 3 + n - 1)),
    nGenes            = as.integer(seq(5, 5 * n, by = 5)),
    `Pathway Genes`   = as.integer(seq(20, 20 * n, by = 20)),
    `Fold Enrichment` = seq(2, 2 * n, by = 2),
    Pathway           = paste0("Pathway_", seq_len(n)),
    URL               = paste0("http://example.com/", seq_len(n)),
    Genes             = replicate(n, paste(sample(LETTERS, 3), collapse = " ")),
    check.names       = FALSE,
    stringsAsFactors  = FALSE
  )
  result <- list(
    x               = df,
    groupings       = data.frame(
      Group   = paste0("G", seq_len(n)),
      Pathway = df$Pathway,
      stringsAsFactors = FALSE
    ),
    categoryChoices = c("GOBP", "GOCC", "GOMF", "KEGG")
  )
  attr(result, ".__real_names") <- c(
    "Enrichment FDR", "nGenes", "Pathway Genes",
    "Fold Enrichment", "Pathway", "URL", "Genes"
  )
  result
}

#' Create a mock sidebar_values reactive list
#'
#' Returns plain functions (not shiny::reactive()) so the object is usable both
#' inside testServer and in ordinary unit tests.
create_mock_sidebar_values <- function(goButton = 1L) {
  list(
    goButton             = function() goButton,
    selectOrg            = function() -9606L,
    selectGO             = function() "GOBP",
    input_text           = function() "BRCA1\nTP53\nEGFR",
    input_text_b         = function() "",
    maxTerms             = function() 20L,
    minFDR               = function() 0.05,
    abbreviatePathway    = function() TRUE,
    removeRedundantSets  = function() TRUE,
    show_pathway_id      = function() FALSE,
    minSetSize           = function() 2L,
    maxSetSize           = function() 2000L,
    gene_count_pathwaydb = function() FALSE
  )
}

#' Create a mock enrichment_values reactive list
create_mock_enrichment_values <- function(result = NULL) {
  if (is.null(result)) result <- create_mock_enrichment_result()
  list(
    significantOverlaps    = function() result,
    significantOverlapsAll = function() result,
    SortPathways           = function() "Select by FDR, sort by Fold Enrichment"
  )
}

#' Create a full mock gene-info dataframe (mirrors geneInfo() output)
#'
#' The `.__real_names` attribute records the columns geneInfo() actually returns.
create_mock_gene_info_full <- function(n = 10) {
  set.seed(42)
  ensembl_ids <- paste0("ENSG", formatC(seq_len(n), width = 11, flag = "0"))
  df <- data.frame(
    ensembl_gene_id  = ensembl_ids,
    entrezgene_id    = 7000L + seq_len(n),
    symbol           = paste0("GENE", seq_len(n)),
    chromosome_name  = sample(c(as.character(1:22), "X", "Y"), n, replace = TRUE),
    gene_biotype     = "protein_coding",
    description      = paste("Mock gene", seq_len(n), "function"),
    start_position   = as.integer(seq(1e6, 10e6, length.out = n)),
    nExons           = sample(5L:20L, n, replace = TRUE),
    transcript_count = sample(2L:8L,  n, replace = TRUE),
    Set              = c(rep("List", ceiling(n / 2)), rep("Genome", floor(n / 2))),
    duplicated       = FALSE,
    stringsAsFactors = FALSE
  )
  attr(df, ".__real_names") <- c(
    "ensembl_gene_id", "entrezgene_id", "symbol", "chromosome_name",
    "gene_biotype", "description", "start_position", "nExons",
    "transcript_count", "Set", "duplicated"
  )
  df
}

#' Create a mock conversionTableData dataframe (mirrors what server.R constructs)
#'
#' This is the merged gene-info + conversion table passed to mod_07_genes and
#' mod_10_string as `conversionTableData`.
create_mock_conversion_table_data <- function(n = 5) {
  set.seed(42)
  ensembl_ids <- paste0("ENSG", formatC(seq_len(n), width = 11, flag = "0"))
  df <- data.frame(
    User_input          = paste0("GENE", seq_len(n)),
    `Ensembl Gene ID`   = ensembl_ids,
    Entrez              = as.character(7000L + seq_len(n)),
    Symbol              = paste0("GENE", seq_len(n)),
    Description         = paste0("Mock gene function; [Source:HGNC]"),
    Chr                 = sample(c(as.character(1:22), "X"), n, replace = TRUE),
    Type                = "protein_coding",
    Species             = "Human",
    Set                 = "List",
    stringsAsFactors    = FALSE,
    check.names         = FALSE
  )
  attr(df, ".__real_names") <- c(
    "User_input", "Ensembl Gene ID", "Entrez", "Symbol",
    "Description", "Chr", "Type", "Species", "Set"
  )
  df
}

#' Create a mock promoter database data frame
#'
#' Mimics the TF × gene score matrix returned by an in-memory promoter DB.
create_mock_promoter_db <- function(n_tf = 5, n_genes = 20) {
  set.seed(42)
  gene_ids <- paste0("ENSG", formatC(seq_len(n_genes), width = 11, flag = "0"))
  tf_ids   <- paste0("TF", seq_len(n_tf))
  # Score matrix: n_genes × n_tf
  scores <- matrix(rnorm(n_genes * n_tf, mean = 5, sd = 2),
                   nrow = n_genes, ncol = n_tf,
                   dimnames = list(gene_ids, tf_ids))
  as.data.frame(scores)
}

#' Create a mock GMT-like named list (mirrors readGMT() output)
create_mock_gmt <- function(n_pathways = 3, genes_per = 5) {
  set.seed(42)
  pathways <- setNames(
    lapply(seq_len(n_pathways), function(i) {
      paste0("GENE", sample(100, genes_per))
    }),
    paste0("PATHWAY_", seq_len(n_pathways))
  )
  pathways
}

cat("✓ Mock helpers loaded\n")
