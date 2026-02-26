# helper-db.R
# Creates a minimal mock SQLite database that mirrors ShinyGO's schema.
# This file is automatically loaded by testthat before any tests run.
# Provides tests/testthat/fixtures/mock_db.sqlite without needing
# the full production database.

local({
  db_path <- "tests/testthat/fixtures/mock_db.sqlite"

  if (!file.exists(db_path)) {
    dir.create(dirname(db_path), showWarnings = FALSE, recursive = TRUE)

    if (!requireNamespace("DBI", quietly = TRUE) ||
        !requireNamespace("RSQLite", quietly = TRUE)) {
      message("DBI/RSQLite not available; mock_db.sqlite not created.")
      return(invisible(NULL))
    }

    con <- DBI::dbConnect(RSQLite::SQLite(), db_path)

    # orgInfo table -------------------------------------------------------
    DBI::dbExecute(con, "
      CREATE TABLE IF NOT EXISTS orgInfo (
        id INTEGER PRIMARY KEY,
        name TEXT,
        name2 TEXT,
        ensembl_dataset TEXT,
        grp TEXT,
        taxon INTEGER,
        genes INTEGER
      )
    ")
    DBI::dbExecute(con, "INSERT INTO orgInfo VALUES
      (-9606, 'Homo_sapiens',         'Human',                  'hsapiens_gene_ensembl',  'Ensembl', 9606,  20000),
      (-10090,'Mus_musculus',          'Mouse',                  'mmusculus_gene_ensembl', 'Ensembl', 10090, 22000),
      (-3702, 'Arabidopsis_thaliana',  'Arabidopsis thaliana',   'athaliana_gene_ensembl', 'Plants',  3702,  27000)
    ")

    # mapping table -------------------------------------------------------
    DBI::dbExecute(con, "
      CREATE TABLE IF NOT EXISTS mapping (
        id TEXT,
        ensembl_gene_id TEXT,
        species_id TEXT,
        idtype TEXT
      )
    ")
    # 10 human Ensembl IDs with HGNC symbols
    human_ensembl <- paste0("ENSG", formatC(1:10, width = 11, flag = "0"))
    human_symbols <- c("TP53","BRCA1","EGFR","MYC","KRAS",
                       "PTEN","CDK2","MDM2","RB1","AKT1")
    for (i in seq_along(human_ensembl)) {
      DBI::dbExecute(con,
        sprintf("INSERT INTO mapping VALUES ('%s', '%s', '-9606', '1')",
                human_symbols[i], human_ensembl[i]))
    }

    # mouse (species -10090)
    mouse_ensembl <- paste0("ENSMUSG", formatC(1:10, width = 11, flag = "0"))
    mouse_symbols <- c("Trp53","Brca1","Egfr","Myc","Kras",
                       "Pten","Cdk2","Mdm2","Rb1","Akt1")
    for (i in seq_along(mouse_ensembl)) {
      DBI::dbExecute(con,
        sprintf("INSERT INTO mapping VALUES ('%s', '%s', '-10090', '1')",
                mouse_symbols[i], mouse_ensembl[i]))
    }

    # idIndex table -------------------------------------------------------
    DBI::dbExecute(con, "
      CREATE TABLE IF NOT EXISTS idIndex (
        id TEXT PRIMARY KEY,
        idType TEXT
      )
    ")
    DBI::dbExecute(con, "INSERT INTO idIndex VALUES
      ('1', 'hgnc_symbol'),
      ('2', 'ensembl_gene_id'),
      ('3', 'entrezgene_id')
    ")

    # quotes table --------------------------------------------------------
    DBI::dbExecute(con, "
      CREATE TABLE IF NOT EXISTS quotes (
        quotes TEXT,
        author TEXT
      )
    ")
    DBI::dbExecute(con, "INSERT INTO quotes VALUES
      ('Test quote one.', 'Author A'),
      ('Test quote two.', 'Author B')
    ")

    # GO levels table (for level2Terms in global.R) -----------------------
    DBI::dbExecute(con, "
      CREATE TABLE IF NOT EXISTS GO (
        id TEXT,
        level INTEGER,
        GO TEXT
      )
    ")
    DBI::dbExecute(con, "INSERT INTO GO VALUES
      ('GO:0008150', 2, 'biological_process'),
      ('GO:0009987', 3, 'biological_process'),
      ('GO:0003674', 2, 'molecular_function'),
      ('GO:0005575', 2, 'cellular_component')
    ")

    DBI::dbDisconnect(con)
    message("Created mock_db.sqlite fixture at: ", db_path)
  }

  # -----------------------------------------------------------------------
  # Create a minimal per-species pathway database for human
  # Used by test-requires_external.R blocks C1, C2, etc.
  # -----------------------------------------------------------------------
  species_dir <- "tests/testthat/fixtures/db"
  dir.create(species_dir, showWarnings = FALSE, recursive = TRUE)
  species_db <- file.path(species_dir, "Homo_sapiens.db")

  if (!file.exists(species_db)) {
    if (!requireNamespace("DBI", quietly = TRUE) ||
        !requireNamespace("RSQLite", quietly = TRUE)) {
      return(invisible(NULL))
    }
    con2 <- DBI::dbConnect(RSQLite::SQLite(), species_db)

    # categories table ----------------------------------------------------
    DBI::dbExecute(con2, "
      CREATE TABLE IF NOT EXISTS categories (category TEXT)
    ")
    DBI::dbExecute(con2, "INSERT INTO categories VALUES
      ('GOBP'), ('GOCC'), ('GOMF'), ('KEGG')
    ")

    # pathwayInfo table ---------------------------------------------------
    DBI::dbExecute(con2, "
      CREATE TABLE IF NOT EXISTS pathwayInfo (
        id TEXT, n INTEGER, description TEXT, memo TEXT, golevel INTEGER
      )
    ")
    # 5 GOBP pathways each containing some of the 10 test genes
    DBI::dbExecute(con2, "INSERT INTO pathwayInfo VALUES
      ('GO:0000001', 5,  'DNA repair',          'http://amigo.geneontology.org/amigo/term/GO:0000001', NULL),
      ('GO:0000002', 4,  'Apoptosis',           'http://amigo.geneontology.org/amigo/term/GO:0000002', 2),
      ('GO:0000003', 3,  'Cell cycle',          'http://amigo.geneontology.org/amigo/term/GO:0000003', 3),
      ('GO:0000004', 6,  'Signal transduction', 'http://amigo.geneontology.org/amigo/term/GO:0000004', NULL),
      ('GO:0000005', 2,  'Transcription',       'http://amigo.geneontology.org/amigo/term/GO:0000005', NULL)
    ")

    # pathway table -------------------------------------------------------
    # Map Ensembl IDs to pathways
    human_ensembl <- paste0("ENSG", formatC(1:10, width = 11, flag = "0"))

    DBI::dbExecute(con2, "
      CREATE TABLE IF NOT EXISTS pathway (
        gene TEXT, pathwayID TEXT, category TEXT
      )
    ")
    pathway_data <- rbind(
      data.frame(gene = human_ensembl[1:5],  pathwayID = "GO:0000001", category = "GOBP"),
      data.frame(gene = human_ensembl[1:4],  pathwayID = "GO:0000002", category = "GOBP"),
      data.frame(gene = human_ensembl[3:5],  pathwayID = "GO:0000003", category = "GOBP"),
      data.frame(gene = human_ensembl[4:9],  pathwayID = "GO:0000004", category = "GOBP"),
      data.frame(gene = human_ensembl[7:8],  pathwayID = "GO:0000005", category = "GOBP")
    )
    DBI::dbWriteTable(con2, "pathway", pathway_data, append = TRUE)

    # geneInfo table ------------------------------------------------------
    DBI::dbExecute(con2, "
      CREATE TABLE IF NOT EXISTS geneInfo (
        ensembl_gene_id TEXT,
        entrezgene_id INTEGER,
        symbol TEXT,
        chromosome_name TEXT,
        gene_biotype TEXT,
        description TEXT,
        start_position INTEGER,
        nExons INTEGER,
        transcript_count INTEGER
      )
    ")
    gene_info <- data.frame(
      ensembl_gene_id  = human_ensembl,
      entrezgene_id    = 7001:7010,
      symbol           = c("TP53","BRCA1","EGFR","MYC","KRAS",
                           "PTEN","CDK2","MDM2","RB1","AKT1"),
      chromosome_name  = c("17","17","7","8","12","10","12","12","13","14"),
      gene_biotype     = "protein_coding",
      description      = paste("Mock gene", 1:10),
      start_position   = seq(1e6, 10e6, length.out = 10),
      nExons           = sample(5:20, 10, replace = TRUE),
      transcript_count = sample(2:8, 10, replace = TRUE),
      stringsAsFactors = FALSE
    )
    DBI::dbWriteTable(con2, "geneInfo", gene_info, append = TRUE)

    DBI::dbDisconnect(con2)
    message("Created fixture species DB at: ", species_db)
  }
})

cat("✓ Mock DB fixture helper loaded\n")
