# Test setup and helper functions
# This file is automatically loaded before tests run

# Suppress unnecessary warnings during testing
options(warn = -1)

# Helper function to check if databases are available
check_databases <- function() {
  datapath <- "./data/data104b/"
  required_files <- c(
    paste0(datapath, "convertIDs.db"),
    paste0(datapath, "db/Homo_sapiens.db"),
    paste0(datapath, "db/Mus_musculus.db")
  )

  missing <- required_files[!file.exists(required_files)]

  if (length(missing) > 0) {
    warning("Missing database files:\n", paste(missing, collapse = "\n"))
    return(FALSE)
  }

  return(TRUE)
}

# Helper to create mock orgInfo data
create_mock_orginfo <- function() {
  data.frame(
    id = c(-9606, -10090, -3702),
    name = c("Homo_sapiens", "Mus_musculus", "Arabidopsis_thaliana"),
    name2 = c("Human", "Mouse", "Arabidopsis thaliana"),
    ensembl_dataset = c("hsapiens_gene_ensembl", "mmusculus_gene_ensembl", "athaliana_gene_ensembl"),
    group = c("Ensembl", "Ensembl", "Plants"),
    taxon = c(9606, 10090, 3702),
    genes = c(20000, 22000, 27000),
    stringsAsFactors = FALSE
  )
}

# Helper to safely disconnect database
safe_disconnect <- function(con) {
  if (!is.null(con) && inherits(con, "SQLiteConnection")) {
    try(DBI::dbDisconnect(con), silent = TRUE)
  }
}

# Helper to clean up global environment after tests
cleanup_globals <- function() {
  vars_to_remove <- c("convert", "datapath", "orgInfo", "speciesChoice",
                      "gmtFiles", "keggSpeciesID", "convert_species")
  for (var in vars_to_remove) {
    if (exists(var, envir = .GlobalEnv)) {
      rm(list = var, envir = .GlobalEnv)
    }
  }
}

cat("Test helpers loaded successfully\n")
