library(RSQLite)

datapath <- "./data/data104b/"
db_file <- "Homo_sapiens.db"
full_path <- paste0(datapath, "db/", db_file)

cat("Connecting to:", full_path, "\n")

con <- try(
  DBI::dbConnect(
    drv = RSQLite::dbDriver("SQLite"),
    dbname = full_path,
    flags = RSQLite::SQLITE_RO
  )
)

if (inherits(con, "try-error")) {
  cat("ERROR connecting:\n")
  print(con)
} else {
  cat("SUCCESS!\n")
  cat("Tables:", dbListTables(con), "\n")

  # Try a query
  result <- try(dbGetQuery(con, "SELECT * FROM categories"))
  if (inherits(result, "try-error")) {
    cat("ERROR querying:\n")
    print(result)
  } else {
    cat("Query successful! Rows:", nrow(result), "\n")
    print(result)
  }

  dbDisconnect(con)
}
