# ShinyGO Test Suite

Comprehensive tests for all ShinyGO modules and core functions.
See [`TESTING.md`](../TESTING.md) in the project root for full documentation.

## Quick Start

```bash
# From project root
Rscript run_tests.R                          # all tests
Rscript run_tests.R test-mod_02_enrichment.R # single file
```

## File Structure

```
tests/
├── README.md                     # This file
├── testthat.R                    # testthat entry point
└── testthat/
    ├── helper-setup.R            # check_databases(), safe_disconnect(), cleanup_globals()
    ├── helper-mocks.R            # mock constructors for all shared data objects
    ├── helper-db.R               # creates SQLite fixture DBs on first run
    ├── fixtures/
    │   ├── mock_db.sqlite        # mirrors convertIDs.db schema (auto-created)
    │   └── db/
    │       └── Homo_sapiens.db   # mirrors species DB schema (auto-created)
    │
    ├── test-mock_contracts.R         # mock drift prevention (.__real_names checks)
    ├── test-fct_statistics.R         # phyper, p.adjust, fold enrichment formula
    ├── test-global_utils.R           # pure utility functions from global.R
    ├── test-fct_helper_functions.R   # fct_08_plots.R and fct_09_genome.R functions
    ├── test-mod_01_sidebar.R         # sidebar module
    ├── test-mod_02_enrichment.R      # enrichment table module
    ├── test-mod_03_chart.R           # enrichment chart module
    ├── test-mod_04_tree.R            # hierarchical tree module
    ├── test-mod_05_network.R         # interactive network module
    ├── test-mod_06_kegg.R            # KEGG pathway module
    ├── test-mod_07_genes.R           # gene table module
    ├── test-mod_08_plots.R           # gene characteristics plots module
    ├── test-mod_09_genome.R          # genome browser module
    ├── test-mod_10_string.R          # STRING-DB module
    ├── test-mod_11_about.R           # About tab
    ├── test-mod_13_groups.R          # gene grouping module
    └── test-requires_external.R      # fixture DB integration + stats pipeline tests
```

## Mock Constructors

All constructors live in `helper-mocks.R` and return objects that mimic real app data:

| Constructor | Mimics | Key attribute |
|-------------|--------|---------------|
| `create_mock_conversion_result(ids)` | `convertID()` output | — |
| `create_mock_enrichment_result(n)` | `FindOverlap()` output (`list(x, groupings, categoryChoices)`) | `.__real_names` on `x` |
| `create_mock_sidebar_values(goButton)` | `mod_01_sidebar_server()` return value | plain functions (not `reactive()`) |
| `create_mock_enrichment_values(result)` | `mod_02_enrichment_server()` return value | plain functions |
| `create_mock_gene_info_full(n)` | `geneInfo()` output | `.__real_names` |
| `create_mock_conversion_table_data(n)` | `conversionTableData` reactive | `.__real_names` |
| `create_mock_gmt(n_pathways, genes_per)` | `readGMT()` output | named list |

The `.__real_names` attribute is checked in `test-mock_contracts.R` to catch drift
between mock data and the column names the modules actually expect.

## Fixture Databases

`helper-db.R` auto-creates minimal SQLite databases on first run:

**`fixtures/mock_db.sqlite`** (mirrors `convertIDs.db`):
- `orgInfo`: human, mouse, Arabidopsis
- `mapping`: 10 human genes (TP53, BRCA1, EGFR, …), 10 mouse genes
- `idIndex`: hgnc_symbol, ensembl_gene_id, entrezgene_id
- `quotes`, `GO` level table

**`fixtures/db/Homo_sapiens.db`** (mirrors species DB):
- `categories`: GOBP, GOCC, GOMF, KEGG
- `pathwayInfo`: 5 GOBP pathways (DNA repair, Apoptosis, Cell cycle, …)
- `pathway`: gene↔pathwayID mappings for 10 human Ensembl IDs
- `geneInfo`: 10 protein-coding genes with chr, position, exon counts

Delete the fixture files and re-run tests to regenerate them.

## Adding Tests for a New Module

1. Create `test-mod_XX_name.R` following the existing pattern:
   ```r
   global_loaded <- tryCatch({
     suppressMessages(suppressWarnings(source("global.R", local = FALSE)))
     source("R/mod_XX_name.R", local = FALSE)
     TRUE
   }, error = function(e) FALSE)

   test_that("mod_XX_name_ui contains expected outputs", {
     skip_if(!global_loaded, "modules not loaded")
     ui   <- mod_XX_name_ui("testXX")
     html <- as.character(ui)
     expect_true(grepl("expectedOutput", html))
   })
   ```

2. For functions that call `bsModal()` (requires shinyBS), guard the UI call:
   ```r
   ui <- tryCatch({
     suppressPackageStartupMessages(library(shinyBS))
     mod_XX_name_modal_ui("id")
   }, error = function(e) NULL)
   skip_if(is.null(ui), "shinyBS::bsModal not available")
   ```

3. For `testServer()` with a module that takes extra args, pass mock reactive lists:
   ```r
   sidebar <- create_mock_sidebar_values(goButton = 1L)
   testServer(mod_XX_name_server,
     args = list(sidebar_values = sidebar),
     { session$setInputs(someInput = "value")
       expect_true(...)
     }
   )
   ```
