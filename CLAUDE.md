# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ShinyGO is an R Shiny web application for functional gene set enrichment analysis. It accepts gene lists, maps them to pathway databases (GO, KEGG, Reactome, etc.), performs hypergeometric enrichment tests, and generates interactive visualizations. It supports 600+ organisms via SQLite databases.

## Files to know
Read the 'README.md' and 'TESTING.md' files in the root for repo and app context. Read the tests/README.md file for testing context.

## Running Tests

```bash
# Run all tests
Rscript run_tests.R

# Run a specific test file
Rscript run_tests.R test-mod_01_sidebar.R
```

Tests live in `tests/testthat/`. Helpers are in `helper-setup.R`, `helper-mocks.R`, and `helper-db.R` (auto-creates SQLite fixtures). See `TESTING.md` for the full test inventory.

## Architecture

### File Structure

- **`global.R`** — Database paths, global constants, and all core analysis functions (`convertID`, `geneInfo`, `FindOverlap`, visualization helpers). This is the largest and most critical file.
- **`server.R`** — Orchestrates module initialization, defines shared reactives (`converted`, `geneInfoLookup`, `converted_background`, `geneInfoLookup_background`, `conversionTableData`), and wires reactive values between modules.
- **`ui.R`** — Top-level tab layout; delegates all tab UI to module functions.
- **`R/mod_0N_*.R`** — Shiny modules (UI + server in same file). Each tab is one module:
  - `mod_01_sidebar` — gene input, species selection, parameters
  - `mod_02_enrichment` — enrichment table
  - `mod_03_chart` — enrichment chart (lollipop/dotplot/barplot)
  - `mod_04_tree` — hierarchical clustering tree
  - `mod_05_network` — interactive pathway network
  - `mod_06_kegg` — KEGG pathway diagrams
  - `mod_07_genes` — gene conversion table
  - `mod_08_plots` — gene characteristic density/bar plots
  - `mod_09_genome` — genome browser / chromosome plot
  - `mod_10_string` — STRING-DB integration and PPI network
  - `mod_11_about` — About tab
  - `mod_12_promoter` — promoter motif analysis (hidden — uncomment in `ui.R` and `server.R` to enable)
  - `mod_13_groups` — gene grouping by high-level GO terms
- **`R/fct_0N_*.R`** — Supporting functions for specific modules:
  - `fct_01_sidebar.R` — sidebar helpers
  - `fct_06_kegg.R` — KEGG pathway mapping
  - `fct_08_plots.R` — plot helpers: `densMode`, `refine_ggplot2`, `fake_plot`, `mark_significance`
  - `fct_09_genome.R` — `genome_sliding_window()` for chromosome enrichment analysis

### Module Communication Pattern

Modules pass reactive values down the chain:

```
mod_01_sidebar  →  sidebar_values  (gene input, species, parameters)
      ↓
mod_02_enrichment  →  enrichment_values  (filtered enrichment results)
      ↓
mod_03_chart  →  chart_values  (chart type, ggplot2 theme)
      ↓
mod_03_chart, mod_04_tree, mod_05_network, mod_06_kegg, mod_07_genes,
mod_08_plots, mod_09_genome, mod_10_string, mod_13_groups
```

Shared reactives defined in `server.R` and passed explicitly as arguments to modules that need them: `converted`, `geneInfoLookup`, `converted_background`, `geneInfoLookup_background`, `conversionTableData`.

Each module follows the pattern:
```r
mod_XX_feature_ui("id")       # Returns UI elements
mod_XX_feature_server("id")   # Returns reactive values (or NULL)
```

Modules with modals (Genome, Network, STRING) define a separate `mod_XX_feature_modal_ui(id)` function; the modal is rendered outside the `tabsetPanel` in `ui.R`.

**Namespacing rule for nested modules**: In the parent module's UI, use `ns("child_id")`; in the parent module's server, use the bare `"child_id"` (moduleServer handles namespacing automatically).

### Core Analysis Pipeline (all in global.R)

1. **`convertID(query, selectOrg)`** — Maps user gene list (any format) to Ensembl IDs via SQLite; auto-detects species.
2. **`geneInfo(converted, selectOrg)`** — Retrieves gene characteristics (CDS length, GC content, biotype).
3. **`FindOverlap(converted, gInfo, GO, selectOrg, ...)`** — Hypergeometric enrichment test across pathway database; returns FDR-corrected results.

### Database Layout

**`convertIDs.db`** — Central mapping DB: `orgInfo` (species metadata), `mapping` (ID↔Ensembl), `idIndex`, `quotes`.

**`db/Homo_sapiens.db`** (per-species) — `pathway` (gene↔pathwayID), `pathwayInfo` (description, URL), `geneInfo`, `categories`.

Default path is set by `IDEP_DATABASE` env var; current database version: `data113`.

### Key Global Constants (global.R)

```r
db_ver <- "data113"
STRING_DB_VERSION <- "12.0"
redundantGeneSetsRatio <- 0.95   # Remove pathways sharing ≥95% genes
minSetSize <- 2
maxSetSize <- 2000
mappingCoverage <- 0.60          # 60% gene mapping required for species auto-detection
```

## Development Notes
- Ensembl and RefSeq IDs should have version numbers stripped (e.g., `ENSG00000211459.2` → `ENSG00000211459`).
- `showTab`/`hideTab` in `server.R` target tab `value` attributes (integers), not namespaced DOM IDs — this works correctly even after modularization because `tabPanel(value = N)` is set inside each module's UI function.

## Known Issues (found during test development)

- **STRING ID hyperlinks** (`mod_07_genes.R:66`): Uses `grepl("ENS", id)` (unanchored) to assign Ensembl hyperlinks. STRING-format IDs like `"9606.ENSP00000269305"` match because `"ENS"` appears inside `"ENSP"`, producing a broken `ensembl.org/id/9606.ENSP...` URL. Fix: use `grepl("^ENS", id)`.
- **`cleanGeneSet` dedup order**: `unique()` runs before `remove_gene_version()`, so submitting `ENSG00000139618.2` and `ENSG00000139618.5` results in two copies of the same gene in the query set. Handled downstream by `convertID` in practice.
