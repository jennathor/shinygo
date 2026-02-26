####################################################
# Author: Steven Ge Xijin.Ge@sdstate.edu
# Lab: Ge Lab
# R version 4.0.5
# Project: ShinyGO v76
# File: server.R
# Purpose of file:main server logic of app
# Start data: NA (mm-dd-yyyy)
# Data last modified: 09-2-2021
#######################################################
server <- function(input, output, session) {
  options(warn = -1)

  # Welcome modal (hidden — uncomment to enable)
  # shiny::showModal(welcome_modal)

  # ============================================================================
  # MODULE SERVERS
  # ============================================================================

  # Sidebar module (ALL sidebar functionality)
  sidebar_values <- mod_01_sidebar_server(
    "sidebar",
    speciesChoice = speciesChoice,
    exampleGeneList = ExampleGeneList1,
    converted = converted,
    orgInfo = orgInfo
  )

  # Enrichment tab module
  enrichment_values <- mod_02_enrichment_server(
    "enrichment",
    sidebar_values = sidebar_values,
    quotes = quotes,
    min_gene_fold = min_gene_fold,
    redundantGeneSetsRatio = redundantGeneSetsRatio,
    converted = converted,
    geneInfoLookup = geneInfoLookup,
    geneInfoLookup_background = geneInfoLookup_background,
    converted_background = converted_background
  )

  # Chart tab module
  chart_values <- mod_03_chart_server(
    "chart",
    sidebar_values = sidebar_values,
    enrichment_values = enrichment_values
  )

  # Tree tab module
  mod_04_tree_server(
    "tree",
    sidebar_values = sidebar_values,
    enrichment_values = enrichment_values
  )

  # Network tab module
  mod_05_network_server(
    "network",
    sidebar_values = sidebar_values,
    enrichment_values = enrichment_values
  )

  # KEGG tab module
  mod_06_kegg_server(
    "kegg",
    sidebar_values = sidebar_values,
    enrichment_values = enrichment_values,
    converted = converted,
    keggSpeciesID = keggSpeciesID
  )

  # Genes tab module
  mod_07_genes_server(
    "genes",
    sidebar_values = sidebar_values,
    conversionTableData = conversionTableData
  )

  # Plots tab module
  mod_08_plots_server(
    "plots",
    sidebar_values = sidebar_values,
    geneInfoLookup = geneInfoLookup,
    geneInfoLookup_background = geneInfoLookup_background,
    converted_background = converted_background,
    ggplot2_theme = reactive(chart_values$ggplot2_theme())
  )

  # Genome tab module
  mod_09_genome_server(
    "genome",
    sidebar_values = sidebar_values,
    geneInfoLookup = geneInfoLookup,
    converted = converted,
    geneInfoLookup_background = geneInfoLookup_background,
    converted_background = converted_background
  )

  # STRING tab module
  mod_10_string_server(
    "string",
    sidebar_values = sidebar_values,
    conversionTableData = conversionTableData
  )

  # Groups tab module
  mod_13_groups_server(
    "groups",
    sidebar_values    = sidebar_values,
    enrichment_values = enrichment_values
  )

  # Promoter tab module (hidden — uncomment to enable)
  # mod_12_promoter_server(
  #   "promoter",
  #   sidebar_values = sidebar_values,
  #   converted = converted
  # )



  # ============================================================================
  # OBSERVERS
  # ============================================================================

  # When barplot_inside is selected, nudge sidebar maxTerms for a better default view
  observeEvent(chart_values$enrichChartType(), {
    req(chart_values$enrichChartType() == "barplot_inside")
    updateSliderInput(session, inputId = "chart-SortPathwaysPlotFontSize", value = 10)
    updateSliderInput(session, inputId = "chart-enrichChartAspectRatio", value = 1.5)
    updateSelectInput(session, inputId = "chart-SortPathwaysPlotLowColor", selected = "yellow")
    updateSelectInput(session, inputId = "sidebar-maxTerms", selected = "15")
    showNotification(
      "To improve the chart, we adjusted the aspect ratio, font size, and color:low. The number of pathways to show is set to 15.",
      type = "message",
      duration = 5
    )
  })

  # click_saved <- reactiveValues(GO = NULL)
  # observeEvent(eventExpr = sidebar_values$selectGO(), handlerExpr = { click_saved$GO <- sidebar_values$selectGO() })

  # Hide Groups tab, show only when GOBP is selected
  observeEvent(sidebar_values$selectGO(), {

    # Show KEGG tab only when KEGG is selected  #disabled as of 4/8/2022. Confused biologists.
    # if(sidebar_values$selectGO() == "KEGG") {
    #  showTab(inputId = "tabs", target = "2")
    # } else {
    # hideTab(inputId = "tabs", target = "2")
    # }

    # Show Groups tab only when GOBP is selected
    if (sidebar_values$selectGO() == "GOBP" | sidebar_values$selectGO() == "GOCC" | sidebar_values$selectGO() == "GOMF") {
      showTab(inputId = "tabs", target = "7")
    } else {
      hideTab(inputId = "tabs", target = "7")
    }
  })

  # Hide genome tab when STRINGdb is matched
  observe({
    showTab(inputId = "tabs", target = "8")

    if (sidebar_values$goButton() != 0 && !is.null(converted()$speciesMatched)) {
      if (grepl("STRING", converted()$speciesMatched[1, 1])) {
        hideTab(inputId = "tabs", target = "8")
      }
    }
  })

  # Species match message, stole from Gavin's code 4/20/22
  observe({
    req(
      isTRUE(sidebar_values$selectOrg() == speciesChoice[[1]]) # best matching species
      && !is.null(converted()) # finished
    )
    showNotification(
      ui = paste(
        gsub("\\(.*", "", converted()$speciesMatched[1, ]),
        ": is the best matching species. If that is incorrect,
                     please use the dropdown to select
                    your species."
      ),
      id = "species_match",
      duration = NULL,
      type = "error"
    )
  })



  # ============================================================================
  # REACTIVES
  # ============================================================================

  # Reactive objects accessessible to all modules:
  # Converted gene IDs (after clicking "GO" button)
  converted <- reactive({
    if (sidebar_values$goButton() == 0 || nchar(sidebar_values$input_text()) < 20) {
      return()
    }

    converted <- convertID(sidebar_values$input_text(), sidebar_values$selectOrg())

    # remove ensembl gene IDs mapped to the same gene (marked as duplicated in gene info)
    if(as.numeric(sidebar_values$selectOrg()) > 0) { # if it is ENSEMBL, not STRING species
      gene_info <- geneInfo(converted, sidebar_values$selectOrg())
      converted$IDs <- gene_info |>
        filter(!duplicated) |>
        filter(ensembl_gene_id %in% converted$IDs) |>
        pull(ensembl_gene_id)
      #conversionTable is not changed. Not unique.
    }
    converted

  })

  # Gene info lookup for converted gene IDs (after clicking "GO" button)
  geneInfoLookup <- reactive({
    if (sidebar_values$goButton() == 0) {
      return()
    }
    geneInfo(converted(), sidebar_values$selectOrg()) # uses converted gene ids thru converted() call
  })

  # Converted background gene IDs (after clicking "GO" button)
  converted_background <- reactive({
    if (sidebar_values$goButton() == 0 || is.null(sidebar_values$input_text_b())) {
      return()
    }
    if (nchar(sidebar_values$input_text_b()) < 10) {
      return()
    }
    
    converted <- convertID(sidebar_values$input_text_b(), sidebar_values$selectOrg())
    if(as.numeric(sidebar_values$selectOrg()) > 0) { # if it is ENSEMBL, not STRING species
      gene_info <- geneInfo(converted, sidebar_values$selectOrg())
      # remove ensembl gene IDs mapped to the same gene (marked as duplicated in gene info)
      converted$IDs <- gene_info |>
        filter(!duplicated) |>
        filter(ensembl_gene_id %in% converted$IDs) |>
        pull(ensembl_gene_id)
      #conversionTable is not changed. Not unique.
    }

    # if more than 100k genes, take samples
    if(length(converted$IDs) > maxGenesBackground + 1) {
      converted$IDs <- sample(converted$IDs, maxGenesBackground)
    }

    converted

  })

  # Gene info lookup for background gene IDs (after clicking "GO" button)
  geneInfoLookup_background <- reactive({
    if (sidebar_values$goButton() == 0 || nchar(sidebar_values$input_text_b()) < 10) {
      return()
    }
    if (is.null(converted_background())) {
      return()
    }
    geneInfo(converted_background(), sidebar_values$selectOrg()) # uses converted gene ids thru converted() call
  })

  # Conversion table data for Genes and STRING tabs (after clicking "GO" button)
  conversionTableData <- reactive({
    if (sidebar_values$goButton() == 0) {
      return()
    } # still have problems when geneInfo is not found!!!!!
    tem <- sidebar_values$selectGO()
    tem <- sidebar_values$selectOrg()
    tem <- sidebar_values$minFDR()
    isolate({
      withProgress(message = sample(quotes, 1), detail = "Looking up gene Info", {
        tem <- converted()
        incProgress(0.1)
        tem2 <- geneInfoLookup()
        incProgress(0.3)
        incProgress(0.6)
        if (is.null(tem)) {
          as.data.frame("ID not recognized.")
        } else {
          # some STRINGdb species has geneInfo, alought incomplete.
          if (dim(tem2)[1] <= 1 | grepl("STRINGdb", converted()$species$name2)) {
            merged <- tem$conversionTable
            ix <- which(colnames(merged) == "ensembl_gene_id")
            colnames(merged)[ix] <- "STRINGdb ID"
          } else { # if gene info is  available
            #         if('chromosome_name' %in% colnames(tem2)) {
            merged <- merge(tem$conversionTable, tem2, by = "ensembl_gene_id")

            desired_cols <- c(
              "User_input", "symbol", "ensembl_gene_id", "entrezgene_id",
              "gene_biotype", "Species", "chromosome_name", "start_position",
              "description", "percentage_gc_content", "transcript_count",
              "genomeSpan", "cds_length", "transcript_length", "FiveUTR",
              "ThreeUTR", "nExons"
            )
            merged <- merged[, desired_cols[desired_cols %in% colnames(merged)], drop = FALSE]

            tem3 <- as.data.frame(tem$originalIDs)
            colnames(tem3) <- "User_input"
            merged <- merge(merged, tem3, all = T)
            merged$ensembl_gene_id[which(is.na(merged$ensembl_gene_id))] <- "Not mapped"
            chrName <- suppressWarnings(as.numeric(as.character(merged$chromosome_name)))
            merged <- merged[order(
              merged$gene_biotype,
              chrName,
              merged$start_position
            ), ]
            merged$start_position <- merged$start_position / 1e6
            colnames(merged)[1:9] <- c(
              "Pasted", "Symbol", "Ensembl Gene ID", "Entrez",
              "Type", "Species", "Chr", "Position (Mbp)", "Description"
            )
          }
        }
        incProgress(0.9)
        return(merged)
      })
    }) # avoid showing things initially
  })



  # ============================================================================
  # MISCELLANEOUS OUTPUTS
  # ============================================================================

  # Welcome modal
  welcome_modal <- shiny::modalDialog(
    title = "Find ShinyGO helpful? Send us an email today so it will be here next year.",
    tags$h3("We are still working on it until June 5th!"),
    tags$p("We need your help to support our NIH grant proposal due June 5th."),
 
    tags$p(
      " Please take a few minutes to send us an email today:  ",
      a(
        " gelabinfo@gmail.com",
        href = "mailto:gelabinfo@gmail.com",
        target = "_blank"
      ), 
      "  Thank you!"
    ),

    easyClose = FALSE,
    size = "l"
  )

  # Gene ID examples modal
  output$showGeneIDs4Species <- renderTable(
    {
      if (input$userSpeciesIDexample == 0) {
        return()
      }
      withProgress(message = "Retrieving gene IDs (2 minutes)", {
        geneIDs <- showGeneIDs(species = input$userSpeciesIDexample, nGenes = 10)
        incProgress(1, detail = paste("Done"))
      })
      geneIDs
    },
    digits = -1,
    spacing = "s",
    striped = TRUE,
    bordered = TRUE,
    width = "auto",
    hover = T
  )


  ### Currently unused components (species table, orgInfo table, and tableDetail)
  output$species <- renderTable(
    {
      if (sidebar_values$goButton() == 0) {
        return()
      }
      tem <- sidebar_values$selectGO()
      tem <- sidebar_values$selectOrg()
      tem <- sidebar_values$minFDR()
      isolate({ # tem <- convertID(sidebar_values$input_text(),sidebar_values$selectOrg() );
        withProgress(message = "Converting gene IDs", {
          tem <- converted()
          incProgress(1, detail = paste("Done"))
        })

        if (is.null(tem)) {
          as.data.frame("ID not recognized.")
        } else {
          tem$speciesMatched
        }
      }) # avoid showing things initially
    },
    digits = -1,
    spacing = "s",
    striped = TRUE,
    bordered = TRUE,
    width = "auto",
    hover = T
  )

  output$orgInfoTable <- DT::renderDataTable({
    df <- orgInfo[, c("ensembl_dataset", "name", "totalGenes")]
    colnames(df) <- c("Ensembl/STRING-db ID", "Name (Assembly)", "Total Genes")
    row.names(df) <- NULL
    df
  })

  output$tableDetail <- renderTable(
    {
      if (sidebar_values$goButton() == 0) {
        return()
      }

      tem <- enrichment_values$significantOverlaps()
      tem$x
    },
    digits = -1,
    spacing = "s",
    striped = TRUE,
    bordered = TRUE,
    width = "auto",
    hover = T
  )

}
