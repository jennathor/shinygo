#' 02_enrichment UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_02_enrichment_ui <- function(id, sidebar_id = "sidebar") {
  ns <- shiny::NS(id)
  sidebar_ns <- shiny::NS(sidebar_id)

  tabPanel("Enrichment",
    value = 1,
    conditionalPanel(
      sprintf("input['%s'] == 0", sidebar_ns("goButton")), # welcome screen
      br(),
      fluidRow(
        column(
          width = 9,
          h4("ShinyGO: a graphical gene-set enrichment tool for animals and plants")
        ),
        column(
          width = 3,
          img(
            src = "shinygo_logo.png",
            width = "43",
            height = "50"
          )
        )
      ),
      p("9/5/25: v.0.85. Database updated to Ensembl Release 113 and STRING-db v12. "),
      p("You can still use the old versions using links on the About tab.", 
        "To support this effort, please cite our paper like ",
          a("these 4000+ papers.", href = "https://scholar.google.com/scholar?oi=bibs&hl=en&cites=4205886424733220184&as_sdt=5"),
          "Just including URL is not enough.",
        a("Email Jenny ", href = "mailto:gelabinfo@gmail.com?Subject=ShinyGO"),
        "(gelabinfo@gmail.com) for questions, suggestions or data contributions.",
        "Follow Dr Ge on ", a("Twitter", href = "https://twitter.com/StevenXGe"), " and ",
        a("LinkedIn", href = "https://www.linkedin.com/in/steven-ge-ab016947/", target = "_blank"),
        " for updates. ",
        "To request to add a new species/genome, fill in this ",
        a("Form.", href = "https://forms.gle/zLtLnqxkW187AgT76"), 
        "We will try to accommodate commonly requested genomes. "
      ),            
      p("For-profit organizations: contact us for local installation or customization services."),
      p("Under active delovelopment with support from NIH. Report bugs or request features on our ", 
          a("GitHub repository.", href = "https://github.com/gexijin/shinygo", target = "_blank")),
      p("NO WARRANTY. Please verify results using other tools.
      Enrichment results may vary depending on gene ID mapping, data sources, 
      database versions, and methods (particularily ranking).
        "),
      h3("GO Enrichment analysis, plus a lot more!"),
      p("Just paste your gene list to get enriched GO terms and othe pathways for over 14,000 species,
      based on annotation from Ensembl and STRING-db."),
      br(), img(src = "enrich.png", align = "center", width = "660", height = "339"),
      br(), img(src = "enrichmentChart.png", align = "center", width = "700", height = "400"),
      br(), br(), img(src = "KEGG2.png", align = "center", width = "541", height = "360"),
      br(), br(), img(src = "GOtree3.png", align = "center", width = "500", height = "258"),
      br(), br(), img(src = "GOnetwork2.png", align = "center", width = "500", height = "248"),
      br(), br(), img(src = "PPInetwork2.png", align = "center", width = "500", height = "391"),
      br(), br(), img(src = "chr.png", align = "center", width = "444", height = "338"),
      br(), br(), img(src = "downSyndrome.png", align = "center", width = "371", height = "276"),
      #br(), br(), img(src = "promoter.png", align = "center", width = "717", height = "288")
    ),
    br(),
    conditionalPanel(
      sprintf("input['%s'] != 0", sidebar_ns("goButton")),
      div(
        style = "display:inline-block",
        selectInput(
          inputId = ns("SortPathways"),
          label = NULL,
          choices = c(
            "Sort by FDR" = "Sort by FDR",
            "Sort by Fold Enrichment" = "Sort by Fold Enrichment",
            "Sort by average ranks(FDR & Fold)" = "Sort by FDR & Fold Enrichment",
            "Select by FDR, sort by Fold Enrichment" = "Select by FDR, sort by Fold Enrichment",
            "Sort by Genes" = "Sort by Genes",
            "Sort by Category Name" = "Sort by Category Name"
          ),
          selected = "Select by FDR, sort by Fold Enrichment"
        ),
        style = "algn:right"
      )
    ),
    tableOutput(ns("EnrichmentTable")),
    conditionalPanel(
      sprintf("input['%s'] != 0", sidebar_ns("goButton")),
      downloadButton(ns("downloadEnrichment"), "Top Pathways shown above"),
      downloadButton(ns("downloadEnrichmentAll"), "Results on all Pathways"),
      br(), br(),
      h3("Methods"),
      p("All query genes are converted to ENSEMBL gene IDs or STRING-db protein IDs,
        as our gene ID mapping and pathway data primarily come from these sources. 
        For model organisms, we manually compile extensive pathway lists from various
        public databases."),

      p("P-values are calculated using the hypergeometric test, and false discovery 
        rates (FDRs) are computed via the Benjamini-Hochberg method to correct for 
        multiple testing. Fold enrichment is defined as the percentage of genes in 
        your list that are in a pathway divided by the corresponding percentage in 
        the background genes. While FDR measures statistical significance, fold 
        enrichment indicates effect size."),

      p("We recommend that users provide their own list of background genes, which 
        could include all genes detected in an experiment, such as genes with probes 
        on a DNA microarray, passed a minimal filter in RNA-seq analysis, or detected 
        in a proteomics experiment. If no background genes are uploaded, the default
        is to use all protein-coding genes. Alternatively, you may select the option
        'Use pathway DB for gene counts,' which calculates the background based on
        the total number of unique genes in the chosen pathway database, limited
        between 5,000 and 30,000 genes. When this option is selected, any genes in 
        the user's original list that are not in the pathway database are excluded."),

      p("Only pathways within specified size limits, as defined by the 
      'Pathway size: (Min, Max)' settings, are considered. Results for smaller
      pathways can be noisy, but some pathways or GO terms have only a few genes. 
      After analysis, pathways are filtered by a user-defined FDR cutoff. Significant 
      pathways are then sorted in different ways, and only the top-ranked are shown 
      in the table above.

      By default, 'Select by FDR, then by Fold Enrichment' is used, where pathways
      are first filtered and sorted by FDR, and then the top 20 are sorted by fold 
      enrichment. In other words, the default setting shows the top 20 most 
      significant pathways ranked by fold enrichment. When the 'Sort by average 
      ranks (FDR & fold enrichment)' option is selected, pathways are sorted by 
      the average of their ranks based on both FDR and fold enrichment. When 
      'Sort by FDR' is selected, pathways are ranked by FDR and only the top 20 are shown.

      The 'Remove redundancy' option eliminates similar pathways that share 95% of 
        their genes and 50% of the words in their names, representing them with the
        pathway that has the highest significance."),

      h3("Interpreting GO Enrichment Results"),

      p("The Gene Ontology (GO) includes tens of thousands of terms (functional 
        categories), each tested individually for enrichment. Hundreds or even 
        thousands of GO terms can be statistically significant. These terms are 
        filtered, ranked, and only the top ones are displayed. Understanding this
        process is crucial for interpreting GO enrichment results."),

      tags$ul(
        tags$li("P-value: Reflects the statistical significance of the enrichment.
                Lower values suggest a lower likelihood of the result occurring by
                chance under the null hypothesis. FDR q-values adjust P-values for
                multiple testing to control the proportion of type I errors."),
        tags$li("Fold Enrichment: Measures the magnitude of enrichment. Higher values
                indicate stronger enrichment and are an important metric of effect size."),
        tags$li("Pathway Genes: The total number of genes in a pathway or GO term."),
        tags$li("nGenes: The number of genes in the pathway that overlap with your gene list.")
      ),

      tags$p("Exercise caution when interpreting FDR values of 0.01 or 0.001 for GO 
            terms, as these levels often represent noise due to the vast number of 
            terms tested. For a gene list of reasonable size, more significant results 
            (FDR < 1E-5) are expected."),

      tags$p("Large pathways, such as the cell cycle, often show smaller FDRs due to 
            increased statistical power, while smaller pathways might have higher 
            FDRs despite their biological relevance. Enrichment analysis tends to 
            favor larger pathways."),

      tags$p("With a default cutoff of FDR < 0.05, thousands of significant GO terms
            may be detected, though only a subset is shown. Therefore, the method of
            filtering and ranking these terms is crucial."),

      tags$p("With large sample sizes, small differences can appear extremely 
            significant. In addition to FDR, fold enrichment should also be 
            considered when prioritizing pathways, as it reflects the strength 
            of the enrichment. We offer several methods that consider both FDR 
            q-values and fold enrichment."),

      tags$p("Many GO terms are closely related (e.g., 'Cell Cycle', 'Regulation of 
            Cell Cycle') and can dominate the top 20, obscuring other pathways. To 
            avoid this, consider examining the top 50 terms. Additionally, use tree
            plots and network plots to identify clusters of related GO terms and 
            uncover overarching themes."),

      tags$p("Discuss the most significant pathways first, even if they do not fit 
            your initial expectations."),
    )
  )
}

#' 02_enrichment Server Functions
#'
#' @noRd
#'
#' @importFrom shiny moduleServer
mod_02_enrichment_server <- function(id, sidebar_values, quotes, min_gene_fold,
                                      redundantGeneSetsRatio, converted, geneInfoLookup,
                                      geneInfoLookup_background, converted_background) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    significantOverlapsAll <- reactive({
      if (sidebar_values$goButton() == 0 || is.null(sidebar_values$selectGO()) || nchar(sidebar_values$input_text()) < 20) {
        return()
      }
      tem <- sidebar_values$selectOrg()
      tem <- sidebar_values$selectGO()
      tem <- sidebar_values$gene_count_pathwaydb()
      tem <- sidebar_values$minSetSize()
      tem <- sidebar_values$maxSetSize()

      isolate({
        withProgress(message = sample(quotes, 1), detail = "enrichment analysis", {
          # gene info is passed to enable lookup of gene symbols
          tem <- geneInfoLookup()
          tem <- tem[which(tem$Set == "List"), ]
          temb <- geneInfoLookup_background()
          if (class(temb) == "data.frame") {
            temb <- temb[which(temb$Set == "List"), ]
          }
          enrichment <- FindOverlap(converted(), tem, sidebar_values$selectGO(), sidebar_values$selectOrg(),
            converted_background(), temb,
            minSetSize = sidebar_values$minSetSize(), maxSetSize = sidebar_values$maxSetSize(),
            gene_count_pathwaydb = sidebar_values$gene_count_pathwaydb()
          )
          return(enrichment)
        })
      })
    })

    observe({
      req(!is.null(significantOverlapsAll() )) # stop if null
      req(sidebar_values$goButton() != 0)
      req(significantOverlapsAll()$x[1,1] == "ID not recognized!" )

      shiny::showModal(
        shiny::modalDialog(
          size = "s",
          p("None of the gene IDs mapped to the IDs of the selected species.
            From ShinyGO 0.80, you have to select the correct species first.
            If you do not select, it defaults to human.
            "),
          easyClose = TRUE
        )
      )

    })

    # Filtering and ranking pathways
    significantOverlaps <- reactive({
      if (sidebar_values$goButton() == 0 || is.null(sidebar_values$selectGO()) || nchar(sidebar_values$input_text()) < 20) {
        return()
      }
      if (is.null(significantOverlapsAll())) {
        return(NULL)
      }

      enrichment <- significantOverlapsAll()
      withProgress(message = sample(quotes, 1), detail = "Sorting and filtering pathways", {
        if (dim(enrichment$x)[2] > 1) { # when there is no overlap, returns a data frame with 1 row and 1 column

          # filter by FDR-------------------------------------------------------------
          enrichment$x <- enrichment$x[enrichment$x[, 1] < sidebar_values$minFDR(), ]

          incProgress(0.1)
          # Sort and keep top pathways -------------------------------------------------------
          if (input$SortPathways == "Select by FDR, sort by Fold Enrichment") {
            # sort by FDR
            enrichment$x <- enrichment$x[order(enrichment$x[, 1]), ]
            # filter/top 20
            if (dim(enrichment$x)[1] > as.integer(sidebar_values$maxTerms())) {
              enrichment$x <- enrichment$x[1:as.integer(sidebar_values$maxTerms()), ]
            }
            # rank by fold
            enrichment$x <- enrichment$x[order(enrichment$x[, 4], decreasing = TRUE), ]
          } else {
            if (input$SortPathways == "Sort by FDR") {
              enrichment$x <- enrichment$x[order(enrichment$x[, 1]), ]
            }
            if (input$SortPathways == "Sort by Fold Enrichment") {
              enrichment$x <- enrichment$x[order(enrichment$x[, 4], decreasing = TRUE), ]
              # when sorting by fold, sometimes tiny pathways on top. Here we require at
              # least 10 genes
              enrichment$x <- enrichment$x[which(enrichment$x[, 3] > min_gene_fold), ]
            }
            if (input$SortPathways == "Sort by Genes") {
              enrichment$x <- enrichment$x[order(enrichment$x[, 2], decreasing = TRUE), ]
            }
            if (input$SortPathways == "Sort by Category Name") {
              enrichment$x <- enrichment$x[order(enrichment$x[, 5]), ]
            }
            if (input$SortPathways == "Sort by FDR & Fold Enrichment") {
              fdr_rank <- rank(enrichment$x[, 1]) # rank by FDR
              fold_rank <- rank(-1 * enrichment$x[, 4]) # rank by fold_enrichment, descending
              average_rank <- (fdr_rank + fold_rank) / 2
              enrichment$x <- enrichment$x[order(average_rank), ]
            }
          }
          incProgress(0.3)

          # preliminary filtering to save time on string manipulations
          if (dim(enrichment$x)[1] > 3 * as.integer(sidebar_values$maxTerms())) {
            enrichment$x <- enrichment$x[1:(3 * as.integer(sidebar_values$maxTerms())), ]
          }


          # remove redundant gene sets-------------------------------------------
          if (sidebar_values$removeRedundantSets()) reduced <- redundantGeneSetsRatio else reduced <- FALSE
          incProgress(0.2)
          # reduced=FALSE no filtering,  reduced = 0.9 filter sets overlap with 90%
          if (reduced != FALSE && dim(enrichment$x)[1] > 5) {
            n <- nrow(enrichment$x)
            flag1 <- rep(TRUE, n)
            # note that it has to be two space characters for splitting
            geneLists <- lapply(
              enrichment$x$Genes,
              function(y) unlist(strsplit(as.character(y), " |  |   "))
            )
            pathways <- lapply(
              enrichment$x$Pathway,
              function(y) unlist(strsplit(as.character(y), " |  |   "))
            )
            for (i in 2:n) {
              for (j in 1:(i - 1)) {
                if (flag1[j]) { # skip if this one is already removed
                  ratio1 <- length(intersect(geneLists[[i]], geneLists[[j]])) /
                    length(union(geneLists[[i]], geneLists[[j]]))

                  # if sufficient genes overlap
                  if (ratio1 > reduced) {
                    # are pathway names similar
                    ratio2 <- length(intersect(pathways[[i]], pathways[[j]])) /
                      length(union(pathways[[i]], pathways[[j]]))
                    # if 50% of the words in the pathway name shared
                    if (ratio2 > 0.5) {
                      flag1[i] <- FALSE
                    }
                  }
                }
              }
            }
            # remove similar pathways
            enrichment$x <- enrichment$x[which(flag1), ]
          }
          incProgress(0.9)

          # keep top pathways
          if (dim(enrichment$x)[1] > as.integer(sidebar_values$maxTerms())) {
            enrichment$x <- enrichment$x[1:as.integer(sidebar_values$maxTerms()), ]
          }

          if (sidebar_values$abbreviatePathway()) {
            enrichment$x[, 5] <- gsub("Positive regulation", "Pos. reg.", enrichment$x[, 5])
            enrichment$x[, 5] <- gsub("Negative regulation", "Neg. reg.", enrichment$x[, 5])
            enrichment$x[, 5] <- gsub("Regulation", "Reg.", enrichment$x[, 5])
            enrichment$x[, 5] <- gsub(" regulation ", " reg. ", enrichment$x[, 5])
            enrichment$x[, 5] <- gsub(" process ", " proc. ", enrichment$x[, 5])
            enrichment$x[, 5] <- substr(enrichment$x[, 5], 1, 100) # maximum 80 characters
          }
        }
      }) # progress bar

      return(enrichment)
    })

    output$EnrichmentTable <- renderTable(
      {
        if (sidebar_values$goButton() == 0) {
          return(NULL)
        }
        tem <- sidebar_values$input_text_b() # just to make it re-calculate if user changes background

        myMessage <- "Analyzing genes."

        if (is.null(significantOverlaps())) {
          return(NULL)
        }
        # this solves an error when there is no significant enrichment

        if (ncol(significantOverlaps()$x) == 1) {
          return(significantOverlaps()$x)
        }

        withProgress(message = sample(quotes, 1), detail = myMessage, {
          pathways <- significantOverlaps()$x
          # remove pathway ID  only in Ensembl species
          if (!sidebar_values$show_pathway_id() && as.integer(sidebar_values$selectOrg()) > 0) {
            pathways$Pathway <- remove_pathway_id(pathways$Pathway, sidebar_values$selectGO())
          }
          pathways$Pathway <- hyperText(pathways$Pathway, pathways$URL)

          pathways <- pathways[, -7]
          # rownames(pathways) <- NULL
          # pathways[, 1] <- as.numeric( format(pathways[, 1], scientific = TRUE, digits = 3 ) )
          pathways[, 4] <- as.character(round(pathways[, 4], 1))
          pathways[, 2] <- as.character(pathways[, 2]) # convert total genes into character 10/21/19
          pathways[, 3] <- as.character(pathways[, 3]) # convert total genes into character 10/21/19
          colnames(pathways)[5] <- "Pathways (click for details)"

          incProgress(1, detail = paste("Done"))
        })

        if (dim(pathways)[2] > 1) pathways[, 2] <- as.character(pathways[, 2])

        if (dim(pathways)[2] == 1) {
          return(pathways)
        } else {
          return(pathways[, 1:5])
        } # If no significant enrichment found x only has 1 column.
      },
      digits = -1,
      spacing = "s",
      striped = TRUE,
      bordered = TRUE,
      width = "auto",
      hover = TRUE,
      sanitize.text.function = function(x) x
    )

    output$downloadEnrichment <- downloadHandler(
      filename = function() {
        "enrichment.csv"
      },
      content = function(file) {
        write.csv(significantOverlaps()$x, file, row.names = FALSE)
      }
    )

    output$downloadEnrichmentAll <- downloadHandler(
      filename = function() {
        "enrichment_all.csv"
      },
      content = function(file) {
        write.csv(significantOverlapsAll()$x, file, row.names = FALSE)
      }
    )

    list(
      SortPathways = reactive(input$SortPathways),
      significantOverlaps = significantOverlaps,
      significantOverlapsAll = significantOverlapsAll
    )
  })
}