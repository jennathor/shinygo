#' Find taxon ID by species ID
#'
#' Find a species in the iDEP database with an ID.
#'
#' @param species_id Species ID to search the database with
#' @param org_info iDEP data org_info file
#'
#' @return Taxon ID for the species
find_taxon_by_id <- function(species_id, org_info) {
  # find species name use id; column name differs between production and local dev
  taxon_col <- if ("taxon_id" %in% colnames(org_info)) "taxon_id" else "taxon"
  return(org_info[which(org_info$id == species_id), taxon_col])
}


#' 10_string UI Function
#'
#' @description A shiny Module for the STRING tab, showing STRING-db based
#'   mapping stats, functional enrichment, and PPI network access.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_10_string_ui <- function(id) {
  ns <- shiny::NS(id)

  tabPanel("STRING",
    value = 11,

    textOutput(ns("STRINGDB_mapping_stat")),
    tags$head(tags$style(paste0("#", ns("STRINGDB_mapping_stat"), "{color: blue;font-size: 15px;}"))),
    br(),
    actionButton(ns("ModalPPI"), "PPI network of DEGs"), br(), br(),
    selectInput(ns("STRINGdbGO"),
      label = "Functional Enrichment",
      choices = list(
        "GO Biological Process" = "Process",
        "GO Cellular Component" = "Component",
        "GO Molecular Function" = "Function",
        "KEGG" = "KEGG",
        "Pfam" = "Pfam",
        "InterPro" = "InterPro"
      ),
      selected = "Process"
    ),
    downloadButton(ns("STRING_enrichmentDownload")),
    tableOutput(ns("stringDB_GO_enrichment")),
    br(), br(),
    h5(
      "To validate your results independent of our algorithm and database,
      your genes are sent to STRING-db website for enrichment analysis.
      This also enables the
       retrieval of a protein-protein network. If it is running,
       please wait until it finishes. The second time it is faster."
    ),
  )
}

#' STRING PPI modal UI — placed outside tabsetPanel in mainPanel
mod_10_string_modal_ui <- function(id) {
  ns <- shiny::NS(id)

  bsModal(ns("ModalExamplePPI"), "Protein-protein interaction networks ", ns("ModalPPI"),
    size = "large",
    h5("By sending your genes to the STRING website,
		shinyGO is retrieving a sub-network, calculating PPI enrichment,
	  and generating custom URLs to the STRING website containing your genes. This can take 5 minutes. Patience will pay off! "),
    sliderInput(ns("nGenesPPI"), label = h5("Genes to include:"), min = 0, max = 400, value = 20, step = 10),
    # ,htmlOutput(ns("stringDB_network_link"))
    # ,tags$head(tags$style(paste0("#", ns("stringDB_network_link"), "{color: blue; font-size: 15px;}")))

    plotOutput(ns("stringDB_network1"))
  ) # bsModal
}


#' 10_string Server Functions
#'
#' @param id                  Module ID string
#' @param sidebar_values      Reactive list returned by mod_01_sidebar_server
#' @param conversionTableData Reactive returning the merged gene info data frame.
#'   Defined in server.R (shared with Genes tab) and passed in here.
#'
#' @noRd
#'
#' @importFrom shiny moduleServer
mod_10_string_server <- function(id, sidebar_values, conversionTableData) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # find Taxonomy ID from species official name
    findTaxonomyID <- reactive({
      if (sidebar_values$goButton() == 0) {
        return(NULL)
      }

      find_taxon_by_id(sidebar_values$selectOrg(), orgInfo)

    })


    STRINGdb_geneList <- reactive({
      if (sidebar_values$goButton() == 0) {
        return(NULL)
      }
      library(STRINGdb, verbose = FALSE)
      tem <- sidebar_values$selectOrg()

      ####################################

      if (is.null(conversionTableData())) {
        return(NULL)
      } # this has to be outside of isolate() !!!
      # if(sidebar_values$selectOrg() == "NEW" && is.null( input$gmtFile) ) return(NULL) # new but without gmtFile
      NoSig <- as.data.frame("No significant enrichment found.")
      taxonomyID <- findTaxonomyID()

      if (is.null(taxonomyID)) {
        return(NULL)
      }

      isolate({
        withProgress(message = sample(quotes, 1), detail = "Mapping gene ids (5 minutes)", {

          # Intialization
          string_db <- STRINGdb$new(
            version = STRING_DB_VERSION, species = taxonomyID,
            score_threshold = 0, input_directory = ""
          )

          # using expression data
          genes <- conversionTableData()
          # STRINGdb species the columns are ensemble_gene_id
          ix <- which(colnames(genes) == "Ensembl Gene ID" | colnames(genes) == "STRINGdb ID")
          colnames(genes)[ix] <- c("gene")
          genes$lfc <- 1
          # remove space character in front of gene symbols. Otherwise STRING won't convert
          genes$gene <- gsub(" ", "", genes$gene)
          mapped <- string_db$map(genes, "gene", removeUnmappedRows = TRUE)

          incProgress(1 / 4, detail = paste("up regulated"))
          up <- subset(mapped, lfc > 0, select = "STRING_id", drop = TRUE)

          incProgress(1 / 2, detail = "Down regulated")
          down <- subset(mapped, lfc < 0, select = "STRING_id", drop = TRUE)

          mappingRatio <- nrow(mapped) / nrow(genes)
          if (nrow(mapped) == 0) {
            return(NULL)
          } else {
            return(list(up = up, down = down, ratio = mappingRatio, geneTable = mapped))
          }
          incProgress(1)
        }) # progress
      }) # isolate
    })

    output$STRINGDB_mapping_stat <- renderText({
      if (sidebar_values$goButton() == 0) {
        return(NULL)
      }

      if (is.null(STRINGdb_geneList())) {
        return("No genes mapped by STRINGdb. Please enter or double-check species name above.")
      }
      if (!is.null(STRINGdb_geneList())) {
        tem <- paste0(100 * round(STRINGdb_geneList()$ratio, 3), "% genes mapped by STRING web server.")
        if (STRINGdb_geneList()$ratio < 0.3) tem <- paste(tem, "Warning!!! Very few gene mapped. Double check if the correct species is selected.")
        return(tem)
      }
    })

    stringDB_GO_enrichmentData <- reactive({
      if (sidebar_values$goButton() == 0) {
        return(-1)
      }
      taxonomyID <- findTaxonomyID()
      if (is.null(taxonomyID)) {
        return(NULL)
      }
      library(STRINGdb, verbose = FALSE)
      withProgress(message = sample(quotes, 1), detail = "Enrichment analysis", {
        tem <- input$STRINGdbGO
        # Intialization
        string_db <- STRINGdb$new(
          version = STRING_DB_VERSION, species = taxonomyID,
          score_threshold = 0, input_directory = ""
        )

        # using expression data
        genes <- conversionTableData()
        minGenesEnrichment <- 1
        if (is.null(genes)) {
          return(-2)
        } else if (dim(genes)[1] <= minGenesEnrichment) {
          return(-2) # if has only few genes
        } else {
          # GO
          ids <- STRINGdb_geneList()$up
          if (length(ids) <= minGenesEnrichment || is.null(ids)) {
            return(-2)
          }
          incProgress(1 / 3)
          result <- string_db$get_enrichment(ids, category = input$STRINGdbGO, methodMT = "fdr", iea = TRUE)
          if (nrow(result) == 0 || is.null(result)) {
            return(-2)
          } else {
            if (min(result$fdr) > sidebar_values$minFDR()) {
              return(-2)
            } else {
              result <- result[which(result$fdr < sidebar_values$minFDR()), ]
              incProgress(1, detail = paste("Done"))
              return(result)
            } # end of check minFDR
          } # check results
        } # end of check genes if
      }) # progress
    }) # end of stringDB_GO_enrichmentData

    output$stringDB_GO_enrichment <- renderTable(
      {
        result <- stringDB_GO_enrichmentData()

        req(!is.null(result))
        if(class(result) == "numeric") {
          if (result == -1) {
            return(NULL)
          } else if (result == -2) {
            return(as.data.frame("No significant enrichment found."))
          }
        } else {
          result <- dplyr::select(
            result,
            c(
              "fdr", "number_of_genes", "term",
              "description"
            )
          )
          colnames(result) <- c(
            "FDR", "nGenes", "GO terms or pathways",
            "Description"
          )
          result$FDR <- as.character(result$FDR)
          if (nrow(result) > 30) {
            result <- result[1:30, ]
          }
          return(result)
        } # end of if else
      },
      digits = 4,
      spacing = "s",
      include.rownames = F,
      striped = TRUE,
      bordered = TRUE,
      width = "auto",
      hover = T
    ) # renderTable

    output$STRING_enrichmentDownload <- downloadHandler(
      filename = function() {
        paste0("STRING_enrichment", input$STRINGdbGO, ".csv")
      },
      content = function(file) {
        write.csv(stringDB_GO_enrichmentData(), file)
      }
    ) # downloadHandler

    output$stringDB_network1 <- renderPlot(
      {
        library(STRINGdb)
        if (sidebar_values$goButton() == 0) {
          return(NULL)
        }


        tem <- input$STRINGdbGO
        tem <- input$nGenesPPI
        taxonomyID <- findTaxonomyID()
        if (is.null(taxonomyID)) {
          return(NULL)
        }
        ####################################

        if (is.null(STRINGdb_geneList())) {
          return(NULL)
        }

        isolate({
          withProgress(message = sample(quotes, 1), detail = "Enrichment analysis", {
            # Intialization
            string_db <- STRINGdb$new(
              version = STRING_DB_VERSION, species = taxonomyID,
              score_threshold = 0, input_directory = ""
            )
            # only up regulated is ploted
            ngenes1 <- input$nGenesPPI
            if (ngenes1 < 2) ngenes1 <- 2
            for (i in c(1:1)) {
              incProgress(1 / 2, detail = paste("Plotting network"))


              ids <- STRINGdb_geneList()[[i]]
              if (length(ids) > ngenes1) { # n of genes cannot be more than 400
                ids <- ids[1:ngenes1]
              }
              incProgress(1 / 3)
              string_db$plot_network(ids, add_link = FALSE)
            }
          }) # progress
        }) # isolate
      },
      width = 1000,
      height = 600
    )

    output$stringDB_network_link <- renderUI({
      library(STRINGdb, verbose = FALSE)

      tem <- input$STRINGdbGO
      tem <- input$nGenesPPI
      taxonomyID <- findTaxonomyID()
      if (is.null(taxonomyID)) {
        return(NULL)
      }

      ####################################
      if (is.null(STRINGdb_geneList())) {
        return(NULL)
      }

      isolate({
        withProgress(message = sample(quotes, 1), detail = "PPI Enrichment and link", {
          # Intialization
          string_db <- STRINGdb$new(
            version = STRING_DB_VERSION, species = taxonomyID,
            score_threshold = 0, input_directory = ""
          )
          # upregulated
          ids <- STRINGdb_geneList()[[1]]

          ngenes1 <- input$nGenesPPI
          if (ngenes1 < 2) ngenes1 <- 2

          if (length(ids) > ngenes1) { # n of genes cannot be more than 400
            ids <- ids[1:ngenes1]
          }
          incProgress(1 / 4)
          link1 <- string_db$get_link(ids)


          tem <- paste("<a href=\"", link1, "\" target=\"_blank\"> Click here for an interactive and annotated network </a>")
          return(HTML(tem))

          incProgress(1)
        }) # progress
      }) # isolate
    })

  })
}
