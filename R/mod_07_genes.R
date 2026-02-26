#' 07_genes UI Function
#'
#' @description A shiny Module for the Genes tab, which displays the gene
#'   ID conversion table with optional detailed descriptions and a download
#'   button.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_07_genes_ui <- function(id) {
  ns <- shiny::NS(id)
  tabPanel("Genes",
    value = 6,
    fluidRow(
      column(3, downloadButton(ns("downloadGeneInfo"), "More info")),
      column(4, checkboxInput(ns("showDetailedGeneInfo"), "Detailed Description", value = FALSE))
    ),
    tableOutput(ns("conversionTable"))
  )
}

#' 07_genes Server Functions
#'
#' @param id               Module ID string
#' @param sidebar_values   Reactive list returned by mod_01_sidebar_server
#' @param conversionTableData  Reactive returning the merged gene info data
#'   frame. Defined in server.R (shared with STRING tab) and passed in here.
#'
#' @noRd
#'
#' @importFrom shiny moduleServer
mod_07_genes_server <- function(id, sidebar_values, conversionTableData) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Gene conversion table ----------------------------------------------------
    output$conversionTable <- renderTable(
      {
        if (sidebar_values$goButton() == 0) {
          return()
        } # still have problems when geneInfo is not found!!!!!
        req(sidebar_values$selectOrg())
        tem <- input$showDetailedGeneInfo # reactive dependency — outside isolate intentionally
        isolate({
          if (dim(conversionTableData())[2] < 9) { # STRINGdb species: only 3 columns
            df <- conversionTableData()
          } else { # ENSEMBL species
            df <- conversionTableData()[, 1:9]
            # show detailed gene info for string species
            if (!input$showDetailedGeneInfo) {
              df$Description <- gsub(";.*|\\[.*", "", df$Description)
            }
            # Remove columns with all missing values; chr and start position in STRINGdb species
            df <- df[, which(!apply(is.na(df), 2, sum) == nrow(df))]
            df$Species <- gsub("STRINGdb", "", df$Species)

            df$Type <- gsub(".*pseudogene", "pseudo", df$Type)
            # coding is not shown
            df$Type <- gsub("protein_coding", "coding", df$Type)
            df$Type <- gsub("_gene", "", df$Type)
            df$Chr[nchar(df$Chr) > 50] <- ""

            # Ensembl gene ID hyperlinks
            ix <- grepl("ENS", df$"Ensembl Gene ID") # Change to grepl("^ENS", ...) to anchor to start of string if needed
            if (sum(ix) > 0) {
              tem <- paste0(
                "<a href='http://www.ensembl.org/id/",
                df$"Ensembl Gene ID",
                "' target='_blank'>",
                df$"Ensembl Gene ID",
                "</a>"
              )
              df$"Ensembl Gene ID"[ix] <- tem[ix]
            }

            # Entrez gene ID hyperlinks
            ix <- !is.na(as.numeric(df$Entrez))
            if (sum(ix) > 0) {
              tem <- paste0(
                "<a href='https://www.ncbi.nlm.nih.gov/gene/",
                df$Entrez,
                "' target='_blank'>",
                df$Entrez,
                "</a>"
              )
              df$Entrez[ix] <- tem[ix]
            }
          }
          return(df)
        }) # isolate
      },
      digits = 4,
      spacing = "s",
      striped = TRUE,
      bordered = TRUE,
      width = "auto",
      hover = T,
      sanitize.text.function = function(x) x
    )

    # Download full gene info CSV ----------------------------------------------
    output$downloadGeneInfo <- downloadHandler(
      filename = function() "geneInfo.csv",
      content = function(file) {
        write.csv(conversionTableData(), file, row.names = FALSE)
      }
    )
  })
}
