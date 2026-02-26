#' 13_groups UI Function
#'
#' @description A shiny Module for the Groups tab, which shows query genes
#'   grouped by high-level GO biological process terms.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_13_groups_ui <- function(id) {
  ns <- shiny::NS(id)

  tabPanel("Groups",
    value = 7,
    downloadButton(ns("downloadGrouping"), "Download"),
    h5("Your genes are grouped by functional categories defined by high-level GO terms. "),
    tableOutput(ns("grouping"))
  )
}


#' 13_groups Server Functions
#'
#' @param id                Module ID string
#' @param sidebar_values    Reactive list returned by mod_01_sidebar_server
#' @param enrichment_values Reactive list returned by mod_02_enrichment_server
#'
#' @noRd
#'
#' @importFrom shiny moduleServer
mod_13_groups_server <- function(id, sidebar_values, enrichment_values) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$grouping <- renderTable(
      {
        if (sidebar_values$goButton() == 0) {
          return()
        }
        myMessage <- "Just a minute. Matching your genes with level 2 and level 3 Gene Ontology biological process terms.
	       This can take up to 1 minute as we have to glue together a large number of gene names. "
        withProgress(message = sample(quotes, 1), detail = myMessage, {
          tem <- enrichment_values$significantOverlaps()

          incProgress(1, detail = paste("Done"))
        })
        tem$groupings
      },
      digits = 1,
      spacing = "s",
      striped = TRUE,
      bordered = TRUE,
      width = "auto",
      hover = T
    )

    output$downloadGrouping <- downloadHandler(
      filename = function() {
        "GO_Groups.csv"
      },
      content = function(file) {
        write.csv(enrichment_values$significantOverlaps()$groupings, file, row.names = FALSE)
      }
    )

  })
}
