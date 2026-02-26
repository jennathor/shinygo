#' 12_promoter UI Function
#'
#' @description A shiny Module for the Promoter tab. Compares transcription
#'   factor binding motifs in gene promoters against genome-wide background.
#'   Tab is currently hidden; uncomment the call in ui.R to enable it.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_12_promoter_ui <- function(id) {
  ns <- shiny::NS(id)

  tabPanel("Promoter",
    value = 10,
    radioButtons(ns("radio"), label = NULL, choices = list(
      "Upstream 300bp as promoter" = 300,
      "Upstream 600bp as promoter" = 600
    ), selected = 300),
    tableOutput(ns("promoter")),
    downloadButton(ns("downloadPromoter"), "Download"),
    h5("The promoter sequences of your genes are compared with those of the
        other genes in the genome in terms of transcription factor (TF) binding motifs.
        \"*Query gene\" indicates a transcription factor coded by a gene included in
        your list.")
  )
}


#' 12_promoter Server Functions
#'
#' @param id          Module ID string
#' @param sidebar_values  Reactive list returned by mod_01_sidebar_server
#' @param converted   Reactive returning the gene ID conversion result
#'
#' @noRd
#'
#' @importFrom shiny moduleServer
mod_12_promoter_server <- function(id, sidebar_values, converted) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    promoterData <- reactive({
      if (sidebar_values$goButton() == 0) {
        return()
      }
      tem <- input$radio
      tem <- sidebar_values$selectOrg()
      isolate({
        myMessage <- "Promoter analysis"
        withProgress(message = sample(quotes, 1), detail = myMessage, {
          tem <- promoter(converted(), sidebar_values$selectOrg(), input$radio)
          incProgress(1, detail = paste("Done"))
        })

        if (is.null(tem)) {
          return(as.data.frame("ID not recognized."))
        } else {
          return(tem)
        }
      }) # avoid showing things initially
    })

    output$promoter <- renderTable(
      {
        if (sidebar_values$goButton() == 0) {
          return()
        }
        tem <- input$radio
        tem <- sidebar_values$selectOrg()
        isolate({
          promoterData()
        }) # avoid showing things initially
      },
      digits = -1,
      spacing = "s",
      striped = TRUE,
      bordered = TRUE,
      width = "auto",
      hover = T
    )

    output$downloadPromoter <- downloadHandler(
      filename = function() {
        "promoterMotif.csv"
      },
      content = function(file) {
        write.csv(promoterData(), file, row.names = FALSE)
      }
    )

  })
}
