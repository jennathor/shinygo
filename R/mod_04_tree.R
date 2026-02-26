#' 04_tree UI Function
#'
#' @description A shiny Module for the hierarchical clustering tree tab.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_04_tree_ui <- function(id) {
  ns <- shiny::NS(id)

  tabPanel("Tree",
    value = 4,
    h5("A hierarchical clustering tree summarizes the correlation among significant pathways
                    listed in the Enrichment tab. Pathways with many shared genes are clustered together.
                      Bigger dots indicate more significant P-values. The width of the plot can be
                      changed by adjusting the width of your browser window."),
    fluidRow(
      column(width = 3, selectInput(
        inputId = ns("treeChartAspectRatio"),
        label = h5("Aspect Ratio"),
        choices = .1 * (5:40),
        selected = 2
      )),
      column(3, style = "margin-top: 25px;", mod_download_images_ui(ns("download_tree"), label = "Download"))
    ),
    plotOutput(ns("GOTermsTree"))
  )
}


#' 04_tree Server Functions
#'
#' @noRd
#'
#' @importFrom shiny moduleServer
mod_04_tree_server <- function(id, sidebar_values, enrichment_values) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    significantOverlaps2 <- reactive({
      if (sidebar_values$goButton() == 0) {
        return()
      }
      tem <- sidebar_values$input_text_b() # just to make it re-calculate if user changes background

      tem <- enrichment_values$significantOverlaps()
      if (dim(tem$x)[2] == 1) {
        return(NULL)
      }
      tem <- tem$x
      colnames(tem) <- c("adj.Pval", "nGenesList", "nGenesCategor", "Fold", "Pathways", "URL", "Genes")
      tem$Pathways <- gsub(".*'_blank'>|</a>", "", tem$Pathways) # remove URL
      tem$Direction <- "Diff"
      # remove pathway ID  only in Ensembl species
      if (!sidebar_values$show_pathway_id() && as.integer(sidebar_values$selectOrg()) > 0) {
        tem$Pathways <- remove_pathway_id(tem$Pathways, sidebar_values$selectGO())
      }
      tem
    })

    tree_plot <- reactive({
      if (sidebar_values$goButton() == 0) {
        return(NULL)
      }
      if (is.null(significantOverlaps2())) {
        return(NULL)
      }
      tem <- sidebar_values$maxTerms()
      p <- enrichmentPlot(significantOverlaps2(), 45)
      return(p)
    })

    output$GOTermsTree <- renderPlot(
      {
        if (sidebar_values$goButton() == 0) {
          return(NULL)
        }
        if (is.null(significantOverlaps2())) {
          return(NULL)
        }
        tem <- sidebar_values$maxTerms()
        # enrichmentPlot(significantOverlaps2(), 56  )
        tree_plot()
      },
      height = function() {
        round(max(350, min(2500, round(18 * as.numeric(sidebar_values$maxTerms())))))
      },
      width = function() {
        width1 <- round(max(350, min(1000, round(18 * as.numeric(sidebar_values$maxTerms())))) * as.numeric(input$treeChartAspectRatio))
        return(min(width1, 1000)) # max width is 1000
      }
    )

    mod_download_images_server(
      "download_tree",
      filename = "tree_plot",
      figure = reactive({
        tree_plot()
      }),
      width = 10,
      height = round(10 / as.numeric(input$treeChartAspectRatio), 1)
    )

  })
}
