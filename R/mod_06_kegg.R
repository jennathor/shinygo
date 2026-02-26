#' 06_kegg UI Function
#'
#' @description A shiny Module for the KEGG pathway diagram tab.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_06_kegg_ui <- function(id) {
  ns <- shiny::NS(id)
  tabPanel("KEGG",
    value = 2,
    conditionalPanel(
      "input['sidebar-selectGO'] != 'KEGG'",
      br(), br(),
      h5("Please select KEGG from the pathway databases to conduct enrichment analysis first.
          Then you can visualize your genes on any of the significant pathways. Only for some species.")
    ),
    conditionalPanel(
      "input['sidebar-selectGO'] == 'KEGG'",
      br(),
      uiOutput(ns("listSigPathways")),
      br(), imageOutput(ns("KeggImage"), width = "100%", height = "100%"),
      h5("Your genes are highlighted in red. Downloading pathway diagram from KEGG can take 3 minutes. ")
    )
  )
}

#' 06_kegg Server Functions
#'
#' @param id          Module ID string
#' @param sidebar_values  Reactive list returned by mod_01_sidebar_server
#' @param enrichment_values  Reactive list returned by mod_02_enrichment_server
#' @param converted   Reactive returning the converted gene ID object
#'                    (list with $IDs and $species)
#' @param keggSpeciesID  Data frame mapping ensembl_dataset -> name -> KEGG code
#'                       (created in global.R from orgInfo)
#'
#' @noRd
#'
#' @importFrom shiny moduleServer
mod_06_kegg_server <- function(id, sidebar_values, enrichment_values, converted, keggSpeciesID) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Dynamic dropdown of significant KEGG pathways ----------------------------
    output$listSigPathways <- renderUI({
      tem <- sidebar_values$selectOrg()
      if (sidebar_values$goButton() == 0 | is.null(enrichment_values$significantOverlaps())) {
        return(NULL)
      }

      tem <- enrichment_values$significantOverlaps()

      if (dim(tem$x)[2] == 1) {
        return(NULL)
      }
      tem$x <- tem$x[tem$x[, 3] < 1000, ] # remove pathways with more than 1000 genes. Very slow.
      tem$x <- tem$x[order(-tem$x[, 4]), ] # sort by fold-enrichment
      choices <- tem$x[, 5]
      selectInput(ns("sigPathways"),
        label = "Select a significant KEGG pathway:",
        choices = choices
      )
    })

    # KEGG pathway image -------------------------------------------------------
    output$KeggImage <- renderImage(
      {
        req(input$sigPathways) # also stops on "", NA — not just NULL

        # First generate a blank image. Otherwise return(NULL) gives us errors.
        outfile <- tempfile(fileext = ".png")
        png(outfile, width = 400, height = 300)
        frame()
        dev.off()
        blank <- list(
          src = outfile,
          contentType = "image/png",
          width = 400,
          height = 300,
          alt = " "
        )

        if (sidebar_values$goButton() == 0) return(blank)
        if (is.null(sidebar_values$selectGO())) return(blank)
        if (sidebar_values$selectGO() != "KEGG") return(blank)
        if (is.null(enrichment_values$significantOverlaps())) return(blank)

        if (!requireNamespace("pathview", quietly = TRUE)) {
          showNotification(
            "The 'pathview' package is required for KEGG pathway diagrams but is not installed.",
            type = "error", duration = NULL
          )
          return(blank)
        }
        library(pathview, verbose = FALSE)

        # modify my.keggview.native to run in pathview's namespace so it can
        # access internal helpers (colorpanel2, render.kegg.node, col.key, etc.)
        # http://stackoverflow.com/questions/23279904
        tmpfun <- get("keggview.native", envir = asNamespace("pathview"))
        environment(my.keggview.native) <- environment(tmpfun)
        attributes(my.keggview.native) <- attributes(tmpfun) # don't know if this is really needed

        isolate({
          withProgress(message = "Rendering KEGG pathway plot", {
            incProgress(1 / 5, "Loading the pathview package")

            Species <- converted()$species[1, 1]
            fold <- convertEnsembl2Entrez(converted()$IDs, Species)

            fold <- fold$entrezgene_id
            keggSpecies <- as.character(keggSpeciesID[which(keggSpeciesID[, 1] == Species), 3])

            if (nchar(keggSpecies) <= 2) {
              return(blank)
            } # not in KEGG

            # kegg pathway id
            incProgress(1 / 2, "Download pathway graph from KEGG.")

            # find pathway id: "Path:hsa04110 Cell cycle" --> "hsa04110"
            pathID <- gsub(" .*", "", input$sigPathways)
            pathID <- gsub("Path:", "", pathID)

            if (nchar(pathID) < 3) {
              return(blank)
            }
            randomString <- gsub(".*file", "", tempfile())
            tempFolder <- tempdir()
            outfile <- paste(tempFolder, "/", pathID, ".", randomString, ".png", sep = "")

            pv.out <- mypathview(
              gene.data = fold, pathway.id = pathID,
              kegg.dir = tempFolder, out.suffix = randomString,
              species = keggSpecies, kegg.native = TRUE
            )

            list(
              src = outfile,
              contentType = "image/png",
              width = "100%",
              height = "100%",
              alt = "KEGG pathway image."
            )
          })
        })
      },
      deleteFile = TRUE
    )
  })
}
