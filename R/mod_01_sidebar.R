#' 01_sidebar UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_01_sidebar_ui <- function(id) {
  ns <- shiny::NS(id)
  tagList(
  sidebarPanel(
    titlePanel("ShinyGO 0.85.1",
      tags$head(tags$link(rel = "icon", type = "image/png", href = "favicon.png"),
            tags$title("ShinyGO 0.85"))
    ),
    # use conditional panel to hide the selectOrg input
    conditionalPanel(
      condition = "0", # hide the selectOrg input, always
      selectInput(
        inputId = ns("selectOrg"),
        label = NULL,
        selectize = TRUE,
        choices = setNames(-9606, "Human"), # Human is selected by default (ID from convertIDs.db)
        selected = setNames(-9606, "Human")
      )
    ),
    fluidRow(
      column(
        width = 6,
        textOutput(ns("selected_species"))
      ),
      column(
        width = 6,
        align = "left",
        # hide the change species button, once submit button is clicked. 
        # this avoids errors where some analyses are not updated when user changes species in the middle of an analysis
        conditionalPanel(
          condition = sprintf("input['%s'] == 0", ns("goButton")),
        # Species list and genome assemblies ----------
          actionButton(
            inputId = ns("genome_assembl_button"),
            label = strong("Change species")
          ),
          tippy::tippy_this(
            ns("genome_assembl_button"),
            "Search and click a row to select a species. ",
            theme = "light-border"
          )
        )
      )
    ),
    tags$head(
      tags$style(
        paste0("#", ns("selected_species"), " {color: red;font-size: 15px;font-style: italic;}")
      )
    ),
    br(),
    fluidRow(
      column(8, 
        conditionalPanel(
          condition = sprintf("input['%s'] == 0", ns("goButton")),
          actionButton(ns("useDemo1"), "Demo genes")
        )
      ),
      # column(4,   actionButton("useDemo2", "Demo 2"),	  	  ),
      column(4, p(HTML("<div align=\"right\"> <A HREF=\"javascript:history.go(0)\">Reset</A></div>")))
    ),
    tags$style(type = "text/css", "textarea {width:100%}"),
    tags$textarea(
      id = ns("input_text"), placeholder = "Change the species if it is not human. Then just paste a list of genes and click Submit. Gene IDs can be NCBI, Ensembl, symbol, or other common types.",
      rows = 8, ""
    ),
    conditionalPanel(
      condition = sprintf("input['%s'] != 0", ns("goButton")),
      textOutput(ns("mapping_stats"))
    ),
    fluidRow(
      column(8, actionButton(ns("backgroundGenes"), "Background (recommended)")),
      column(4, actionButton(ns("goButton"), strong("Submit")))
    ),
    br(),
    htmlOutput(ns("selectGO1")),
    fluidRow(
      column(
        6,
        numericInput(
          inputId = ns("minFDR"),
          label = h5("FDR cutoff"),
          value = 0.05, step = 0.01
        ),
        tippy::tippy_this(
          ns("minFDR"),
          "Minimum  P-value, ajusted using the FDR (false discovery rate)
        method. P-value is derived from hypergeometric distribution.
          Really significant FDR are between 1E-5 to 1E-20. Be cautious
          when you get an FDR of 1E-2 or 1E-3, as thousands of
          gene sets are tested.",
          theme = "light-border"
        )
      ),
      column(
        6,
        selectInput(ns("maxTerms"), h5("# pathways to show"),
          choices = list(
            "10" = 10,
            "15" = 15,
            "20" = 20,
            "25" = 25,
            "30" = 30,
            "40" = 40,
            "50" = 50,
            "60" = 60,
            "80" = 80,
            "100" = 100,
            "200" = 200,
            "500" = 500
          ),
          selected = "20",
          selectize = FALSE
        ),
        tippy::tippy_this(
          ns("maxTerms"),
          "How many top pathways to show.
          You can download nearly all significant ones.
          We typically recommend focusing on the top 10 to 20 pathways.
          If you go down the list,
          you can always find the one that help you tell the story you
          want to tell.",
          theme = "light-border"
        )
      )
    ),
    # tags$style(type='text/css', "#minFDR { width:100%;   margin-top:-15px}"),
    # selectInput("selectOrg", label = NULL,"Best matching species",width='100%'),


    fluidRow(
      column(
        width = 6,
        numericInput(ns("minSetSize"),
          label = h5("Pathway size: Min."),
          min   = 2,
          max   = 30,
          value = 2,
          step  = 1
        ),
        tippy::tippy_this(
          ns("minSetSize"),
          "Smaller pathways can introduce noise. Generally safe to incrase to 10 or 15.
          It is automatically raised to 10 when \"Sort by Fold Enrichment \" is selected.",
          theme = "light-border"
        )
      ),
      column(
        width = 6,
        numericInput(ns("maxSetSize"),
          label = h5("Max."),
          min   = 1000,
          max   = 20000,
          value = 5000,
          step  = 200
        ),
        tippy::tippy_this(
          ns("maxSetSize"),
          "Big gene sets, such as those associated with top-level GO term
          \"Cellular Process\", are less informative, but
          tend to have small P values due to increased power.",
          theme = "light-border"
        )
      )
    ), # fluidRow
    # tags$style(type='text/css', "#minSetSize { width:100%;   margin-top:-12px}"),
    # tags$style(type='text/css', "#maxSetSize { width:100%;   margin-top:-12px}"),
    fluidRow(
      column(
        width = 6,
        checkboxInput(
          ns("removeRedundantSets"),
          "Remove redundancy",
          value = TRUE
        ),
        tippy::tippy_this(
          ns("removeRedundantSets"),
          "Similar pathways sharing 95% of genes are represented by the most significant pathway.",
          theme = "light-border"
        )
      ),
      column(
        width = 6,
        checkboxInput(
          ns("abbreviatePathway"),
          "Abbreviate pathways",
          value = TRUE
        ),
        tippy::tippy_this(
          ns("abbreviatePathway"),
          "Positive regulation --> Pos. reg.",
          theme = "light-border"
        )
      )
    ), # fluidRow

    fluidRow(
      column(
        width = 6,
        checkboxInput(
          ns("gene_count_pathwaydb"),
          "Use pathway DB for gene counts",
          value = FALSE
        ),
        tippy::tippy_this(
          ns("gene_count_pathwaydb"),
          "If turned on, a gene must match at least one pathway in the selected pathway database.
          Otherwise, this gene is ignored when calculating enrichment. Be cautious
          when the selected pathway database is small, such as KEGG. ",
          theme = "light-border"
        )
      ),
      column(
        width = 6,
        checkboxInput(
          inputId = ns("show_pathway_id"),
          label = "Show pathway IDs",
          value = FALSE
        ),
        tippy::tippy_this(
          ns("show_pathway_id"),
          "If selected, pathway IDs, such as Path:mmu04115 and GO:0042770,  will be appended to pathway name.",
          theme = "light-border"
        )
      )
    ), # fluidRow

    actionButton(ns("MGeneIDexamples"), "Gene IDs examples"),
    tippy::tippy_this(
      ns("MGeneIDexamples"),
      "Show some example gene IDs in our database for a specific species.",
      theme = "light-border"
    ),
    h5("Try ", a(" iDEP", href = "https://bioinformatics.sdstate.edu/idep/", target = "_blank"), "for RNA-Seq data analysis")#,
    #tableOutput("species")
  ), #sidebarPanel

  # Background genes modal
  bsModal("BackgroundGenes", "Customized background genes (recommended)", "backgroundGenes",
    size = "large",
    tags$textarea(
      id = ns("input_text_b"),
      placeholder = "
Paste all genes from which the gene list is derived. These are all
genes whose expression or other activity that you measured.
This could be all the genes on a DNA microarray or all the genes
detected by a proteomics experiment.

By default, we compare your gene list with a background of all
protein-coding genes in the genome. When your genes are not selected
from genome-wide data, customized background genes might yield more
accurate results for enrichment analysis. For gene lists derived from
a typical RNA-seq dataset, many  use the subset of genes with detectable
expression, typically the genes passed a minimum filter.
We can also customize background genes to overcome bias in selection.
Currently only less than 30,000 genes are accepted.",
      rows = 20,
      ""
    )
  ) # bsModal 3
  ) # tagList
}


#' 01_sidebar Server Functions
#'
#' @noRd
#'
#' @importFrom shiny moduleServer
mod_01_sidebar_server <- function(id, speciesChoice, exampleGeneList, converted, orgInfo) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # connect to species specific database
    observeEvent(input$selectOrg, {

      # Query idIndex from the main convert database (convertIDs.db)
      # not from species-specific database, as idIndex is a global table
      idIndex <- dbGetQuery(
        convert,
        "select * from idIndex;"
      )
      # Note: 'convert' is the global database connection from global.R
      # We don't need to connect to convert_species just for idIndex
    })

    # Pop-up modal for gene assembl information ----
    observeEvent(input$genome_assembl_button, {
      shiny::showModal(
        shiny::modalDialog(
          size = "l",
          title = "Click on a row to select a species",
          p("Search annotated species by common or scientific names,
            or NCBI taxonomy id. Click on a row to select.
            Use annotation in STRING-db as a last resort.
            "),
          easyClose = TRUE,

          DT::renderDataTable({
            # column names differ between production (academicName, taxon_id)
            # and mock/local databases (name, name2, taxon)
            academic_col <- if ("academicName" %in% colnames(orgInfo)) "academicName" else "name"
            name_col     <- if ("academicName" %in% colnames(orgInfo)) "name" else "name2"
            taxon_col    <- if ("taxon_id"     %in% colnames(orgInfo)) "taxon_id"     else "taxon"
            df <- orgInfo[, c("ensembl_dataset", academic_col, name_col, taxon_col, "group")]
            colnames(df) <- c(
              "Ensembl/STRING-db ID",
              "Academic Name",
              "Name (Assembly)",
              "Taxonomy ID",
              "Source"
            )
            row.names(df) <- NULL
            DT::datatable(
              df,
              selection = "single",
              options = list(
                lengthChange = FALSE,
                pageLength = 20,
                scrollY = "400px"
              ),
              callback = DT::JS(
                paste0(
                  "table.on('click', 'tr', function() {
                      var data = table.row(this).data();
                      if (data) {
                        Shiny.setInputValue('", ns("clicked_row"), "', data[0]);
                      }
                    });"
                )
              ),
              rownames = FALSE
            )
          })
        )
      )
    })

    observeEvent(input$clicked_row, {
      # find species ID from ensembl_dataset
      selected <- find_species_id_by_ensembl(
        input$clicked_row,
        orgInfo
      )
      # assign name
      selected <- setNames(
        selected,
        find_species_by_id_name(selected, orgInfo)
      )

      updateSelectizeInput(
        session = session,
        inputId = "selectOrg",
        choices = selected,
        selected = selected,
        server = TRUE
      )
      # update species name
      selected_species_name(find_species_by_id_name(selected, orgInfo))
    })

    # default species name
    selected_species_name <- reactiveVal("Human")

    output$selected_species <- renderText({
      tem <- input$clicked_row
      selected_species_name()
    })

    # dynamic UI for pathway database selection
    output$selectGO1 <- renderUI({ # gene set for pathway analysis
      if (input$goButton == 0) {
        return(NULL)
      }

      choices <- gmtCategory(converted(), input$selectOrg)
      if (length(choices) > 12) { # more than 12 categories in human and mouse, we default to GOBP
        selected <- "GOBP"
      } else { # otherwise all gene sets
        selected <- "All"
      }
      if ("KEGG" %in% choices) {
        selected <- "KEGG"
      }


      selectInput(ns("selectGO"),
        label = h5("Pathway database:"),
        choices = choices,
        selected = selected
      )
    })

    # Gene ID examples modal trigger
    observe({
      # for gene ID example
      updateSelectizeInput(session, "userSpeciesIDexample", choices = speciesChoice, selected = speciesChoice[1])

      # load demo data when clicked
      if (isTRUE(input$useDemo1 > 0)) {
        updateTextInput(session, "input_text", value = ExampleGeneList1)
      }

      # update species for STRING-db related API access
      # tried to solve the double reflashing problems
      # https://stackoverflow.com/questions/30991900/avoid-double-refresh-of-plot-in-shiny
    })

    # Mapping stats
    output$mapping_stats <- renderText({
      req(input$goButton)
      req(converted())
      n_genes <- length(converted()$originalIDs)
      n_mapped <- length(converted()$IDs)

      paste0(n_genes, " IDs mapped to ", n_mapped, " (", round(n_mapped / n_genes * 100, 0), "%) ", tolower(converted()$species$name2), " genes.")
    })

    # Background genes modal trigger
    observeEvent(input$backgroundGenes, {
      toggleModal(session, "BackgroundGenes", toggle = "open")
    })

    # Gene ID examples modal trigger
    observeEvent(input$MGeneIDexamples, {
      toggleModal(session, "geneIDexamples", toggle = "open")
    })

    list(
      # Species
      selectOrg = reactive({ val <- input$selectOrg; if (is.null(val)) -9606L else val }),

      # Gene input
      input_text = reactive({ val <- input$input_text; if (is.null(val)) "" else val }),
      input_text_b = reactive({ val <- input$input_text_b; if (is.null(val)) "" else val }),
      goButton = reactive({ val <- input$goButton; if (is.null(val)) 0L else val }),
      mapping_stats = reactive(output$mapping_stats),

      # Pathway database
      selectGO = reactive(input$selectGO),  # from dynamic UI

      # Parameters
      minFDR = reactive(input$minFDR),
      maxTerms = reactive(input$maxTerms),
      minSetSize = reactive(input$minSetSize),
      maxSetSize = reactive(input$maxSetSize),

      # Display options
      removeRedundantSets = reactive(input$removeRedundantSets),
      abbreviatePathway = reactive(input$abbreviatePathway),
      gene_count_pathwaydb = reactive(input$gene_count_pathwaydb),
      show_pathway_id = reactive(input$show_pathway_id)
    )

  })
}