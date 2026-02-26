###################################################
# Author: Steven Ge Xijin.Ge@sdstate.edu
# Lab: Ge Lab
# R version 4.1
# Project: ShinyGO
# File: ui.R
# Purpose of file: ui logic of app
# Start data: NA (mm-dd-yyyy)
# Data last modified: 04-8-2022
#######################################################
library(shiny, verbose = FALSE)
library(shinyBS, verbose = FALSE) # for popup figures
library(plotly) # interactive network plot
library(visNetwork)
library("reactable", verbose = FALSE)

ui <- fluidPage(
  # Reduce the space between label and widgets, globally
  tags$head(
    tags$style(HTML(
      "label { font-size:100%; font-family:Times New Roman; margin-bottom:-15px; }"
    ))
  ),
  shinybusy::add_busy_spinner(spin = "fading-circle"), # Add spinner
  sidebarLayout(
    # Sidebar module (ALL sidebar functionality)
    mod_01_sidebar_ui(id = "sidebar"), # sidebarPanel

    mainPanel(
      tabsetPanel(
        id = "tabs", type = "tabs",
        mod_02_enrichment_ui(id = "enrichment"), # Enrichment tab
        mod_03_chart_ui(id = "chart"), # Chart tab
        mod_04_tree_ui(id = "tree"), # Tree tab
        mod_05_network_ui(id = "network"), # Network tab
        mod_06_kegg_ui(id = "kegg"), # KEGG tab
        mod_07_genes_ui(id = "genes"), # Genes tab
        mod_13_groups_ui(id = "groups"), # Groups tab
        mod_08_plots_ui(id = "plots"), # Plots tab
        mod_09_genome_ui(id = "genome"), # Genome tab
        # mod_12_promoter_ui(id = "promoter"), # Promoter tab (hidden — uncomment to enable)
        mod_10_string_ui(id = "string"), # STRING tab
        mod_11_about_ui(id = "about"), # About tab

      ), # tabsetPanel

      # bsModals (pop-up figures)
      # STRING PPI modal
      mod_10_string_modal_ui(id = "string"),
      # Static Network Plot Modal
      mod_05_network_modal_ui(id = "network"),
      # Gene ID Examples Modal
      bsModal("geneIDexamples", "What the gene IDs in our database look like?", "MGeneIDexamples",
        size = "large",
        selectizeInput(
          inputId = "userSpeciesIDexample",
          label = "Select or search for species", choices = NULL
        ),
        tableOutput("showGeneIDs4Species")
      ),
      # Genome Static Plot Modal
      mod_09_genome_modal_ui(id = "genome")

    ) # mainPanel
  ) # sidebarLayout

  #, tags$head(includeScript("google_analytics.js")), # tracking usage
  #tags$head(includeHTML(("google_analytics_GA4.html")))
) # fluidPage
