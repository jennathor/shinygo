#' 08_plots UI Function
#'
#' @description A shiny Module for the Plots tab, showing gene characteristic
#'   density plots and barplots comparing query genes vs the genome/background.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_08_plots_ui <- function(id) {
  ns <- shiny::NS(id)

  tabPanel("Plots",
    value = 8,
    h5("The characteristics of your genes are compared with the rest in the genome. Chi-squared and Student's
        t-tests are run to see if your genes have special characteristics when compared with all the other genes or, if uploaded, a customized background."),
    fluidRow(
      column(
        width = 4,
        mod_download_images_ui(ns("download_gene_plot_dist"), "Download density plots")
      ),
      column(
        width = 4,
        mod_download_images_ui(ns("download_gene_barplot"), "Download barplots")
      )
    ),
    br(),
    plotOutput(ns("genePlot2"), inline = TRUE, width = "auto", height = "auto"),
    plotOutput(ns("gene_barplot"), inline = TRUE, width = "auto", height = "auto")
  )
}


#' 08_plots Server Functions
#'
#' @param id                      Module ID string
#' @param sidebar_values          Reactive list returned by mod_01_sidebar_server
#' @param geneInfoLookup          Reactive returning gene info for query genes
#' @param geneInfoLookup_background  Reactive returning gene info for background genes
#' @param converted_background    Reactive returning background gene conversion result
#' @param ggplot2_theme           Reactive returning the current ggplot2 theme string
#'
#' @noRd
#'
#' @importFrom shiny moduleServer
mod_08_plots_server <- function(id, sidebar_values, geneInfoLookup,
                                geneInfoLookup_background, converted_background,
                                ggplot2_theme) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # barplots using R base graphics
    gene_barplot_object <- reactive({
      if (sidebar_values$goButton() == 0) {
        return()
      }
      tem <- sidebar_values$selectOrg()
      isolate({
        withProgress(message = "Ploting gene characteristics", {
          x <- geneInfoLookup()
          x2 <- x[which(x$gene_biotype == "protein_coding"), ] # only coding for some analyses


          # background genes	--------------------------------------------------------------
          xB <- geneInfoLookup_background()
          convertedB <- converted_background()
          if (!is.null(xB) &&
            !is.null(convertedB) ) {
            x <- x[x$Set == "List", ] # remove background from selected genes
            xB <- xB[xB$Set == "List", ] # remove Genome genes from background
            xB$Set <- "Background"
            x <- rbind(x, xB)
            x2 <- x[which(x$gene_biotype == "protein_coding"), ] # only coding for some analyses
          }
          # end background genes


          if (dim(x)[1] >= minGenes) # only making plots if more than 20 genes
            { # only plot when there 10 genes or more   # some columns have too many missing values
              par(mfrow = c(4, 1))
              par(mar = c(8, 6, 8, 2))
              # chromosomes
              if (sum(!is.na(x$chromosome_name)) >= minGenes && length(unique(x$chromosome_name)) > 2 && length(which(x$Set == "List")) > minGenes) {
                freq <- table(x$chromosome_name, x$Set)
                freq <- as.matrix(freq[which(nchar(row.names(freq)) < 50), ]) # remove unmapped chromosomes
                if (dim(freq)[2] > 1 && dim(freq)[1] > 1) { # some organisms do not have fully seuqence genome: chr. names: scaffold_99816
                  Pval <- chisq.test(freq)$p.value
                  sig <- paste("Distribution of query genes on chromosomes \nChi-squared test P=", formatC(Pval, digits = 2, format = "G"))

                  if (Pval < PvalGeneInfo2) {
                    sig <- paste(sig, "***")
                  } else
                  if (Pval < PvalGeneInfo1) {
                    sig <- paste(sig, "**")
                  } else
                  if (Pval < PvalGeneInfo) sig <- paste(sig, "*")

                  freq <- freq[order(as.numeric(row.names(freq))), ]
                  freq[, 1] <- freq[, 1] * colSums(freq)[2] / colSums(freq)[1] # expected
                  freq <- freq[, c(2, 1)] # reverse order
                  barplot(t(freq),
                    beside = TRUE, las = 3, col = c("red", "lightgrey"), ylab = "Number of Genes", main = sig,
                    cex.lab = 1.5, cex.axis = 2, cex.names = 2, cex.main = 1.5
                  )

                  legend("topright", c("List", "Expected"), pch = 15, col = c("red", "lightgrey"), bty = "n", cex = 2)
                }
              } else { # Create empty plot
                plot(x = 0:1, y = 0:1, ann = F, bty = "n", type = "n", xaxt = "n", yaxt = "n")
                text(
                  x = 0.5, y = 0.5, # Add text to empty plot
                  "Chromosome plot not available.",
                  cex = 1.8
                )
              }
              incProgress(1 / 8)


              # gene type
              if (sum(!is.na(x$gene_biotype)) >= minGenes && length(unique(x$gene_biotype)) > 2 && length(which(x$Set == "List")) > minGenes) {
                freq <- table(x$gene_biotype, x$Set)
                freq <- as.matrix(freq[which(freq[, 1] / colSums(freq)[1] > .01), ])
                if (dim(freq)[2] > 1 && dim(freq)[1] > 1) {
                  Pval <- chisq.test(freq)$p.value
                  sig <- paste("Distribution by gene type \nChi-squared test P=", formatC(Pval, digits = 2, format = "G"))
                  if (Pval < PvalGeneInfo2) {
                    sig <- paste(sig, "***")
                  } else
                  if (Pval < PvalGeneInfo1) {
                    sig <- paste(sig, "**")
                  } else
                  if (Pval < PvalGeneInfo) sig <- paste(sig, "*")
                  freq <- freq[order(freq[, 1], decreasing = T), ]
                  freq[, 1] <- freq[, 1] * colSums(freq)[2] / colSums(freq)[1]
                  tem <- gsub("protein_coding", "Coding", rownames(freq))
                  tem <- gsub("pseudogene", "pseudo", tem)
                  tem <- gsub("processed", "proc", tem)
                  row.names(freq) <- tem
                  par(mar = c(20, 6, 4.1, 2.1))
                  freq <- freq[, c(2, 1)] # reverse order

                  barplot(t(freq),
                    beside = TRUE, las = 2, col = c("red", "lightgrey"), ylab = "Number of Genes",
                    main = sig, cex.lab = 1.2, cex.axis = 1.2, cex.names = 1.2, cex.main = 1.2
                  )
                  legend("topright", c("List", "Expected"), pch = 15, col = c("red", "lightgrey"), bty = "n", cex = 2)
                }
              } else { # Create empty plot
                plot(x = 0:1, y = 0:1, ann = F, bty = "n", type = "n", xaxt = "n", yaxt = "n")
                text(
                  x = 0.5, y = 0.5, # Add text to empty plot
                  "Gene type plot not available.",
                  cex = 1.8
                )
              }


              incProgress(1 / 8)
              par(mar = c(12, 6, 4.1, 2.1))
              # N. exons

              if (sum(!is.na(x2$nExons)) >= minGenes && length(unique(x2$nExons)) > 2 && length(which(x2$Set == "List")) > minGenes) {
                freq <- table(x2$nExons, x2$Set)
                freq <- as.matrix(freq[which(freq[, 1] / colSums(freq)[1] > .02), ])
                if (dim(freq)[2] > 1 && dim(freq)[1] > 1) {
                  Pval <- chisq.test(freq)$p.value
                  sig <- paste("Number of exons (coding genes only) \nChi-squared test P=", formatC(Pval, digits = 2, format = "G"))
                  if (Pval < PvalGeneInfo2) {
                    sig <- paste(sig, "***")
                  } else
                  if (Pval < PvalGeneInfo1) {
                    sig <- paste(sig, "**")
                  } else
                  if (Pval < PvalGeneInfo) sig <- paste(sig, "*")
                  # freq <- freq[order(    freq[,1], decreasing=T), ]
                  freq[, 1] <- freq[, 1] * colSums(freq)[2] / colSums(freq)[1]
                  freq <- freq[, c(2, 1)] # reverse order
                  barplot(t(freq),
                    beside = TRUE, las = 2, col = c("red", "lightgrey"), ylab = "Number of Genes",
                    main = sig, xlab = c("Number of exons"), cex.lab = 1.5, cex.axis = 2, cex.names = 1.5, cex.main = 1.5
                  )
                  legend("topright", c("List", "Expected"), pch = 15, col = c("red", "lightgrey"), bty = "n", cex = 2)
                }
              } else { # Create empty plot
                plot(x = 0:1, y = 0:1, ann = F, bty = "n", type = "n", xaxt = "n", yaxt = "n")
                text(
                  x = 0.5, y = 0.5, # Add text to empty plot
                  "Exon plot not available.",
                  cex = 1.8
                )
              }
              incProgress(1 / 8)

              # Transcript count
              if (sum(!is.na(x2$transcript_count)) >= minGenes && length(unique(x2$transcript_count)) > 2 && length(which(x2$Set == "List")) > minGenes) {
                freq <- table(x2$transcript_count, x2$Set)
                freq <- as.matrix(freq[which(freq[, 1] / colSums(freq)[1] > .02), ])
                if (dim(freq)[2] > 1 && dim(freq)[1] > 1) {
                  Pval <- chisq.test(freq)$p.value
                  sig <- paste("Number of transcript isoforms per coding gene \nChi-squared test P=", formatC(Pval, digits = 2, format = "G"))
                  if (Pval < PvalGeneInfo2) {
                    sig <- paste(sig, "***")
                  } else
                  if (Pval < PvalGeneInfo1) {
                    sig <- paste(sig, "**")
                  } else
                  if (Pval < PvalGeneInfo) sig <- paste(sig, "*")
                  freq <- freq[order(freq[, 1], decreasing = T), ]
                  freq[, 1] <- freq[, 1] * colSums(freq)[2] / colSums(freq)[1]
                  freq <- freq[, c(2, 1)] # reverse order
                  barplot(t(freq),
                    beside = TRUE, las = 2, col = c("red", "lightgrey"), ylab = "Number of Genes",
                    main = sig, xlab = c("Number of transcripts per gene"), cex.lab = 1.5, cex.axis = 2, cex.names = 1.5, cex.main = 1.5
                  )
                  legend("topright", c("List", "Expected"), pch = 15, col = c("red", "lightgrey"), bty = "n", cex = 2)
                }
              } else { # Create empty plot
                plot(x = 0:1, y = 0:1, ann = F, bty = "n", type = "n", xaxt = "n", yaxt = "n")
                text(
                  x = 0.5, y = 0.5, # Add text to empty plot
                  "Transcript plot not available.",
                  cex = 1.8
                )
              }
              incProgress(1 / 8)
            } # if minGenes
          incProgress(1 / 8, detail = paste("Done"))
          return(recordPlot())
        })
      }) # isolate
    })

    output$gene_barplot <- renderPlot(
      {
        gene_barplot_object()
      },
      width = 600,
      height = 1500
    )

    mod_download_images_server(
      "download_gene_barplot",
      filename = "gene_characteristics_barplot",
      figure = reactive({
        gene_barplot_object()
      }),
      width = 8,
      height = 20
    )

    # density plots using ggplot2
    gene_density_plot <- reactive({
      if (sidebar_values$goButton() == 0) {
        return()
      }
      req(ggplot2_theme())
      req(sidebar_values$selectOrg())
      isolate({
        withProgress(message = "Ploting gene characteristics", {
          x <- geneInfoLookup()
          x2 <- x[which(x$gene_biotype == "protein_coding"), ] # only coding for some analyses

          # background genes	--------------------------------------------------------------
          xB <- geneInfoLookup_background()
          convertedB <- converted_background()
          if (!is.null(xB) &&
            !is.null(convertedB)) { # if more than 30k genes, ignore background genes.

            x <- x[x$Set == "List", ] # remove background from selected genes
            xB <- xB[xB$Set == "List", ] # remove Genome genes from background
            xB$Set <- "Background"
            x <- rbind(x, xB)
            x2 <- x[which(x$gene_biotype == "protein_coding"), ] # only coding for some analyses
          }
          # end background genes


          if (dim(x)[1] >= minGenes) # only making plots if more than 20 genes
            { # only plot when there 10 genes or more   # some columns have too many missing values
              # par(mfrow=c(10,1))
              # par(mar=c(8,6,8,2))

              # increase fonts
              theme_set(theme_gray(base_size = 20))

              # Coding Sequence length
              if (sum(!is.na(x2$cds_length)) >= minGenes && length(unique(x2$cds_length)) > 2 &&
                length(which(x2$Set == "List")) > minGenes) {
                Pval <- t.test(log(cds_length) ~ Set, data = x2)$p.value
                sig <- mark_significance(Pval, PvalGeneInfo2, PvalGeneInfo1, PvalGeneInfo)


                p1 <- ggplot(x2, aes(cds_length, fill = Set, colour = Set)) +
                  geom_density(alpha = 0.1) +
                  scale_x_log10() +
                  labs(x = "Coding sequence length (bp)", y = "Density") +
                  annotate("text", x = min(x2$cds_length) + 50, y = .5, label = sig, size = 6) +
                  # annotate("text",x= max(x2$cds_length), y = densMode(x2$cds_length)$y, label=sig, size=8, hjust=1) +
                  guides(color = guide_legend(nrow = 2)) +
                  theme(
                    legend.key = element_rect(color = NA, fill = NA),
                    legend.key.size = unit(1.2, "line")
                  ) +
                  theme(plot.margin = unit(c(0, 0, 1, 0), "cm"))
              } else {
                p1 <- fake_plot("Coding Sequence length plot not available.")
              }

              incProgress(1 / 8)

              # Transcript length------------
              if (sum(!is.na(x2$transcript_length)) >= minGenes &&
                length(unique(x2$transcript_length)) > 2 &&
                length(which(x2$Set == "List")) > minGenes) {
                Pval <- t.test(log(transcript_length) ~ Set, data = x2[which(!is.na(x2$transcript_length)), ])$p.value
                sig <- mark_significance(Pval, PvalGeneInfo2, PvalGeneInfo1, PvalGeneInfo)

                p2 <- ggplot(x2, aes(transcript_length, fill = Set, colour = Set)) +
                  geom_density(alpha = 0.1) +
                  scale_x_log10() +
                  annotate("text", x = min(x2$transcript_length) + 100, y = .5, label = sig, size = 6) +
                  # annotate("text",x= max(x2$transcript_length), y = densMode(x2$transcript_length)$y, label=sig, size=8, hjust=1) +
                  labs(x = "Transcript length (bp)", y = "Density") +
                  guides(color = guide_legend(nrow = 2)) +
                  theme(
                    legend.key = element_rect(color = NA, fill = NA),
                    legend.key.size = unit(1.2, "line")
                  ) +
                  theme(plot.margin = unit(c(0, 0, 1, 0), "cm"))
              } else {
                p2 <- fake_plot("Transcript length plot not available.")
              }
              incProgress(2 / 8)

              # Genome span ------------

              if (sum(!is.na(x2$genomeSpan)) >= minGenes && length(unique(x2$genomeSpan)) > 2 && length(which(x2$Set == "List")) > minGenes) {
                Pval <- t.test(log(genomeSpan) ~ Set, data = x2[which(!is.na(x2$genomeSpan)), ])$p.value
                sig <- mark_significance(Pval, PvalGeneInfo2, PvalGeneInfo1, PvalGeneInfo)
                p3 <- ggplot(x2, aes(genomeSpan, fill = Set, colour = Set)) +
                  geom_density(alpha = 0.1) +
                  scale_x_log10() +
                  annotate("text", x = min(x2$genomeSpan) + 200, y = .5, label = sig, size = 6) +
                  # annotate("text",x= max(x2$genomeSpan), y = densMode(x2$genomeSpan)$y, label=sig, size=8, hjust=1) +
                  labs(x = "Genome span (bp)", y = "Density") +
                  guides(color = guide_legend(nrow = 2)) +
                  theme(
                    legend.key = element_rect(color = NA, fill = NA),
                    legend.key.size = unit(1.2, "line")
                  ) +
                  theme(plot.margin = unit(c(0, 0, 1, 0), "cm"))
              } else {
                p3 <- fake_plot("Genome span plot not available.")
              }

              incProgress(3 / 8)

              # 5' UTR ------------

              if (sum(!is.na(x2$FiveUTR)) >= minGenes && length(unique(x2$FiveUTR)) > 2 && length(which(x2$Set == "List")) > minGenes) {
                Pval <- t.test(log(FiveUTR) ~ Set, data = x2[which(!is.na(x2$FiveUTR) & x2$FiveUTR > 0), ])$p.value
                sig <- mark_significance(Pval, PvalGeneInfo2, PvalGeneInfo1, PvalGeneInfo)

                p4 <- ggplot(x2, aes(FiveUTR, fill = Set, colour = Set)) +
                  geom_density(alpha = 0.1) +
                  scale_x_log10() +
                  annotate("text",
                    x = min(x2[which(!is.na(x2$FiveUTR) & x2$FiveUTR > 0), "FiveUTR"]) + 5,
                    y = .5, label = sig, size = 6
                  ) +
                  # annotate("text",x= max(x2$FiveUTR), y = densMode(x2$FiveUTR)$y, label=sig, size=8, hjust=1) +
                  labs(x = "5' UTR length (bp)", y = "Density") +
                  guides(color = guide_legend(nrow = 2)) +
                  theme(
                    legend.key = element_rect(color = NA, fill = NA),
                    legend.key.size = unit(1.2, "line")
                  ) +
                  theme(plot.margin = unit(c(0, 0, 1, 0), "cm"))
              } else {
                p4 <- fake_plot("5' UTR plot not available.")
              }

              incProgress(4 / 8)

              # 3' UTR ------------
              if (sum(!is.na(x2$ThreeUTR)) >= minGenes && length(unique(x2$ThreeUTR)) > 2 && length(which(x2$Set == "List")) > minGenes) {
                Pval <- t.test(log(ThreeUTR) ~ Set, data = x2[which(!is.na(x2$ThreeUTR) & x2$ThreeUTR > 0), ])$p.value
                sig <- mark_significance(Pval, PvalGeneInfo2, PvalGeneInfo1, PvalGeneInfo)

                p5 <- ggplot(x2, aes(ThreeUTR, fill = Set, colour = Set)) +
                  geom_density(alpha = 0.1) +
                  scale_x_log10() +
                  annotate("text", x = min(x2[which(!is.na(x2$ThreeUTR) & x2$ThreeUTR > 0), "ThreeUTR"]) + 5, y = .5, label = sig, size = 6) +
                  # annotate("text",x= max(x2$ThreeUTR), y = densMode(x2$ThreeUTR)$y, label=sig, size=8, hjust=1) +
                  labs(x = "3' UTR length (bp)", y = "Density") +
                  guides(color = guide_legend(nrow = 2)) +
                  theme(
                    legend.key = element_rect(color = NA, fill = NA),
                    legend.key.size = unit(1.2, "line")
                  ) +
                  theme(plot.margin = unit(c(0, 0, 1, 0), "cm"))
              } else {
                p5 <- fake_plot("3' UTR plot not available.")
              }
              incProgress(5 / 8)

              # GC content ------------
              if (sum(!is.na(x2$percentage_gc_content)) >= minGenes &&
                length(unique(x2$percentage_gc_content)) > 2 &&
                length(which(x2$Set == "List")) > minGenes) {
                Pval <- t.test(percentage_gc_content ~ Set,
                  data = x2[which(!is.na(x2$percentage_gc_content) & x2$percentage_gc_content > 0), ]
                )$p.value
                sig <- mark_significance(Pval, PvalGeneInfo2, PvalGeneInfo1, PvalGeneInfo)

                p6 <- ggplot(x2, aes(percentage_gc_content, fill = Set, colour = Set)) +
                  geom_density(alpha = 0.1) +
                  # annotate("text",x= min(x2$percentage_gc_content)+5, y = .02, label=sig, size=8)+
                  annotate("text", x = max(x2$percentage_gc_content), y = densMode(x2$percentage_gc_content)$y, label = sig, size = 6, hjust = 1) +
                  labs(x = "GC content (%)", y = "Density") +
                  guides(color = guide_legend(nrow = 2)) +
                  theme(
                    legend.key = element_rect(color = NA, fill = NA),
                    legend.key.size = unit(1.2, "line")
                  ) +
                  theme(plot.margin = unit(c(0, 0, 1, 0), "cm"))
              } else {
                p6 <- fake_plot("GC content plot not available.")
              }

              incProgress(6 / 8)
              p1 <- refine_ggplot2(
                p = p1,
                gridline = FALSE,
                ggplot2_theme = ggplot2_theme()
              )
              p2 <- refine_ggplot2(
                p = p2,
                gridline = FALSE,
                ggplot2_theme = ggplot2_theme()
              )
              p3 <- refine_ggplot2(
                p = p3,
                gridline = FALSE,
                ggplot2_theme = ggplot2_theme()
              )
              p4 <- refine_ggplot2(
                p = p4,
                gridline = FALSE,
                ggplot2_theme = ggplot2_theme()
              )
              p5 <- refine_ggplot2(
                p = p5,
                gridline = FALSE,
                ggplot2_theme = ggplot2_theme()
              )
              p6 <- refine_ggplot2(
                p = p6,
                gridline = FALSE,
                ggplot2_theme = ggplot2_theme()
              )
              incProgress(7 / 8, detail = paste("Done"))
              gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 1)
            }
        })
      }) # isolate
    })

    # density plots using ggplot2
    output$genePlot2 <- renderPlot(
      {
        req(ggplot2_theme())
        req(sidebar_values$selectOrg())
        req(gene_density_plot())

        gene_density_plot()
      },
      width = 600,
      height = 2400
    )

    mod_download_images_server(
      "download_gene_plot_dist",
      filename = "gene_plot_dist",
      figure = reactive({
        ggpubr::as_ggplot(gene_density_plot())
      }),
      width = 6,
      height = 24
    )

  })
}
