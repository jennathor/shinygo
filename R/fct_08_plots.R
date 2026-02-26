# find peak values in density plots
# for adding annotation texts
# http://ianmadd.github.io/pages/PeakDensityDistribution.html
densMode <- function(x) {
  td <- density(x, na.rm = TRUE)
  maxDens <- which.max(td$y)
  list(x = td$x[maxDens], y = td$y[maxDens])
}

#' Change ggplot2 plots
#'
#'
#' @param p ggplot2 object
#' @param gridline TRUE of FALSE
#'
#' @export
#' @return ggplot2 object
refine_ggplot2 <- function(p, gridline, ggplot2_theme = "light") {

  # apply theme based on selection
  p <- switch(ggplot2_theme,
    "linedraw" = p + ggplot2::theme_linedraw(),
    "classic" = p + ggplot2::theme_classic(),
    "gray" = p + ggplot2::theme_gray(),
    "light" = p + ggplot2::theme_light(),
    "dark" = p + ggplot2::theme_dark(),
    "bw" = p + ggplot2::theme_bw(),
    p # default, no change
  )

  if (ggplot2_theme != "Add grid") { # keep grid
    if (!gridline) { # by default it has gridlines
      p <- p +
        ggplot2::theme(panel.grid = ggplot2::element_blank())
    }
  }

  return(p)
}

# generates a fake ggplot2, with some message like: "Not available."
fake_plot <- function(some_text) {
  p <- ggplot2::ggplot() +
    geom_point() +
    xlim(-10, 10) +
    ylim(-10, 10) +
    annotate("text",
      x = 0,
      y = 0,
      label = some_text
    ) +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
  return(p)
}

# 0.000234   <- 2.3E-4 *
mark_significance <- function(Pval, PvalGeneInfo2, PvalGeneInfo1, PvalGeneInfo) {
  sig <- paste("P=", formatC(Pval, digits = 2, format = "G"), sep = "")
  if (Pval < PvalGeneInfo2) {
    sig <- paste(sig, "***")
  } else
  if (Pval < PvalGeneInfo1) {
    sig <- paste(sig, "**")
  } else
  if (Pval < PvalGeneInfo) sig <- paste(sig, "*")
  return(sig)
}
