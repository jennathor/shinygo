#' Genome sliding window enrichment analysis
#'
#' Scans the genome with a sliding window and applies a hypergeometric test
#' to determine if query genes are significantly overrepresented in each window.
#' Returns FDR-corrected enriched regions.
#'
#' @param x0        Data frame with columns: start_position (Mbp), chNum,
#'                  chromosome_name, Fold (1 = list gene, 0 = background)
#' @param listN     Number of genes in the query list
#' @param totalN    Total number of genes (list + background)
#' @param windowSize  Window size in Mbp
#' @param steps     Number of steps per window
#' @param pvalCutoff  FDR cutoff for filtering enriched windows
#' @param chD       Y-axis spacing between chromosomes (default 30)
#'
#' @return Data frame of enriched regions with columns chNum, x, n, k, pval, y
genome_sliding_window <- function(x0, listN, totalN, windowSize, steps, pvalCutoff, chD = 30) {
  for (i in 0:(steps - 1)) {
    # step size is  windowSize/steps
    # If windowSize=10 and steps = 2; then step size is 5Mb
    # 1.3 becomes 5, 11.2 -> 15 for step 1
    # 1.3 -> -5
    x0$x <- (floor((x0$start_position - i * windowSize / steps) / windowSize)
    + 0.5 + i / steps) * windowSize

    movingAverage1 <- x0 %>%
      select(chNum, x, Fold) %>%
      filter(x >= 0) %>% # beginning bin can be negative for first bin in the 2nd step
      group_by(chNum, x) %>%
      summarize(n = n(), k = sum(Fold)) %>%
      filter(k > 0) %>%
      filter(k / n > listN / totalN) %>%
      mutate(pval = phyper(k - 1,
        n,
        totalN - n,
        listN,
        lower.tail = FALSE
      ))

    if (i == 0) {
      movingAverage <- movingAverage1
    } else {
      movingAverage <- rbind(movingAverage, movingAverage1)
    }
  }

  # translate fold to y coordinates
  movingAverage <- movingAverage %>%
    filter(n >= 3) %>%
    mutate(pval = p.adjust(pval, method = "fdr")) %>%
    filter(pval < pvalCutoff) %>%
    mutate(y = chNum * chD - 4)

  return(movingAverage)
}
