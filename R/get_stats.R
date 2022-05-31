#' Get statistics on how many SNPs have been correctly identified
#' 
#' This functions gives statistics on the count and frequency of false negatives, true negatives,
#' false positives and true positives in analysed data
#' 
#' @param pvalues A vector of the pvalues calculated from analysis.
#' @param true_causal A vector indicating causal SNPs (1 where causal, 0 where not causal).
#' @param a An integer specifying significance level and thereby the cutoff for a SNP being causal.
#' @return A dataframe containing the statistics
#' @export
#'

get_stats <- function(pvalues,  true_causal, a) {
  if (length(pvalues) != length(true_causal)) stop("Length of pvalues vector and true_causal vector must be the same")
  if (a < 0) stop("threshold pvalue must be non-negative")
  
  # find predicted causal SNPs 
  pred_causal <- (pvalues < a) + 0

  # calculate true negatives, false negatives, true positives and false positives
  false_negatives <- sum(pred_causal == 0 & true_causal == 1)
  true_negatives <- sum(pred_causal == 0 & true_causal == 0)
  false_positives <- sum(pred_causal == 1 & true_causal == 0)
  true_positives <- sum(pred_causal == 1 & true_causal == 1)
  total <- length(pred_causal)
  
  # Collect stats in dataframe
  stats <- c(false_negatives, true_negatives, false_positives, true_positives, total)
  data <- matrix(c(stats, stats / total), ncol = 5, byrow = T)
  rownames(data) <- c("Count", "Frequency")
  colnames(data) <- c("false_negatives", "true_negatives", "false_positives", "true_positives", "total")
  data.frame(data)
}