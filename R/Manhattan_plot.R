#' Description - Manhattan_plot
#' 
#' @param x
#' @return 
#' @example 
#' 

Manhattan_plot <- function(x) {
  ggplot(x, aes(x=1:length(pval), y=-log10(pval), size=-log10(pval))) + 
    geom_point(color="blue") + 
    ylim(0,15) +
    geom_hline(yintercept=-log10(5e-7), linetype=2) + xlab("SNP") + 
    ylab("-log10(P-value)")  
}