#' Make a Manhattan plot from analysis data 
#' 
#' This function produces a Manhattan plot from p-values found in an association analysis. 
#' 
#' @param x List object containing p-values from an analysis. 
#' @return A Manhattan plot with each regressor (SNP) on the x-axis and their -log10(P-values) on the y-axis. 
#' @export
#' 
Manhattan_plot <- function(x) {
  ggplot2::ggplot(x, aes(x=1:length(p.value), y=-log10(p.value), size=-log10(p.value))) + 
    ggplot2::geom_point(color="blue") + 
    ggplot2::ylim(0,15) +
    ggplot2::geom_hline(yintercept=-log10(5e-7), linetype=2) +
    ggplot2::xlab("SNP") + 
    ggplot2::ylab("-log10(P-value)")  
}