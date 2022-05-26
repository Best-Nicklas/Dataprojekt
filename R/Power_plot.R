#' Create a power plot from analysis
#' 
#' This function plots the power of different analysis methods.
#' 
#' @param x placeholder
#' @return placeholder
#' @export 
Power_plot <- function(x) {
  powerplot <- x %>% dplyr::mutate (causal_snp = p.value < 5e-3) %>%
    dplyr::arrange(abs(estim)) %>%
    dplyr::mutate (cpower = cumsum(causal_snp))
  
  return(powerplot %>% ggplot2::ggplot(aes(x = estim, y = cpower)) + ggplot2::geom_line())
}