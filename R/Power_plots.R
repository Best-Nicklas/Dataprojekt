#' Create power plots from GWAS, GWAX and LTFH data
#' 
#' This function is used to calculate and plot the power for each of the models GWAS, GWAX and LTFH.
#' This function is only design to plot all three plots against each other to get a visual comparison of the power plots for each model. 
#' 
#' @param Gwas_data A list with exactly three entrances, with the GWAS-values gained from from each model (GWAS, GWAX and LTFH), where the data from GWAS has to be the first entrance, the data from Gwax has to be the second entrance and the data from LTFH has to be on the last entrance. 
#' @param a significance level used to determine causal SNPs
#' @return Ggplot with power plots from GWAS, GWAX and LTFH, to get a visual comparison of the prediction power for each model. 
#' @export
#' 

Power_plots <- function(Gwas_data, a){
   
   if (!(length(Gwas_data) == 3)) stop("Needs Gwas_data from all three analysis in list format")
  
    Gwas <- Gwas_data[[1]]
    Gwax <- Gwas_data[[2]]
    Ltfh <- Gwas_data[[3]]
    
    # Calculate correct input data, causal SNPs
    T1 <- Gwas %>% dplyr::mutate(causal_snp = p.value < a) %>%
      dplyr::arrange(abs(estim)) %>%
      dplyr::mutate(cpower = cumsum(causal_snp)) %>%
      dplyr::mutate(Method = "GWAS")
    
    T2 <- Gwax %>% dplyr::mutate(causal_snp = p.value < a) %>%
      dplyr::arrange(abs(estim)) %>%
      dplyr::mutate (cpower = cumsum(causal_snp)) %>%
      dplyr::mutate (Method = "GWAX")
    
    T3 <- Ltfh %>% dplyr::mutate(causal_snp = p.value < a) %>%
      dplyr::arrange(abs(estim)) %>%
      dplyr::mutate(cpower = cumsum(causal_snp)) %>%
      dplyr::mutate(Method = "LTFH")
    
    # Plot powers
    ggplot2::ggplot(rbind(T1, T2, T3)) + 
      ggplot2::geom_line(ggplot2::aes(x = estim, y = cpower, group = Method, colour = Method)) +
      ggplot2::xlab("Estimated Effect Size") +
      ggplot2::ylab("Cumulative Power (Cumulative number of SNPs found Causal)") 
}
