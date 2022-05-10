#' Description - 
#' 
#' @param p1
#' @param p2
#' @param disease
#' @param path
#' @param n_blocks
#' @param seed
#' @return


simulate_child_with_sibs <- function(p1, p2, disease, path, n_blocks = 20, seed = NULL) {
  #Simulates children from p1 and p2
  child <- Simulate_children(p1, p2, disease, path)
  
  rows <- nrow(child$genotypes)
  cols <- ncol(child$genotypes)
  
  threshold <- qnorm(1 - disease$PREVALENS)
  h2 <- disease$H2
  
  # Calculates the family information - liability,status, config etc.
  child$FAM$n_sib = 2
  
  #Creates parents according to the number of sibs to generate children and their g_liabs from
  child$FAM$sibs_g_liabs = purrr::map2(1:rows, child$FAM$n_sib, .f = ~ {sibs_liab_calc(.y, 
                                                                                       matrix(rep(p1$genotypes[.x, ], .y), .y, cols, byrow = T),
                                                                                       matrix(rep(p2$genotypes[.x, ], .y), .y, cols, byrow = T), 
                                                                                       disease)})
  #calculates full liabilities
  child$FAM$sibs_full_liabs = purrr::map(child$FAM$sibs_g_liabs, .f = ~ 
                                           {if ( is.null(.x)) NULL 
                                             else .x + rnorm(length(.x), 0, sqrt(1 - h2))})
  #calculates status
  child$FAM$sibs_status = purrr::map(child$FAM$sibs_full_liabs, .f = ~ 
                                       {if ( is.null(.x)) NULL
                                         else(.x > threshold) + 0})
  
  
  
  return(child)
  
}