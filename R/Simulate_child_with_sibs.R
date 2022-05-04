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
  child <- Simulate_children(p1, p2, disease, path)
  rows <- nrow(child$genotypes)
  cols <- ncol(child$genotypes)
  
  threshold <- qnorm(1 - disease$PREVALENS)
  h2 <- disease$H2
  
  
  family_info <- tibble(id = 1:rows, 
                        
                        n_sib = 2,
                        
                        sibs_g_liabs = purrr::map2(id, n_sib, .f = ~ {sibs_liab_calc(.y, 
                                                                                     matrix(rep(p1$genotypes[.x, ], .y), .y, cols, byrow = T),
                                                                                     matrix(rep(p2$genotypes[.x, ], .y), .y, cols, byrow = T), disease)}),
                        
                        sibs_full_liabs = purrr::map(sibs_g_liabs , .f = ~ 
                                                       {if ( is.null(.x)) NULL 
                                                         else .x + rnorm(length(.x), 0, sqrt(1 - h2))}),
                        
                        sibs_status = purrr::map(sibs_full_liabs, .f = ~ 
                                                   {if ( is.null(.x)) NULL
                                                     else(.x > threshold) + 0}),
                        
                        p1_status = p1$FAM$Status,
                        
                        p2_status = p2$FAM$Status,
                        
                        config = mapply(function(s_child, s_p1,
                                                 s_p2, s_sibs)
                        {paste(toString(s_child),
                               toString(max(s_p1, s_p2)),
                               toString(min(s_p1, s_p2)),
                               strrep(1, sum(s_sibs)),
                               strrep(0,length(s_sibs) - sum(s_sibs)),
                               sep = "")},
                        child$FAM$Status, p1_status, p2_status, sibs_status))
  
  child$family_info <- family_info
  return(child)
  
}