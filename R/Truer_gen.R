#' Description - 
#' 
#' @param p1
#' @param p2
#' @param seed
#' @return 
#' @example 
#' 

truer_gen = function(p1, p2, seed = NULL){
  if (!is.null(seed)) {set.seed(seed)}
  ph1 = which(p1 == 1, arr.ind = T) %>% as_tibble()
  ph2 = which(p2 == 1, arr.ind = T) %>% as_tibble()
  ind_11 = bind_rows(inner_join(ph1, ph2, by = c("row", "col")), 
                     inner_join(ph1, ph2, by = c("row", "col"))) %>% distinct()
  
  temp <- (p1 + p2) / 2
  temp[temp == 0.5] <- sample(0:1,length(temp[temp == 0.5]), replace = TRUE)
  temp[temp == 1.5] <- sample(1:2,length(temp[temp == 1.5]), replace = TRUE)
  temp[as.matrix(ind_11)] <- sample(0:2, size = nrow(ind_11), replace = T)  
  return(temp)
}