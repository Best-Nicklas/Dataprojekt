#' Predict Status for new data with pre-existing model
#' 
#' This function is used to predict the status of one or more new persons.
#' 
#' @param rds.obj A rds object containing an FBM.code256.
#' @param model A model to use for predicting.
#' @param configs A data frame with configurations and their liabilities. Only used for LTFH.
#' @return A vector of 1's and 0's indicating predicted status.
#' @export
#' 


Predict_status <- function(rds.obj, model, configs = NULL, prevalence = 0.05){
  pred_value <- bigsnpr::snp_PRS(G = rds.obj$genotypes, betas.keep = model$Regression$estim, lpS.keep = -log10(model$Regression$p.value), thr.list = model$Pvalue)
  normalized_pred_value <- (pred_value - mean(pred_value))/sd(pred_value)
  if (is.null(configs)){
    predicted_status <- as.vector((normalized_pred_value > qnorm(prevalence, lower.tail = FALSE)) + 0)
  }
  else {
    configs <- configs[order(configs$Liability),]
    intervals <- findInterval(normalized_pred_value, configs$Liability)
    intervals[intervals == 0] <- 1
    predicted_status <- as.numeric(substr(configs$Config[intervals],1,1))
  }
  return(predicted_status)
}
