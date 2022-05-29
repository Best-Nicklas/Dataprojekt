#' Perform Prediction Cross-validation to determine best model
#' 
#' This function is used to calculate predictive powers of different models.
#' 
#' @param rds.obj A .rds file with an FBM.code256 and accompanying FAM and MAP.
#' @param model A model to use for predicting.
#' @param configs Data frame with configurations and their liabilities. Only used for LTFH.
#' @return A list
#' @export
#' 


Predict_status <- function(rds.obj, model, configs = NULL){
  pred_value <- bigsnpr::snp_PRS(G = rds.obj$genotypes, betas.keep = model$Regression$estim, lpS.keep = -log10(model$Regression$p.value), thr.list = model$Pvalue)
  normalized_pred_value <- (pred_value - mean(pred_value))/sd(pred_value)
  
  #lwr <- normalized_pred_value
  #upr
  return(test)
}
