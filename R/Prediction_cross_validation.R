#' Perform Prediction Cross-validation to determine best model
#' 
#' This function is used to calculate predictive powers of different models at different thresholds.
#' 
#' @param rds.obj A .rds file with an FBM.code256 and accompanying FAM and MAP.
#' @param k Number of folds to be used in cross-validation. The number of rows in the FBM must be at least twice as large as k. Highly recommended to choose k to be at most ~1% of the number of rows, unless working with a very small dataset, as errors may occur.
#' @param threshold Vector of P-values to be used in thresholding. Default does not use thresholding.
#' @param method Method to use for prediction. Possible methods are "GWAS", "GWAX", "LTFH". Default is "GWAS".
#' @param liabilities Vector of liabilities used for prediction with "LTFH" method. If not specified, uses "GWAS" method instead.
#' @return A list with 2 values: a tibble with average and best scores for each threshold, and a data.frame with the best model, fitted values, residuals, best p-value and its R^2.
#' @export
#' 


Prediction_cross_validation <- function(rds.obj, k, threshold = 1, method = "GWAS", liabilities = rds.obj$FAM$Status) {
  if (k < 1) stop("k must be positive")
  if (k%%1 != 0) stop("k must be an integer")
  n <- nrow(rds.obj$genotypes)
  if (n < k*2) stop("number of rows must be at least twice as large as k")
  block_size <- n%/%k
  bestest_score <- 0
  bestest_model <- NULL
  bestest_pval <- 0
  bestest_block_start <- 0
  bestest_block_end <- 0
  best_scores <- numeric(length(threshold))
  avg_scores <- numeric(length(threshold))
  # Get proxy status from GWAX for ease
  if(method == "GWAX"){
    liabilities <- GWAX(rds.obj)$Proxy_Status
  }
  
  for (i in 1:length(threshold)){
    scores <- numeric(k)
    best_score <- 0
    for (j in 0:(k-1)){
      #Create blocks
      block_start <- j*block_size + 1
      if (j == k-1){
        block_end <- n} else {
          block_end <- (j+1)*block_size}
      
      #Train on k-1 folds
      if(method == "GWAS"){
        regr <- GWAS(rds.obj, rds.obj$FAM$Status, include = c(1:n)[-(block_start:block_end)])
        regr$estim[is.nan(regr$estim)] <- 0
        regr$p.value[is.nan(regr$p.value)] <- 1
      }
      else if(method == "GWAX"){
        regr <- GWAS(rds.obj, liabilities, include = c(1:n)[-(block_start:block_end)])
        regr$estim[is.nan(regr$estim)] <- 0
        regr$p.value[is.nan(regr$p.value)] <- 1
      }
      else if(method == "LTFH"){
        regr <- GWAS(rds.obj, liabilities, include = c(1:n)[-(block_start:block_end)])
        regr$estim[is.nan(regr$estim)] <- 0
        regr$p.value[is.nan(regr$p.value)] <- 1
      }
      
      #Calculate PRS on 1 fold
      PRS <- bigsnpr::snp_PRS(G = rds.obj$genotypes, betas.keep = regr$estim, ind.test = block_start:block_end, lpS.keep = -log10(regr$p.value), thr.list = -log10(threshold)[i])
      
      #Normalize PRS with mean = 0 and sd = 1
      PRS <- (PRS - mean(PRS))/sd(PRS)
      
      #Correlation between PRS and status
      score <- cor(PRS, rds.obj$FAM$Status[block_start:block_end])
      
      #Update best models, scores
      scores[j+1] <- score
      if(score > best_score){
        best_score <- score
        best_model <- regr
        best_block_start <- block_start
        best_block_end <- block_end
      }
    }
    #Update best models, scores
    best_scores[i] <- best_score[1,1]
    avg_scores[i] <- mean(scores)
    if(best_score > bestest_score){
      bestest_score <- best_score[1,1]
      bestest_model <- best_model
      bestest_pval <- threshold[i]
      bestest_block_start <- best_block_start
      bestest_block_end <- best_block_end
    }
  }
  results <- tibble::tibble(Pvalue = threshold, Average_Score = avg_scores, Best_Score = best_scores, R2 = best_scores^2)
  fittedvals <- bigsnpr::snp_PRS(G = rds.obj$genotypes, betas.keep = bestest_model$estim, ind.test = c(1:n)[-(bestest_block_start:bestest_block_end)], lpS.keep = -log10(bestest_model$p.value), thr.list = bestest_pval)
  resids <- liabilities[-(bestest_block_start:bestest_block_end)]-fittedvals
  return(list(Results = results, Best_Model = list(Regression = data.frame(bestest_model), Fittedvalues = as.vector(fittedvals), Residuals = as.vector(resids), Score = bestest_score, Pvalue = bestest_pval)))
}
