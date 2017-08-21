constructRiskStats <- function(block_dat){
  risk_mat <- matrix(0, nrow=6, ncol=ncol(block_dat))
  row.names(risk_mat) <- c("Mean", "Volatility", "Skewness", "Kurtosis",  
                           "Expected Shortfall", "Sharpe Ratio")
  colnames(risk_mat) <- colnames(block_dat)
                     
  # Compute the risk and performance measures
  risk_mat["Mean", ] <- colMeans(block_dat)
  risk_mat["Volatility", ] <- apply(block_dat, 2, sd)
  risk_mat["Skewness", ] <- moments::skewness(block_dat)
  risk_mat["Kurtosis", ] <- moments::kurtosis(block_dat)
  risk_mat["Expected Shortfall", ] <- apply(block_dat, 2, expectedShortfall)
  risk_mat["Sharpe Ratio", ] <- risk_mat["Mean", ] / risk_mat["Volatility", ]
  #risk_mat["Sharpe Ratio", ] <- apply(block_dat, 2, SharpeRatio, FUN="StdDev")
  
  return(risk_mat)
}
