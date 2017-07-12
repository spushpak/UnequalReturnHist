constructRiskStats <- function(block.dat){
  risk.mat <- matrix(0, nrow=6, ncol=ncol(block.dat))
  row.names(risk.mat) <- c("Skewness", "Kurtosis", "Mean", "Volatility", 
                           "Sharpe Ratio", "Expected Shortfall")
  colnames(risk.mat) <- colnames(block.dat)
                     
  # Compute the risk and performance measures
  risk.mat["Skewness", ] <- moments::skewness(block.dat)
  risk.mat["Kurtosis", ] <- moments::kurtosis(block.dat)
  risk.mat["Mean", ] <- colMeans(block.dat)
  risk.mat["Volatility", ] <- apply(block.dat, 2, sd)
  risk.mat["Sharpe Ratio", ] <- risk.mat["Mean", ] / risk.mat["Volatility", ]
  #risk.mat["Sharpe Ratio", ] <- apply(block.dat, 2, SharpeRatio, FUN="StdDev")
  risk.mat["Expected Shortfall", ] <- apply(block.dat, 2, expectedShortfall)
  
  return(risk.mat)
}
