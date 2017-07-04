############################ Page MI backfill ###############################
constructMIData <- function(fitted.xts, resid.mat, miss.hist.var, miss.itr, na_count, saveReps, M){
  
  # No. of bootstrap samples
  M <- M
  
  
  # Create an empty list
  risk.metrics <- vector("list", M)
  
  # Create an empty matrix to hold the risk measures
  risk.vals <- matrix(0, nrow=6, ncol=length(miss.hist.var),
                      dimnames=list(c("Skewness", "Kurtosis", "Mean", "Volatility", "Sharpe Ratio", "Expected Shortfall"), miss.hist.var))
  
  
  temp.skew <- matrix(0, nrow = M, ncol=length(miss.hist.var))
  colnames(temp.skew) <- miss.hist.var
  temp.kurt <- matrix(0, nrow = M, ncol=length(miss.hist.var))
  colnames(temp.kurt) <- miss.hist.var
  temp.mean <- matrix(0, nrow = M, ncol=length(miss.hist.var))
  colnames(temp.mean) <- miss.hist.var
  temp.vol <- matrix(0, nrow = M, ncol=length(miss.hist.var))
  colnames(temp.vol) <- miss.hist.var
  temp.mean <- matrix(0, nrow = M, ncol=length(miss.hist.var))
  colnames(temp.mean) <- miss.hist.var
  temp.SR <- matrix(0, nrow = M, ncol=length(miss.hist.var))
  colnames(temp.SR) <- miss.hist.var
  temp.ES <- matrix(0, nrow = M, ncol=length(miss.hist.var))
  colnames(temp.ES) <- miss.hist.var
  
  
  # Bootstrap sampling - draw a random sample of size k (no. of missing obs) from 
  # (n-k) residuals with replacement. These k residuals are added to the k fitted 
  # values.
  for(i in 1:M){
    risk.mat <- matrix(0, nrow=6, ncol=length(miss.hist.var), 
                       dimnames=list(c("Skewness", "Kurtosis", "Mean", "Volatility", "Sharpe Ratio", "Expected Shortfall"), miss.hist.var))
    
    for (j in miss.itr) {
      short.hist.var <- rownames(na_count)[na_count[, "num.miss"]==j]
      miss.val <- which(is.na(resid.mat[, short.hist.var[1]]))
      num.miss <- length(miss.val)
      sample.indx <- which(!is.na(resid.mat[, short.hist.var[1]]))
      
      boot.resid <- as.matrix(resid.mat[sample(sample.indx, size=num.miss, replace=TRUE), short.hist.var])
      colnames(boot.resid) <- short.hist.var
      boot.resid <- as.xts(boot.resid, order.by = index(fitted.xts)[miss.val])
      fitted.xts[miss.val, short.hist.var] <- fitted.xts[miss.val, short.hist.var, drop=F] + boot.resid[, short.hist.var, drop=F]
    }
    
    risk.mat["Skewness", ] <- moments::skewness(fitted.xts[, miss.hist.var, drop=F])
    risk.mat["Kurtosis", ] <- moments::kurtosis(fitted.xts[, miss.hist.var, drop=F])
    risk.mat["Mean", ] <- colMeans(fitted.xts[, miss.hist.var, drop=F])
    risk.mat["Volatility", ] <- apply(fitted.xts[, miss.hist.var, drop=F], 2, sd)
    risk.mat["Sharpe Ratio", ] <- risk.mat["Mean", ] / risk.mat["Volatility", ]
    #risk.mat["Sharpe Ratio", ] <- apply(fitted.xts[, miss.hist.var, drop=F], 2, SharpeRatio, FUN="StdDev")
    risk.mat["Expected Shortfall", ] <- apply(fitted.xts[, miss.hist.var, drop=F], 2, expectedShortfall)
    
    risk.metrics[[i]] <- risk.mat
    
    temp.skew[i, ] <- risk.mat["Skewness", ]  
    temp.kurt[i, ] <- risk.mat["Kurtosis", ]
    temp.mean[i, ] <- risk.mat["Mean", ]
    temp.vol[i, ] <- risk.mat["Volatility", ]
    temp.SR[i, ] <- risk.mat["Sharpe Ratio", ]
    temp.ES[i, ] <- risk.mat["Expected Shortfall", ]
  }
  
  risk.vals["Skewness", ] <- colMeans(temp.skew)
  risk.vals["Kurtosis", ] <- colMeans(temp.kurt)
  risk.vals["Mean", ] <- colMeans(temp.mean)
  risk.vals["Volatility", ] <- colMeans(temp.vol)
  risk.vals["Sharpe Ratio", ] <- colMeans(temp.SR)
  risk.vals["Expected Shortfall", ] <- colMeans(temp.ES)
  
  if(saveReps==TRUE) 
    return(risk.metrics) 
  else 
    return(risk.vals)
}