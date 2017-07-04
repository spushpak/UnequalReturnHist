############################ Page MI backfill ###############################
constructMIData <- function(fitted.xts, resid.mat, miss.hist.var, miss.itr, na_count, saveReps, M){
  
  # No. of bootstrap samples
  M <- M
  
  # Create an empty list
  risk.metrics <- vector("list", M)
  
  # Create an empty matrix to hold the risk measures
  risk.vals <- matrix(0, nrow=6, ncol=length(miss.hist.var),
                      dimnames=list(c("Skewness", "Kurtosis", "Mean", "Volatility", "Sharpe Ratio", "Expected Shortfall"), miss.hist.var))
  
  new.dat <- fitted.xts
  
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
      
      r.sample <- sample(sample.indx, size=num.miss, replace=TRUE)
      boot.resid <- resid.mat[r.sample, short.hist.var, drop=F]
      colnames(boot.resid) <- short.hist.var
      boot.resid <- as.xts(boot.resid, order.by = index(fitted.xts)[miss.val])
      new.dat[miss.val, short.hist.var] <- fitted.xts[miss.val, short.hist.var, drop=F] + boot.resid[, short.hist.var, drop=F]
    }
    
    block.dat <- new.dat[, miss.hist.var, drop=F]
      
    risk.mat["Skewness", ] <- moments::skewness(block.dat)
    risk.mat["Kurtosis", ] <- moments::kurtosis(block.dat)
    risk.mat["Mean", ] <- colMeans(block.dat)
    risk.mat["Volatility", ] <- apply(block.dat, 2, sd)
    risk.mat["Sharpe Ratio", ] <- risk.mat["Mean", ] / risk.mat["Volatility", ]
    #risk.mat["Sharpe Ratio", ] <- apply(block.dat, 2, SharpeRatio, FUN="StdDev")
    risk.mat["Expected Shortfall", ] <- apply(block.dat, 2, expectedShortfall)
    
    risk.metrics[[i]] <- risk.mat
  }
  
  risk.vals["Skewness", ] <- colMeans(do.call("rbind", lapply(risk.metrics, "[", "Skewness", TRUE)))
  risk.vals["Kurtosis", ] <- colMeans(do.call("rbind", lapply(risk.metrics, "[", "Kurtosis", TRUE)))
  risk.vals["Mean", ] <- colMeans(do.call("rbind", lapply(risk.metrics, "[", "Mean", TRUE)))
  risk.vals["Volatility", ] <- colMeans(do.call("rbind", lapply(risk.metrics, "[", "Volatility", TRUE)))
  risk.vals["Sharpe Ratio", ] <- colMeans(do.call("rbind", lapply(risk.metrics, "[", "Sharpe Ratio", TRUE)))
  risk.vals["Expected Shortfall", ] <- colMeans(do.call("rbind", lapply(risk.metrics, "[", "Expected Shortfall", TRUE)))
  

  if(saveReps==TRUE) 
    return(risk.metrics) 
  else 
    return(risk.vals)
}
