#' @title Constructs the new dataset by bootstrapping the residuals following Page(2013)
#' 
#' @description Constructs the new dataset by bootstrapping the residuals following Page(2013)
#' 
#' 
#' @importFrom moments skewness kurtosis
#' 
#' 
#' @param fitted.xts is the fitted data after regression.
#' @param resid.mat is the matrix containg residuals from linear regression.
#' @param miss.hist.var is the vector contating the names of the missing 
#' history assets.
#' @param miss.itr is the sorted vector containg the count of missing observations 
#' from the smallest to the largest.
#' @param na_count is the matrix containg the count of missing observations.
#' @param sReps is a logical flag, to indicate whether the risk measures  
#' need to be computed for all the replicates.
#' 
#' @return
#' \item{risk.metrics}{A list containing a matrix of risk measures for each replicate.}
#' \item{risk.vals}{A matrix of risk measures for the whole combined backfilled dataset.}
#' 
#'
#'    
#' @author Pushpak Sarkar
#' 
#' 
#' 
#' 
#' 
#' 
#' 

############################ Page MI backfill ###############################
constructPageData <- function(fitted.xts, resid.mat, miss.hist.var, miss.itr, na_count, sReps){
  
  # No. of bootstrap samples
  #M = 100000
  M = 100
  
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
      
      boot.resid <- as.matrix(resid.mat[sample(sample.indx, size=num.miss, 
                                               replace=TRUE), short.hist.var])
      colnames(boot.resid) <- short.hist.var
      fitted.xts[miss.val, short.hist.var] <- fitted.xts[miss.val, short.hist.var] + boot.resid[, short.hist.var]
    }
    
    risk.mat["Skewness", ] <- moments::skewness(fitted.xts[, miss.hist.var])
    risk.mat["Kurtosis", ] <- moments::kurtosis(fitted.xts[, miss.hist.var])
    risk.mat["Mean", ] <- colMeans(fitted.xts[, miss.hist.var])
    risk.mat["Volatility", ] <- apply(fitted.xts[, miss.hist.var], 2, sd)
    #risk.mat["Sharpe Ratio", ] <- SharpeRatio(fitted.xts[, miss.hist.var], FUN="StdDev")
    risk.mat["Sharpe Ratio", ] <- risk.mat["Mean", ]/risk.mat["Volatility", ]
    #risk.mat["Expected Shortfall", ] <- ES(fitted.xts[, miss.hist.var], p=.95, method="historical")
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
  
  if(sReps==TRUE) 
    return(risk.metrics) 
  else 
    return(risk.vals)
}
