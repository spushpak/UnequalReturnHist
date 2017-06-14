#' @title Implements Multiple Imputation Method folowing Page(2013)
#' 
#' @description This function implements Multiple Imputation Method folowing 
#' Page(2013)
#' 
#' @importFrom xts xts
#' @importFrom moments skewness kurtosis
#' 
#' 
#' @param dat.mat is the returns data for multiple assets with unequal return history.
#' @param saveReps is a logical flag, to indicate whether the risk measures  
#' need to be computed for all the replicates. Default is FALSE.
#' 
#' @return
#' \item{risk.list}{A list or matrix containing risk measures for either the whole dataset or for each replicate.}
#'    
#' @author Pushpak Sarkar
#' 
#' 
#' 
#' 
#' 
#' 
#' @export
#' 


#### Do the regression and find fitted value for the missing returns
fitValue <- function(short.hist.var, long.hist.var, ret.dat, num.miss){
  
  regressor.list <- paste(long.hist.var, collapse = "+")
  
  est.coef <- matrix(0, nrow=length(long.hist.var)+1, ncol=length(short.hist.var),
                     dimnames=list(c("Intercept", long.hist.var), short.hist.var))
  
  num.miss <- num.miss
  len.diff <- nrow(ret.dat) - num.miss
  
  err.mat <- matrix(0, nrow=len.diff, ncol=length(short.hist.var))
  colnames(err.mat) <- short.hist.var
  
  for (j in short.hist.var) {
    reg.eqn <- paste(j, "~", regressor.list)
    reg <- lm(as.formula(reg.eqn), data=ret.dat)
    
    est.coef[ , j] <- as.matrix(coef(reg))
    err.mat[ , j] <- as.matrix(resid(reg))
  }
  
  miss.val <- which(is.na(ret.dat[, short.hist.var[1]]))
  
  for (k in 1:num.miss) {
    X.mat <- cbind(1, ret.dat[miss.val[k], long.hist.var])
    ret.dat[miss.val[k], short.hist.var] <- X.mat%*%est.coef
  }
  
  # Remove the long.hist.var
  ret.dat <- ret.dat[, ! names(ret.dat) %in% long.hist.var, drop=FALSE]
  
  # Create a list of fitted values, error matrix and missing values
  reg.dat <- list(fitted.dat=ret.dat, err.mat=err.mat, miss.val=miss.val)
  
  return(reg.dat)
}


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




uneqhistPage <- function(dat.mat, saveReps=FALSE){

# Convert the data to xts object (probably it's not necessary)
dat.xts <- xts(dat.mat[, -1], order.by = as.Date(dat.mat[, 1], "%m/%d/%Y"))

# Keep a copy of original data - will be updated with fitted values for the
# missing portion
fitted.xts <- dat.xts 

# Find which columns have 'NA' values; so these columns have shorter histtory
miss.hist.var <- colnames(dat.xts)[apply(dat.xts, 2, anyNA)]

# Find which column has full history
long.hist.var <- setdiff(colnames(dat.xts), miss.hist.var)

# Find how many observatins are missing for each short history columns
na_count <- apply(dat.xts, 2, function(x) sum(is.na(x)))
na_count <- as.matrix(na_count)
colnames(na_count) <- "num.miss"

# Alternative way to find the full history column
#long.hist.var <- rownames(na_count)[na_count[, "num.miss"]==0]

# Drop full history group from 'na_count' as 'NA' count for this group is zero
na_count <- na_count[miss.hist.var, "num.miss", drop=F]

# Sort count of missing obs from smallest to largest
miss.itr <- sort(unique(na_count[,"num.miss"]))

# Full length of the dataset
full.length <- nrow(dat.xts)

resid.mat <- matrix(NA, nrow = full.length, ncol = length(miss.hist.var))
colnames(resid.mat) <- miss.hist.var

# Start with the short history columns which have the smallest number of 
# observations missing i.e. this group is the longest short history group. Then 
# move on to the next shorter history group and so on.
for (i in miss.itr) {
  short.hist.var <- rownames(na_count)[na_count[, "num.miss"]==i]
  temp.dat <- dat.xts[, c(long.hist.var, short.hist.var), drop=F]
  num.miss <- i
  reg.dat <- fitValue(short.hist.var, long.hist.var, temp.dat, num.miss)

  fitted.xts[reg.dat$miss.val, short.hist.var] <- reg.dat$fitted.dat[reg.dat$miss.val, short.hist.var]
  resid.mat[(length(reg.dat$miss.val)+1):nrow(resid.mat), short.hist.var] <- reg.dat$err.mat[, short.hist.var]
} ############# end of for loop (for each short history group)

risk.list <- constructPageData(fitted.xts, resid.mat, miss.hist.var, miss.itr, 
                               na_count, sReps=saveReps)

return(round(risk.list, digits = 3))

}
