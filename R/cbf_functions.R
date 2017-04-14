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

#################################################################################
### Form the new dataset based on the num of residuals

uneqhistCBF <- function(reg.dat, new.dat){
  
  num.resid <- nrow(reg.dat$err.mat)
  miss.val <- reg.dat$miss.val
  err.mat <- reg.dat$err.mat
  num.miss <- length(miss.val)
  
  # Stack the original block 'new.dat' equal to number of residual times
  temp.newdat <- matrix(rep(t(coredata(new.dat)), num.resid), ncol= ncol(new.dat), 
                    byrow=TRUE)
  colnames(temp.newdat) <- colnames(new.dat)
  
  # Number of times the 'fitted.dat' i.e. short.hist.var should be repeated
  rep.num <- nrow(temp.newdat)/full.length
  temp.fitdat <- matrix(rep(t(coredata(reg.dat$fitted.dat)), rep.num), 
                        ncol= ncol(reg.dat$fitted.dat), byrow=TRUE)
  colnames(temp.fitdat) <- colnames(reg.dat$fitted.dat)
  
  # Merge the existing new.dat and fitted.dat
  new.dat <- cbind(temp.newdat, temp.fitdat)
  
  # Crate the start index of each repeating block
  rep.indx <- seq(from = miss.val[1]-1, by = nrow(new.dat)/num.resid, 
                  length.out = num.resid)
  
  # Add the  respective residual to the fitted values of short history vars
  for (m in 1:length(rep.indx)) {
    for (n in 1:num.miss) {
      new.dat[rep.indx[m]+n, short.hist.var] <- new.dat[rep.indx[m]+n, short.hist.var] + err.mat[m, short.hist.var]
    }
  }
  return(new.dat)
}


#############################################################################
############################ Page MI backfill ###############################
uneqhistPage <- function(fitted.xts, resid.mat, miss.hist.var, miss.itr, na_count, saveReps=FALSE){
  
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

  if(saveReps==TRUE) 
    return(risk.metrics) 
  else 
    return(risk.vals)
}


#############################################################################
############################ Page (2013) ###################################
page_2013 <- function(reg.dat, short.hist.var){
  
  num.resid <- nrow(reg.dat$err.mat)
  miss.val <- reg.dat$miss.val
  err.mat <- reg.dat$err.mat
  num.miss <- length(miss.val)
  fitted.dat <- reg.dat$fitted.dat
  
  # No. of bootstrap samples
  #M = 100000
  M = 100
  
  # Make matrices to hold the skewness and kurtosis for each bootstrap sample.
  skew.mat <- matrix(0, nrow = M, ncol = length(short.hist.var))
  colnames(skew.mat) <- short.hist.var
  
  kurt.mat <- matrix(0, nrow = M, ncol = length(short.hist.var))
  colnames(kurt.mat) <- short.hist.var
  
  
  # Bootstrap sampling - draw a random sample of size k (no. of missing obs) from 
  # (n-k) residuals with replacement. These k residuals are added to the k fitted 
  # values.
  for(i in 1:M){
    boot.resid <- as.matrix(err.mat[sample(nrow(err.mat), size=num.miss, replace=TRUE), ])
    colnames(boot.resid) <- short.hist.var
    
    fitted.dat[miss.val, short.hist.var] <- fitted.dat[miss.val, short.hist.var] + boot.resid[, short.hist.var]
    
    skew.mat[i, short.hist.var] <- skewness(fitted.dat[miss.val, short.hist.var])
    kurt.mat[i, short.hist.var] <- kurtosis(fitted.dat[miss.val, short.hist.var])
  }
  
  skew.kurt <- matrix(0, nrow = 2, ncol = length(short.hist.var), 
                      dimnames = list(c("skewness", "kurtosis"), short.hist.var))
  
  # Mean values of skewness and kurtosis are calculated and saved
  skew.kurt["skewness", ] <- apply(skew.mat, 2, mean)
  skew.kurt["kurtosis", ] <- apply(kurt.mat, 2, mean)
  
  return(skew.kurt)
}

################################################################################
##################### Function for descriptive statistics ######################
desc.stats <- function(dataset){
  
  library(moments)
  stats.vec <- c("Mean","Median","Maximum","Minimum","Std. Dev.","Skewness","Kurtosis")
  dstats <- matrix(0, nrow = length(stats.vec), ncol = ncol(dataset), 
                   dimnames = list(stats.vec, colnames(dataset)))
  
  # dstats[1, ] <- apply(dataset, 2, mean, na.rm=T)
  # dstats[2, ] <- apply(dataset, 2, median, na.rm=T)
  # dstats[3, ] <- apply(dataset, 2, max, na.rm=T)
  # dstats[4, ] <- apply(dataset, 2, min, na.rm=T)
  # dstats[5, ] <- apply(dataset, 2, sd, na.rm=T)
  dstats[6, ] <- apply(dataset, 2, skewness, na.rm=T)
  dstats[7, ] <- apply(dataset, 2, kurtosis, na.rm=T)
  
  return(dstats)
} 


#################################################################################
######## Skewness and Kurtosis calculation for missing portion of data ##########

calcMoM.missRet <- function(new.dat, full.length, miss.itr, na_count){
  
  # Create the start index of each repeating block in the new dataset. Each repeating
  # block has length equal to the row numbers of initial dataset.
  row.grp <- seq(1, nrow(new.dat), by=full.length)

  # Create a matrix for holding skewness and kurtosis of missing portion of the 
  # short history returns
  mom.mat <- matrix(0, nrow=nrow(na_count), ncol=2, 
                    dimnames=list(miss.hist.var, c("Skewness", "Kurtosis")))

  for (i in miss.itr) {
    # Pick the missing history group staring with the shortest
    short.hist.var <- rownames(na_count)[na_count[, "num.miss"]==i]
    num.miss <- i
  
    # Start from the first row; and select those many rows equal to the number of 
    # missing obs (i.e. num.miss) for this missing history group
    mom.dat <- new.dat[row.grp[1]:(row.grp[1]+num.miss-1), short.hist.var, drop=F]
  
    # Repeat the above step for each repaeting block 
    for (j in 2:length(row.grp)) {
      temp <- new.dat[row.grp[j]:(row.grp[j]+num.miss-1), short.hist.var, drop=F]
      mom.dat <- rbind(mom.dat, temp)
    }
  
    mom.mat[short.hist.var, "Skewness"] <- t(moments::skewness(mom.dat))
    mom.mat[short.hist.var, "Kurtosis"] <- t(moments::kurtosis(mom.dat))
  }

  return(mom.mat)
}




#################################################################################
#################################################################################
combined.backfill.sim <- function(dat.xts) {
  
  short.hist.var <- colnames(dat.xts)[apply(dat.xts, 2, anyNA)]
  long.hist.var <- setdiff(colnames(dat.xts), short.hist.var)
  
  miss.val <- which(is.na(dat.xts[, short.hist.var[1]]))
  num.miss <- length(miss.val)
  
  regressor.list <- paste(long.hist.var, collapse = "+")
  
  est.coef <- matrix(0, nrow=length(long.hist.var)+1, ncol=length(short.hist.var),
                     dimnames=list(c("Intercept", long.hist.var), short.hist.var))
  
  len.diff <- nrow(dat.xts) - num.miss
  err.mat <- matrix(0, nrow=len.diff, ncol=length(short.hist.var))
  colnames(err.mat) <- short.hist.var
  
  for (k in short.hist.var) {
    reg.eqn <- paste(k, "~", regressor.list)
    #print(reg.eqn)
    reg <- lm(as.formula(reg.eqn), data=dat.xts)
    est.coef[ , k] <- as.matrix(coef(reg))
    err.mat[ , k] <- as.matrix(resid(reg))
  }
  
  for (j in 1:num.miss) {
    X.mat <- cbind(1, dat.xts[miss.val[j], long.hist.var])
    dat.xts[miss.val[j], short.hist.var] <- X.mat%*%est.coef
  }
  
  # Save a copy of the updated dataset
  fitted.dat <- dat.xts
  
  ############################ Combined Backfill #########################
  # Stack the block (which had missing returns for short history asset) 
  # equal to number of rep.num
  rep.num <- 5000

  temp <- dat.xts[miss.val, short.hist.var]
  
  new.dat <- matrix(rep(temp, rep.num), ncol= 1, byrow=TRUE)

  # Crate the start index of each repeating block; each block has num.miss obs
  rep.indx <- seq(from = miss.val[1]-1, by = num.miss, length.out = rep.num)
  
  # Draw a random sample of 5000 residuals from 90,000 residuals 
  err.sample <- sample(as.numeric(err.mat), rep.num, replace=T)
  
  # Add the  respective residual to the fitted values of short history vars
  for (i in 1:rep.num) {
    for (j in 1:num.miss) {
      new.dat[rep.indx[i]+j, 1] <- new.dat[rep.indx[i]+j, 1] + err.sample[i]
    }
  }
  
  # This new.dat is the new long data series after combined backfill
  
  return(new.dat)
}

################################################################################
################################################################################
page.backfill.sim <- function(dat.xts){
  
  short.hist.var <- colnames(dat.xts)[apply(dat.xts, 2, anyNA)]
  long.hist.var <- setdiff(colnames(dat.xts), short.hist.var)
  
  miss.val <- which(is.na(dat.xts[, short.hist.var[1]]))
  num.miss <- length(miss.val)
  
  regressor.list <- paste(long.hist.var, collapse = "+")
  
  est.coef <- matrix(0, nrow=length(long.hist.var)+1, ncol=length(short.hist.var),
                     dimnames=list(c("Intercept", long.hist.var), short.hist.var))
  
  len.diff <- nrow(dat.xts) - num.miss
  err.mat <- matrix(0, nrow=len.diff, ncol=length(short.hist.var))
  colnames(err.mat) <- short.hist.var
  
  for (k in short.hist.var) {
    reg.eqn <- paste(k, "~", regressor.list)
    #print(reg.eqn)
    reg <- lm(as.formula(reg.eqn), data=dat.xts)
    est.coef[ , k] <- as.matrix(coef(reg))
    err.mat[ , k] <- as.matrix(resid(reg))
  }
  
  for (j in 1:num.miss) {
    X.mat <- cbind(1, dat.xts[miss.val[j], long.hist.var])
    dat.xts[miss.val[j], short.hist.var] <- X.mat%*%est.coef
  }
  
  # Save a copy of the updated dataset
  fitted.dat <- dat.xts
  
  
  # No. of bootstrap samples
  M=5000
  # Make a matrix to hold the skewness and kurtosis for each bootstrap sample.
  mom.page.mat <- matrix(0, nrow = M, ncol = 2)
  colnames(mom.page.mat) <- c("skewness", "kurtosis")
  
  # Bootstrap sampling - draw a random sample of size k (no. of missing obs) from 
  # (n-k) residuals with replacement. These k residuals are added to the k fitted 
  # values of MSCI. Then skweness and kurtosis are computed for whole dataset
  for(i in 1:M){
    boot.resid <- sample(as.numeric(err.mat), num.miss, replace=T)
    
    for (k in miss.val) {
      dat.xts[k, 2] <- dat.xts[k, 2] + as.numeric(boot.resid)[k]
    }
    
    mom.page.mat[i, 1] <- skewness(dat.xts[miss.val,2])
    mom.page.mat[i, 2] <- kurtosis(dat.xts[miss.val,2])
    
    dat.xts <- fitted.dat
  }
  
  # Mean values of skewness and kurtosis are calculated
  skew.page <- mean(mom.page.mat[, 1])
  kurt.page <- mean(mom.page.mat[, 2])
  
  page.mat <- matrix(c(skew.page, kurt.page), nrow=1, ncol=2)
  colnames(page.mat) <- c("skewness", "kurtosis")
  
  return(page.mat)
}