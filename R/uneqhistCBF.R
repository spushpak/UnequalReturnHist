#' @title Implements Combined Backfill for Multiple Asset Groups
#' 
#' @description This function implements Combined Backfill for Multiple Asset Groups
#' 
#' @importFrom xts xts
#' @importFrom moments skewness kurtosis
#' @importFrom stats as.formula coef lm resid sd
#' 
#' @param dat.mat is the returns data for multiple assets with unequal return history.
#' @param saveReps is a logical flag, to indicate whether the risk measures  
#' need to be computed for all the replicates. Default is FALSE.
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
#' @rdname uneqhistCBF
#' @export

uneqhistCBF <- function(dat.mat, saveReps=FALSE){
  
  # Convert the data to xts object (probably it's not necessary)
  dat.xts <- xts(dat.mat[, -1], order.by = as.Date(dat.mat[, 1], "%m/%d/%Y"))
  
  # Keep a copy of original data
  old.xts <- dat.xts 
  
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
  
  new.dat <- dat.xts[, long.hist.var]
  
  # Start with the short history columns which have the smallest number of 
  # observations missing i.e. this group is the longest short history group. Then 
  # move on to the next shorter history group and so on.
  for (i in miss.itr) {
    short.hist.var <- rownames(na_count)[na_count[, "num.miss"]==i]
    temp.dat <- dat.xts[, c(long.hist.var, short.hist.var), drop=F]
    num.miss <- i
    reg.dat <- fitValue(short.hist.var, long.hist.var, temp.dat, num.miss)
    new.dat <- constructCBFData(reg.dat, new.dat, full.length)  
  } 
  
  
  ######################################################################
  # Create the start index of each repeating block in the new dataset. Each repeating
  # block has length equal to the row numbers of initial dataset.
  row.grp <- seq(1, nrow(new.dat), by=full.length)
  
  # No. of repeatng blocks
  num.block <- nrow(new.dat)/full.length
  
  # Create an empty list to contain risk measures of each repeating block
  risk.metrics <- vector("list", num.block)
  
  # Create an empty matrix to hold the aggregate risk measures of the CBF data 
  risk.vals <- matrix(0, nrow=6, ncol=length(miss.hist.var),
                      dimnames=list(c("Skewness", "Kurtosis", "Mean", "Volatility", "Sharpe Ratio", "Expected Shortfall"), miss.hist.var))
  
  
  # Repeat the above step for each repeating block
  for (j in 0:(num.block-1)) {
    risk.mat <- matrix(0, nrow=6, ncol=length(miss.hist.var),
                       dimnames=list(c("Skewness", "Kurtosis", "Mean", "Volatility", "Sharpe Ratio", "Expected Shortfall"), miss.hist.var))
    
    start.indx <- j*full.length+1  
    end.indx <- (j+1)*full.length
    
    # Extract the block based on start and end index    
    block.dat <- new.dat[start.indx:end.indx, miss.hist.var, drop=F]
    
    # Compute risk measures for the extracted block
    risk.mat["Skewness", ] <- moments::skewness(block.dat)
    risk.mat["Kurtosis", ] <- moments::kurtosis(block.dat)
    risk.mat["Mean", ] <- colMeans(block.dat)
    risk.mat["Volatility", ] <- apply(block.dat, 2, sd)
    risk.mat["Sharpe Ratio", ] <- apply(block.dat, 2, SharpeRatio, FUN="StdDev")
    risk.mat["Expected Shortfall", ] <- apply(block.dat, 2, ES, p=0.95, method="historical")
    
    risk.metrics[[j+1]] <- risk.mat
  }
  
  
  # Constrcut the risk matrix for the aggregate CBF data
  
  risk.vals["Skewness", ] <- moments::skewness(new.dat[, miss.hist.var])
  risk.vals["Kurtosis", ] <- moments::kurtosis(new.dat[, miss.hist.var])
  risk.vals["Mean", ] <- colMeans(new.dat[, miss.hist.var])
  risk.vals["Volatility", ] <- apply(new.dat[, miss.hist.var], 2, sd)
  risk.vals["Sharpe Ratio", ] <- apply(new.dat[, miss.hist.var], 2, SharpeRatio, 
                                       FUN="StdDev")
  risk.vals["Expected Shortfall", ] <- apply(new.dat[, miss.hist.var], 2, ES, 
                                             p=0.95, method="historical")
  
  
  if(saveReps==TRUE) 
    return(round(risk.metrics, digits = 3)) 
  else 
    return(round(risk.vals, digits = 3))
}


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


### Constrcut the combined backfilled dataset based on the number of residuals
constructCBFData <- function(reg.dat, new.dat, full.length){
  
  num.resid <- nrow(reg.dat$err.mat)
  miss.val <- reg.dat$miss.val
  err.mat <- reg.dat$err.mat
  num.miss <- length(miss.val)
  short.hist.var <- colnames(reg.dat$err.mat)
  
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

