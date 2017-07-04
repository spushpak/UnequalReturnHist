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
  #old.xts <- dat.xts 
  
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
    risk.mat["Sharpe Ratio", ] <- risk.mat["Mean", ] / risk.mat["Volatility", ]
    #risk.mat["Sharpe Ratio", ] <- apply(block.dat, 2, SharpeRatio, FUN = "StdDev")
    risk.mat["Expected Shortfall", ] <- apply(block.dat, 2, expectedShortfall)
    
    risk.metrics[[j+1]] <- risk.mat
  }
  
  
  # Constrcut the risk matrix for the aggregate CBF data
  risk.vals["Skewness", ] <- colMeans(do.call("rbind", lapply(risk.metrics, "[", "Skewness", TRUE)))
  risk.vals["Kurtosis", ] <- colMeans(do.call("rbind", lapply(risk.metrics, "[", "Kurtosis", TRUE)))
  risk.vals["Mean", ] <- colMeans(do.call("rbind", lapply(risk.metrics, "[", "Mean", TRUE)))
  risk.vals["Volatility", ] <- colMeans(do.call("rbind", lapply(risk.metrics, "[", "Volatility", TRUE)))
  risk.vals["Sharpe Ratio", ] <- colMeans(do.call("rbind", lapply(risk.metrics, "[", "Sharpe Ratio", TRUE)))
  risk.vals["Expected Shortfall", ] <- colMeans(do.call("rbind", lapply(risk.metrics, "[", "Expected Shortfall", TRUE)))
  
  
  
  if(saveReps==TRUE) 
    return(round(risk.metrics, digits = 4)) 
  else 
    return(round(risk.vals, digits = 4))
}
