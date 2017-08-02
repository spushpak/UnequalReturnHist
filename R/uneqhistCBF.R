#' @title Implements Combined Backfill for Multiple Asset Groups
#'   
#' @description This function implements Combined Backfill for Multiple Asset
#'   Groups
#'   
#' @importFrom xts xts
#' @importFrom moments skewness kurtosis
#' @importFrom stats as.formula coef lm resid sd
#'   
#' @param dat_mat returns data for multiple assets with unequal return history.
#' @param FUN indicates whether risk measures or the Global Minimum Variance
#'   (GMV) portfolio statistics such as portfolio weights, portfolio reurn,
#'   portfolio standard deviation, portfolio sharpe ratio need to be computed.
#'   Two possible values are "riskMeasures" and "gmvPortfolio".
#'   
#'   
#' @return Based on the value of the \code{FUN} argument the function either
#' returns the risk measures or the GMV portfolio statistics. 
#' \item{risk_metrics}{A matrix of risk measures for the whole combined
#' backfilled dataset.} \item{gmvPortfolio_list}{A list containing the Global
#' Minimum Variance portfolio weights and portfolio statistics such as portfolio
#' return, portfolio standard deviation and portfolio sharpe ratio.}
#' 
#' 
#' @author Pushpak Sarkar
#'   
#'   
#' @rdname uneqhistCBF
#' @export

uneqhistCBF <- function(dat_mat, FUN){
  
  # Convert the data to xts object (probably it's not necessary)
  dat_xts <- xts(dat_mat[, -1], order.by = as.Date(dat_mat[, 1], "%m/%d/%Y"))
  
  # Keep a copy of original data - will be updated with fitted values for the
  # missing portion
  fitted_xts <- dat_xts 
  
  # Find which columns have 'NA' values; so these columns have shorter histtory
  miss_hist_var <- colnames(dat_xts)[apply(dat_xts, 2, anyNA)]
  
  # Find which column has full history
  long_hist_var <- setdiff(colnames(dat_xts), miss_hist_var)
  
  # Find how many observatins are missing for each short history columns
  na_count <- apply(dat_xts, 2, function(x) sum(is.na(x)))
  na_count <- as.matrix(na_count)
  colnames(na_count) <- "numberMissing"
  
  # Alternative way to find the full history column
  #long_hist_var <- rownames(na_count)[na_count[, "numberMissing"]==0]
  
  # Drop full history group from 'na_count' as 'NA' count for this group is zero
  na_count <- na_count[miss_hist_var, "numberMissing", drop=F]
  
  # Sort count of missing obs from smallest to largest
  miss_itr <- sort(unique(na_count[,"numberMissing"]))
  
  # Full length of the dataset
  full_length <- nrow(dat_xts)
  
  new_dat <- dat_xts[, long_hist_var]
  
  # Start with the short history columns which have the smallest number of 
  # observations missing i.e. this group is the longest short history group. Then 
  # move on to the next shorter history group and so on.
  for (i in miss_itr) {
    dep_var <- rownames(na_count)[na_count[, "numberMissing"]==i]
    indep_var <- rownames(na_count)[na_count[, "numberMissing"] < i]
    indep_var <- c(long_hist_var, indep_var)
    temp_dat <- fitted_xts[, c(indep_var, dep_var), drop=F]
    num_miss <- i
    reg_dat <- fitValue(dep_var, indep_var, temp_dat, num_miss)
    fitted_xts[reg_dat$miss_val, dep_var] <- reg_dat$fitted_dat[reg_dat$miss_val, dep_var]
    
    num_resid <- nrow(reg_dat$err_mat)
    miss_val <- reg_dat$miss_val
    err_mat <- reg_dat$err_mat
    num_miss <- length(miss_val)
    
    # Stack the original block 'new_dat' equal to number of residual times
    temp_newdat <- matrix(rep(t(coredata(new_dat)), num_resid), ncol= ncol(new_dat), 
                          byrow=TRUE)
    colnames(temp_newdat) <- colnames(new_dat)
    
    # Number of times the 'fitted_dat' i.e. short_hist_var should be repeated
    rep_num <- nrow(temp_newdat)/full_length
    temp_fitdat <- matrix(rep(t(coredata(reg_dat$fitted_dat)), rep_num), 
                          ncol= ncol(reg_dat$fitted_dat), byrow=TRUE)
    colnames(temp_fitdat) <- colnames(reg_dat$fitted_dat)
    
    # Merge the existing new_dat and fitted_dat
    new_dat <- cbind(temp_newdat, temp_fitdat)
    
    # Crate the start index of each repeating block
    rep_indx <- seq(from = miss_val[1]-1, by = nrow(new_dat)/num_resid, 
                    length.out = num_resid)
    
    # Add the  respective residual to the fitted values of short history vars
    for (m in 1:length(rep_indx)) {
      for (n in 1:num_miss) {
        new_dat[rep_indx[m]+n, dep_var] <- new_dat[rep_indx[m]+n, dep_var] + err_mat[m, dep_var]
      }
    }
  } 
  
######################################################################
  risk_metrics <- constructRiskStats(new_dat)  
   
  if (FUN == "gmvPortfolio") {
    w_gmv <- constructGMVPortfolio(new_dat)
    colnames(w_gmv) <- "portfolio_wts"
    
    portfolio_ret <-  as.numeric(risk_metrics["Mean", ]%*%w_gmv)
    portfolio_sd <- sqrt(as.numeric(t(w_gmv)%*%cov(new_dat)%*%w_gmv))
    portfolio_SR <- portfolio_ret / portfolio_sd
    
    portfolio_stats <- cbind(portfolio_ret, portfolio_sd, portfolio_SR)
    colnames(portfolio_stats) <- c("Return", "Std Dev", "Sharpe Ratio")
    row.names(portfolio_stats) <- "CBF"
    
    gmvPortfolio_list <- list("Portfolio_Weights" = w_gmv, 
                              "Portfolio_Stats" = portfolio_stats)
  }  
  
  if(FUN == "gmvPortfolio") 
    return(gmvPortfolio_list) 
  else if(FUN == "riskMeasures")
    return(risk_metrics)
  else
    return("Not a valid function name.")
}
