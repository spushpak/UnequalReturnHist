#' @title Implements Combined Backfill for Multiple Asset Groups
#'   
#' @description This function implements Combined Backfill for Multiple Asset
#'   Groups
#'  
#' @importFrom zoo coredata index 
#' @importFrom xts xts
#' @importFrom moments skewness kurtosis
#' @importFrom stats as.formula coef lm resid sd cov
#'   
#' @param dat_xts xts object containing the returns data for multiple assets with unequal return history.
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
#' @references
#' Jiang, Y. and Martin, R. D. (2016). "Turning Long and Short Return Histories into Equal Histories: A Better Way to Backfill Returns", https://ssrn.com/abstract=2833057.    
#'   
#' @rdname uneqhistCBF
#' @export

uneqhistCBF <- function(dat_xts, FUN){
  
  # Convert the data to xts object
  # dat_xts <- xts(dat_mat[, -1], order.by = as.Date(dat_mat[, 1], "%m/%d/%Y"))
  
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
  
  est_coef <- matrix(NA, nrow = 1+ncol(dat_xts), 
                     ncol = ncol(dat_xts) - length(long_hist_var))
  rownames(est_coef) <- c("Intercept", colnames(dat_xts))
  colnames(est_coef) <- miss_hist_var
  
  err_mat <- matrix(NA, nrow = nrow(dat_xts), ncol = ncol(dat_xts) - length(long_hist_var)) 
  colnames(err_mat) <- miss_hist_var
  
  # Start with the short history columns which have the smallest number of 
  # observations missing i.e. this group is the longest short history group. Then 
  # move on to the next shorter history group and so on.
  for (i in miss_itr) {
    dep_var <- rownames(na_count)[na_count[, "numberMissing"]==i]
    indep_var <- rownames(na_count)[na_count[, "numberMissing"] < i]
    indep_var <- c(long_hist_var, indep_var)
    
    regressor_list <- paste(indep_var, collapse = "+")
    reg_dat <- dat_xts[, c(indep_var, dep_var), drop=F]
    miss_val <- which(is.na(dat_xts[, dep_var[1]]))
    
    for (j in dep_var) {
      reg_eqn <- paste(j, "~", regressor_list)
      reg <- lm(as.formula(reg_eqn), data=reg_dat)
      betas <- as.matrix(coef(reg))
      row.names(betas)[1] <- "Intercept"
      
      est_coef[which(row.names(est_coef) %in% row.names(betas)), j] <- betas
      err_mat[(i+1):nrow(err_mat), j] <- as.matrix(resid(reg))
    }
    
    num_resid <- nrow(as.matrix(resid(reg)))
    num_miss <- i
    
    # Stack the 'new_dat' equal to number of residual times
    temp_newdat <- matrix(rep(t(coredata(new_dat)), num_resid), ncol= ncol(new_dat), 
                           byrow=TRUE)
    colnames(temp_newdat) <- colnames(new_dat)
     
    # Number of times the short_hist_var should be repeated
    rep_num <- nrow(temp_newdat)/full_length
    dep_dat <- matrix(rep(t(coredata(dat_xts[, dep_var, drop=F])), rep_num), 
                      ncol= length(dep_var), byrow=TRUE)
    colnames(dep_dat) <- dep_var
         
    # Merge the existing new_dat and temp_newdat
    new_dat <- cbind(temp_newdat, dep_dat)

    for (k in dep_var){
      X_mat <- cbind(1, new_dat[which(is.na(new_dat[, k])), indep_var, drop=F])
      betas <- est_coef[c("Intercept", indep_var), k, drop=F]
      new_dat[which(is.na(new_dat[, k])), k] <- X_mat%*%betas
    }

    # Crate the start index of each repeating block
    rep_indx <- seq(from = miss_val[1]-1, by = nrow(new_dat)/num_resid, 
                   length.out = num_resid)
     
    # Add the  respective residual to the fitted values of short history vars
    for (m in 1:length(rep_indx)) {
     for (n in 1:num_miss) {
       new_dat[rep_indx[m]+n, dep_var] <- new_dat[rep_indx[m]+n, dep_var] + err_mat[m+num_miss, dep_var]
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
