#' @title Implements Multiple Imputation Method folowing Page(2013)
#' 
#' @description This function implements Multiple Imputation Method folowing 
#' Page(2013)
#' 
#' @importFrom xts xts
#' @importFrom moments skewness kurtosis
#' 
#' 
#' @param dat_xts xts object containing the returns data for multiple assets with unequal return history.
#' @param FUN indicates whether risk measures or the Global Minimum Variance (GMV) 
#' portfolio statistics such as portfolio weights, portfolio reurn, portfolio standard
#' deviation, portfolio sharpe ratio need to be computed. Two possible values
#' are "riskMeasures" and "gmvPortfolio".
#' @param M is the number of replicates for multiple imputation. Default is 100.
#' @param saveReps is a logical flag, to indicate whether the risk measures  
#' and GMV portfolio statistics need to be computed for all the replicates. Default is FALSE.
#' 
#' 
#' 
#' @return
#' Based on the value of the \code{FUN} argument the function either returns the risk 
#' measures or the GMV portfolio statistics. If \code{saveReps} = TRUE, risk measures or 
#' GMV portfolio statistics are returned for each replicate else the average over all
#' the replicates are returned. 
#' \item{risk_metrics}{If \code{FUN} = "riskMeasures" and \code{saveReps} = TRUE, a list 
#'  of length \code{M} containing risk measures for each replicate are returned.}
#' \item{risk_vals}{If \code{FUN} = "riskMeasures" and \code{saveReps} = FALSE, a matrix 
#'  containing average of the risk measures over all replicates are returned.}
#' \item{gmvPortfolio_list}{If \code{FUN} = "gmvPortfolio" and \code{saveReps} = TRUE, a list 
#'  of length \code{M} containing GMV portfolio statistics for each replicate are retruned.}
#' \item{gmvPortfolio}{If \code{FUN} = "gmvPortfolio" and \code{saveReps} = FALSE, a list 
#'  containing average of the GMV portfolio statistics over all replicates are returned.}
#'    
#' @author Pushpak Sarkar
#' 
#' 
#' @rdname uneqhistMI
#' @export

uneqhistMI <- function(dat_mat, FUN, M=100, saveReps = FALSE){
  
  # Convert the data to xts object
  #dat_xts <- xts(dat_mat[, -1], order.by = as.Date(dat_mat[, 1], "%m/%d/%Y"))
  
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
    #temp_dat <- fitted_xts[, c(indep_var, dep_var), drop=F]
    
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
  }   ############# end of for loop (for each short history group)
  
  # No. of bootstrap samples
  M <- M
  
  # Create an empty list to hold the risk and performance measures
  risk_metrics <- vector("list", M)
  
  # Create a list to hold portfolio weights for all replicates
  portfolio_wtlist <- vector("list", M)
  
  # Create a list to hold the covariance of the assets for all replicates 
  cov_gmvport <- vector("list", M)
  
  # Create a list to hold portfolio return, sd, sharpe ratio for replicates
  portfolio_stats <- vector("list", M)
  gmvPortfolio_list <- vector("list", M)
  
  # Bootstrap sampling - draw a random sample of size k (no. of missing obs) from 
  # (n-k) residuals with replacement. These k residuals are added to the k fitted 
  # values.
  for(i in 1:M){
    new_dat <- dat_xts
    for (j in miss_itr) {
      dep_var <- rownames(na_count)[na_count[, "numberMissing"]==j]
      indep_var <- rownames(na_count)[na_count[, "numberMissing"] < j]
      indep_var <- c(long_hist_var, indep_var)
      
      miss_val <- which(is.na(err_mat[, dep_var[1]]))
      num_miss <- j
      sample_indx <- which(!is.na(err_mat[, dep_var[1]]))
      
      r_sample <- sample(sample_indx, size=num_miss, replace=TRUE)
      boot_resid <- err_mat[r_sample, dep_var, drop=F]
      colnames(boot_resid) <- dep_var
      boot_resid <- as.xts(boot_resid, order.by = index(dat_xts)[miss_val])
      
      # Compute the fitted values by using estimated beta coefficients
      for (k in dep_var){
        X_mat <- cbind(1, new_dat[which(is.na(new_dat[, k])), indep_var, drop=F])
        betas <- est_coef[c("Intercept", indep_var), k, drop=F]
        new_dat[which(is.na(dat_xts[, k])), k] <- X_mat%*%betas
        new_dat[miss_val, k] <- new_dat[miss_val, k, drop=F] + boot_resid[, k, drop=F]  
      }
    }  
    
    risk_metrics[[i]] <- constructRiskStats(new_dat)
    
    w_gmv <- constructGMVPortfolio(new_dat)
    colnames(w_gmv) <- "portfolio_wts"
    portfolio_wtlist[[i]] <- t(w_gmv) 
    
    portfolio_ret <-  as.numeric(risk_metrics[[i]]["Mean", ]%*%w_gmv)
    portfolio_sd <- sqrt(as.numeric(t(w_gmv)%*%cov(new_dat)%*%w_gmv))
    portfolio_SR <- portfolio_ret / portfolio_sd
    port_stats <-  cbind(portfolio_ret, portfolio_sd, portfolio_SR)
    colnames(port_stats) <- c("Return", "Std Dev", "Sharpe Ratio")
    row.names(port_stats) <- "MI"
    portfolio_stats[[i]] <- port_stats
    
    gmvPortfolio_list[[i]] <- list("Portfolio_Weights" = w_gmv, 
                              "Portfolio_Stats" = port_stats)
    
    cov_gmvport[[i]] <- cov(new_dat)
  }
  
  risk_vals <-  Reduce("+", risk_metrics) / length(risk_metrics)
  
  # Compute the standard error
  bootstrap_list <- vector("list", M)
  SR_list <- vector("list", M)
  
  for(k in 1:M) {
    bootstrap_list[[k]] <- (risk_metrics[[k]] - risk_vals)^2
    SR <- risk_metrics[[k]]["Sharpe Ratio", ]
    sk <- risk_metrics[[k]]["Skewness", ]
    kurt <- risk_metrics[[k]]["Kurtosis", ]
    SR_list[[k]] <- sqrt(1 - sk*SR + ((kurt+2)/4)*SR^2)/sqrt(full_length)
  }
  
  boot_std_error <- sqrt(Reduce("+", bootstrap_list) / (M - 1))
  row.names(boot_std_error) <- paste(row.names(boot_std_error), 
                                     ".bootstrap_stdErr", sep = "")
  Sharpe_Ratio_stdErr <- Reduce("+", SR_list)/M
  
  risk_vals <- rbind(risk_vals, boot_std_error["Sharpe Ratio.bootstrap_stdErr", , drop=F], 
                     Sharpe_Ratio_stdErr)
  
  risk_vals <- risk_vals[, order(colnames(risk_vals))]
  risk_vals <- round(risk_vals, digits = 4) 
  
  if (FUN == "gmvPortfolio") { 
    # Portfolio weights and stats are found for each replicate (Method 1)
    portfolio_stats_avg <- Reduce("+", portfolio_stats ) / length(portfolio_stats)
    w_gmv1 <- Reduce("+", portfolio_wtlist ) / length(portfolio_wtlist)
  
    portfolio_ret1 <- as.numeric(portfolio_stats_avg[, "Return"])
    portfolio_sd1 <- as.numeric(portfolio_stats_avg[, "Std Dev"])
    portfolio_SR1 <- as.numeric(portfolio_stats_avg[, "Sharpe Ratio"])

    # Method 2 - taking the average of the covariance matrices
    avg_cov <- Reduce("+", cov_gmvport) / length(cov_gmvport)
    inv_avg_cov <- solve(avg_cov)
  
    ones <- matrix(rep(1, ncol(dat_xts)), nrow=ncol(dat_xts), ncol=1)
    numerator <- inv_avg_cov%*%ones
    denominator <- t(ones)%*%inv_avg_cov%*%ones
    w_gmv2 <- numerator / as.numeric(denominator)
  
    # Using the average of the mean return vector and the portfolio vector
    portfolio_ret2 <- as.numeric(risk_vals["Mean", , drop=F]%*%w_gmv2)
    portfolio_sd2 <- sqrt(as.numeric(t(w_gmv2)%*%avg_cov%*%w_gmv2))
    portfolio_SR2 <- portfolio_ret2 / portfolio_sd2

    # Method 3 - invert each cov matrix in the list
    inv_cov_gmvport <- lapply(cov_gmvport, solve)
    avg_inv_cov <- Reduce("+", inv_cov_gmvport) / length(inv_cov_gmvport)
  
    ones <- matrix(rep(1, ncol(dat_xts)), nrow=ncol(dat_xts), ncol=1)
    numerator <- avg_inv_cov%*%ones
    denominator <- t(ones)%*%avg_inv_cov%*%ones
    w_gmv3 <- numerator / as.numeric(denominator)

    # Using the average of the mean return vector and the portfolio vector
    portfolio_ret3 <- as.numeric(risk_vals["Mean", , drop=F]%*%w_gmv3)
    portfolio_sd3 <- sqrt(as.numeric(t(w_gmv3)%*%solve(avg_inv_cov)%*%w_gmv3))
    portfolio_SR3 <- portfolio_ret3 / portfolio_sd3

    port_wts <- rbind(w_gmv1, t(w_gmv2), t(w_gmv3))
    row.names(port_wts) <- c("MI 1", "MI 2", "MI 3")
  
    port_returns <- rbind(portfolio_ret1, portfolio_ret2, portfolio_ret3)
    row.names(port_returns) <- c("MI 1", "MI 2", "MI 3")
  
    port_sds <- rbind(portfolio_sd1, portfolio_sd2, portfolio_sd3)
    row.names(port_returns) <- c("MI 1", "MI 2", "MI 3")
  
    port_SRs <- rbind(portfolio_SR1, portfolio_SR2, portfolio_SR3)
    row.names(port_returns) <- c("MI 1", "MI 2", "MI 3")
  
    comparison_mat <- cbind(port_returns, port_sds, port_SRs)
    colnames(comparison_mat) <- c("Return", "Std Dev", "Sharpe Ratio")
  
    gmvPortfolio <- list("Portfolio_Weights" = port_wts, 
                            "Portfolio_Stats" = comparison_mat)
  }

  if(FUN == "gmvPortfolio") {
    if(saveReps == TRUE) 
      return(gmvPortfolio_list) 
    else
      return(gmvPortfolio)
  }
  else if(FUN == "riskMeasures") {
    if(saveReps == TRUE) 
      return(risk_metrics)
    else 
      return(risk_vals)
  } else
      return("Not a valid function name.")
}
