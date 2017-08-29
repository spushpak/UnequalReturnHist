#' @title Computes the standard error of estimated risk measures by the bootstrap method.
#'   
#' @description This function computes the standard error of estimated risk measures by the bootstrap method.
#'   
#' @param dat_xts xts object containing the returns data for multiple assets with unequal returns history.
#' 
#' @return The function returns the standard error of estimated risk measures by the bootstrap method. 
#' \item{boot_std_error}{A matrix containing standard error of estimated risk measures for different assets.}
#' 
#' @author Pushpak Sarkar
#'

bootstrapStdError <- function(dat_xts) {

  B <- 1000

  risk_metrics <- vector("list", B)
  bootstrap_list <- vector("list", B)

  for(b in 1:B) {

    boot_sample <- sample(1:nrow(dat_xts), size=nrow(dat_xts), replace=TRUE)
    boot_dat <- coredata(dat_xts)[sort(boot_sample), ]
  
    # Find which columns have 'NA' values; so these columns have shorter histtory
    miss_hist_var <- colnames(boot_dat)[apply(boot_dat, 2, anyNA)]
  
    # Find which column has full history
    long_hist_var <- setdiff(colnames(boot_dat), miss_hist_var)
  
    # Find how many observatins are missing for each short history columns
    na_count <- apply(boot_dat, 2, function(x) sum(is.na(x)))
    na_count <- as.matrix(na_count)
    colnames(na_count) <- "numberMissing"
  
    # Drop full history group from 'na_count' as 'NA' count for this group is zero
    na_count <- na_count[miss_hist_var, "numberMissing", drop=F]
  
    # Sort count of missing obs from smallest to largest
    miss_itr <- sort(unique(na_count[,"numberMissing"]))
  
    # Full length of the dataset
    full_length <- nrow(boot_dat)
  
    est_coef <- matrix(NA, nrow = 1+ncol(boot_dat), 
                      ncol = ncol(boot_dat) - length(long_hist_var))
    rownames(est_coef) <- c("Intercept", colnames(boot_dat))
    colnames(est_coef) <- miss_hist_var
    
    err_mat <- matrix(NA, nrow = nrow(boot_dat), ncol = ncol(boot_dat) - length(long_hist_var)) 
    colnames(err_mat) <- miss_hist_var
    
    # Start with the short history columns which have the smallest number of 
    # observations missing i.e. this group is the longest short history group. Then 
    # move on to the next shorter history group and so on.
    for (i in miss_itr) {
      dep_var <- rownames(na_count)[na_count[, "numberMissing"]==i]
      indep_var <- rownames(na_count)[na_count[, "numberMissing"] < i]
      indep_var <- c(long_hist_var, indep_var)
      
      regressor_list <- paste(indep_var, collapse = "+")
      reg_dat <- boot_dat[, c(indep_var, dep_var), drop=F]
      miss_val <- which(is.na(boot_dat[, dep_var[1]]))
      
      for (j in dep_var) {
        reg_eqn <- paste(j, "~", regressor_list)
        reg <- lm(as.formula(reg_eqn), data=as.data.frame(reg_dat))
        betas <- as.matrix(coef(reg))
        row.names(betas)[1] <- "Intercept"
        
        est_coef[which(row.names(est_coef) %in% row.names(betas)), j] <- betas
        err_mat[(i+1):nrow(err_mat), j] <- as.matrix(resid(reg))
      }
    }   ############# end of for loop (for each short history group)
    
      new_dat <- boot_dat
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
        
        # Compute the fitted values by using estimated beta coefficients
        for (k in dep_var){
          X_mat <- cbind(1, new_dat[which(is.na(new_dat[, k])), indep_var, drop=F])
          betas <- est_coef[c("Intercept", indep_var), k, drop=F]
          new_dat[which(is.na(boot_dat[, k])), k] <- X_mat%*%betas
          new_dat[miss_val, k] <- new_dat[miss_val, k, drop=F] + boot_resid[, k, drop=F]  
        }
      }
      risk_metrics[[b]] <- constructRiskStats(new_dat)
    }  # End of bootstrap loop
  
  risk_vals <-  Reduce("+", risk_metrics) / length(risk_metrics)

  for(k in 1:B) {
    bootstrap_list[[k]] <- (risk_metrics[[k]] - risk_vals)^2
  }

  boot_std_error <- sqrt(Reduce("+", bootstrap_list) / (B - 1))
  row.names(boot_std_error) <- paste(row.names(boot_std_error), 
                                     ".bootstrap_stdErr", sep = "")
  return(boot_std_error)
}

