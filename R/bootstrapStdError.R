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

  # Find which columns have 'NA' values; so these columns have shorter histtory
  miss_hist_var <- colnames(dat_xts)[apply(dat_xts, 2, anyNA)]
  
  # Find which column has full history
  long_hist_var <- setdiff(colnames(dat_xts), miss_hist_var)
  
  # Find how many observatins are missing for each short history columns
  na_count <- apply(dat_xts, 2, function(x) sum(is.na(x)))
  na_count <- as.matrix(na_count)
  colnames(na_count) <- "numberMissing"
  
  # Drop full history group from 'na_count' as 'NA' count for this group is zero
  na_count <- na_count[miss_hist_var, "numberMissing", drop=F]
  
  # Sort count of missing obs from smallest to largest
  miss_itr <- sort(unique(na_count[,"numberMissing"]))
  
  # Full length of the dataset
  full_length <- nrow(dat_xts)

  # Form the combined backfilled dataset   
  new_dat <- constructCBFData(dat_xts)
  
  B <- 1000
  
  risk_metrics <- vector("list", B)
  bootstrap_list <- vector("list", B)

  for(b in 1:B) {
    boot_sample <- sample(1:nrow(new_dat), size=full_length, replace=TRUE)
    boot_dat <- coredata(new_dat)[sort(boot_sample), ]
    boot_dat <- as.data.frame(boot_dat)
    
    for(i in miss_itr){
      column_names <- row.names(na_count)[na_count[, "numberMissing"] == i]
      boot_dat[1:i, column_names] <- NA
    }
    
    # Form the combined backfilled dataset   
    cbf_dat <- constructCBFData(boot_dat)
    
    risk_metrics[[b]] <- constructRiskStats(cbf_dat)
  }  # End of bootstrap loop
  
  risk_vals <-  Reduce("+", risk_metrics) / length(risk_metrics)

  for(k in 1:B) {
    bootstrap_list[[k]] <- (risk_metrics[[k]] - risk_vals)^2
  }

  boot_std_error <- sqrt(Reduce("+", bootstrap_list) / (B - 1))
  row.names(boot_std_error) <- paste(row.names(boot_std_error), 
                                     "._stdErr_CBF_btstrp", sep = "")
  return(boot_std_error)
}

