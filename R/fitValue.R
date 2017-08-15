#### Do the regression and find fitted value for the missing returns
fitValue <- function(dep_var, indep_var, ret_dat, num_miss){
  
  regressor_list <- paste(indep_var, collapse = "+")
  
  est_coef <- matrix(0, nrow=length(indep_var)+1, ncol=length(dep_var),
                     dimnames=list(c("Intercept", indep_var), dep_var))
  
  num_miss <- num_miss
  len_diff <- nrow(ret_dat) - num_miss
  
  err_mat <- matrix(0, nrow=len_diff, ncol=length(dep_var))
  colnames(err_mat) <- dep_var
  
  for (j in dep_var) {
    reg_eqn <- paste(j, "~", regressor_list)
    reg <- lm(as.formula(reg_eqn), data=ret_dat)
    
    est_coef[, j] <- as.matrix(coef(reg))
    err_mat[, j] <- as.matrix(resid(reg))
  }
  
  miss_val <- which(is.na(ret_dat[, dep_var[1]]))
  
  for (k in 1:num_miss) {
    X_mat <- cbind(1, ret_dat[miss_val[k], indep_var])
    ret_dat[miss_val[k], dep_var] <- X_mat%*%est_coef
  }
  
  # Remove the indep_var
  ret_dat <- ret_dat[, ! names(ret_dat) %in% indep_var, drop=FALSE]
  
  # Create a list of fitted values, error matrix and missing values
  reg_output <- list(fitted_dat=ret_dat, err_mat=err_mat, miss_val=miss_val)
  
  return(reg_output)
}

