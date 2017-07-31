constructGMVPortfolio <- function(ret_dat) {

  ones <- matrix(rep(1, ncol(ret_dat)), nrow=ncol(ret_dat), ncol=1)
  
  sigma_mat <- cov(na.omit(ret_dat))
  #cor_mat <- cor(na.omit(ret_dat))
  
  numerator <- solve(sigma_mat)%*%ones
  denominator <- (t(ones)%*%solve(sigma_mat))%*%ones
  
  w_gmv <- numerator / as.numeric(denominator)
  
  return(w_gmv)
}
























