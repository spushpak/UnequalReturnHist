constructGMVPortfolio <- function(ret.dat) {

  ones <- matrix(c(1, 1), nrow=ncol(ret.dat), ncol=1)
  
  sigma.mat <- cov(na.omit(ret.dat))
  #cor.mat <- cor(na.omit(ret.dat))
  
  numerator <- solve(sigma.mat)%*%ones
  denominator <- (t(ones)%*%solve(sigma.mat))%*%ones
  
  w_gmv <- numerator / as.numeric(denominator)
  
  return(w_gmv)
}
























