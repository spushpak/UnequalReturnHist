#' @title Constructs Global Minimum Variance (GMV) Portfolio.
#'   
#' @description This function Constructs Global Minimum Variance (GMV) Portfolio.
#'
#' @importFrom stats na.omit
#'   
#' @param ret_dat xts object containing the returns data for multiple assets.
#' 
#' @return The function returns the GMV portfolio weights. 
#' \item{w_gmv}{A vector of weights of different assets in the Global
#' Minimum Variance portfolio.}
#' 
#' @author Pushpak Sarkar
#'

constructGMVPortfolio <- function(ret_dat) {

  ones <- matrix(rep(1, ncol(ret_dat)), nrow=ncol(ret_dat), ncol=1)
  
  sigma_mat <- cov(na.omit(ret_dat))
  numerator <- solve(sigma_mat)%*%ones
  denominator <- (t(ones)%*%solve(sigma_mat))%*%ones
  
  w_gmv <- numerator / as.numeric(denominator)
  
  return(w_gmv)
}
























