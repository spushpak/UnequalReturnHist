#' @title Implements the regression of long history asset on short history asset.
#' 
#' @description This function implements the regression of long history asset on short history asset.
#' 
#' 
#' 
#' @param short.hist.var is the asset name with missing history.
#' @param long.hist.var is the asset name with full history.
#' @param ret.dat is the asset returns data containg the full history asset returns
#' and one or more of the short history asset returns.
#' @param num.miss is the number of missing observations.
#' 
#' @return
#' \item{reg.dat}{A list containing a fitted values, residuals etc.}
#' 
#'
#'    
#' @author Pushpak Sarkar
#' 
#' 
#' 
#' 
#' 
#' 
#'  


#### Do the regression and find fitted value for the missing returns
fitValue <- function(short.hist.var, long.hist.var, ret.dat, num.miss){
  
  regressor.list <- paste(long.hist.var, collapse = "+")
  
  est.coef <- matrix(0, nrow=length(long.hist.var)+1, ncol=length(short.hist.var),
                     dimnames=list(c("Intercept", long.hist.var), short.hist.var))
  
  num.miss <- num.miss
  len.diff <- nrow(ret.dat) - num.miss
  
  err.mat <- matrix(0, nrow=len.diff, ncol=length(short.hist.var))
  colnames(err.mat) <- short.hist.var
  
  for (j in short.hist.var) {
    reg.eqn <- paste(j, "~", regressor.list)
    reg <- lm(as.formula(reg.eqn), data=ret.dat)
    
    est.coef[ , j] <- as.matrix(coef(reg))
    err.mat[ , j] <- as.matrix(resid(reg))
  }
  
  miss.val <- which(is.na(ret.dat[, short.hist.var[1]]))
  
  for (k in 1:num.miss) {
    X.mat <- cbind(1, ret.dat[miss.val[k], long.hist.var])
    ret.dat[miss.val[k], short.hist.var] <- X.mat%*%est.coef
  }
  
  # Remove the long.hist.var
  ret.dat <- ret.dat[, ! names(ret.dat) %in% long.hist.var, drop=FALSE]
  
  # Create a list of fitted values, error matrix and missing values
  reg.dat <- list(fitted.dat=ret.dat, err.mat=err.mat, miss.val=miss.val)
  
  return(reg.dat)
}
