#' @title Constructs the new dataset by Combined Backfill Method.
#' 
#' @description Constructs the new dataset by Combined Backfill Method.
#' 
#' @importFrom zoo coredata
#' 
#' @param reg.dat is the regression data i.e. fitted data, residuals from regression.
#' @param new.dat is the backfilled dataset constructed for the previous 
#' missing history asset.
#' @param full.length is the length of the full history asset. 
#' 
#' @return
#' \item{new.dat}{is the backfilled dataset constructed by CBF.}
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
#' 


###########################################################################
### Constrcut the combined backfilled dataset based on the number of residuals

constructCBFData <- function(reg.dat, new.dat, full.length){
  
  num.resid <- nrow(reg.dat$err.mat)
  miss.val <- reg.dat$miss.val
  err.mat <- reg.dat$err.mat
  num.miss <- length(miss.val)
  short.hist.var <- colnames(reg.dat$err.mat)
  
  # Stack the original block 'new.dat' equal to number of residual times
  temp.newdat <- matrix(rep(t(coredata(new.dat)), num.resid), ncol= ncol(new.dat), 
                        byrow=TRUE)
  colnames(temp.newdat) <- colnames(new.dat)
  
  # Number of times the 'fitted.dat' i.e. short.hist.var should be repeated
  rep.num <- nrow(temp.newdat)/full.length
  temp.fitdat <- matrix(rep(t(coredata(reg.dat$fitted.dat)), rep.num), 
                        ncol= ncol(reg.dat$fitted.dat), byrow=TRUE)
  colnames(temp.fitdat) <- colnames(reg.dat$fitted.dat)
  
  # Merge the existing new.dat and fitted.dat
  new.dat <- cbind(temp.newdat, temp.fitdat)
  
  # Crate the start index of each repeating block
  rep.indx <- seq(from = miss.val[1]-1, by = nrow(new.dat)/num.resid, 
                  length.out = num.resid)
  

  # Add the  respective residual to the fitted values of short history vars
  for (m in 1:length(rep.indx)) {
    for (n in 1:num.miss) {
      new.dat[rep.indx[m]+n, short.hist.var] <- new.dat[rep.indx[m]+n, short.hist.var] + err.mat[m, short.hist.var]
    }
  }
  return(new.dat)
}
