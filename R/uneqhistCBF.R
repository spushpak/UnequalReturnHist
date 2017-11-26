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
  
  new_dat <- constructCBFData(dat_xts)

  risk_metrics <- constructRiskStats(new_dat)  
  row.names(risk_metrics) <- paste(row.names(risk_metrics), "_CBF", sep = "")
   
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
