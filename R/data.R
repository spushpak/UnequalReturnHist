#' @name hfunds5_ue_ts
#' @title Monthly retrun data of five different assets.
#' @description A dataset containing the monthly retrun data of five different assets.
#'
#' @format An xts object with 123 rows and 6 columns:
#' \describe{
#'   \item{E.M.}{The asset with full history from Jan 1993 to Mar 2003}
#'   \item{EQUITY}{This asset has data from Jan 1994 to Mar 2003}
#'   \item{EVENT.1}{This asset has data from Jan 1994 to Mar 2003}
#'   \item{EVENT.2}{This asset has data from Jan 1994 to Mar 2003}
#'   \item{HEALTH}{This asset has data from Jan 1997 to Mar 2003}
#'   \item{HIGH.YIELD}{This asset has data from Jan 1997 to Mar 2003}
#' }
#' 
#' @docType data
#' @keywords data
#' @source TBA
#' @usage data("hfunds5_ue_ts")
NULL

#' @name sim_dat
#' @title Simulated monthly retrun data of four different assets.
#' @description A dataset containing the simulated monthly retrun data generated from multivariate normal distribution.
#' These data have constant correlation of rho = 0.2 and mean return mu = c(0.01, 0.007, 0.008, 0.005) and volatility 
#' vol = c(0.009, 0.017, 0.020, 0.014).
#' @format An xts object with 120 rows and 4 columns:
#' \describe{
#'   \item{V1}{The asset has full history}
#'   \item{V2}{This asset has first 12 values missing}
#'   \item{V3}{This asset has first 12 values missing}
#'   \item{V4}{This asset has first 24 values missing}
#' }
#' 
#' @docType data
#' @keywords data
#' @usage data("sim_dat")
NULL

