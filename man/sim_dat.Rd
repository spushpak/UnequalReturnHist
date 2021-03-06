% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{sim_dat}
\alias{sim_dat}
\title{Simulated monthly retrun data of four different assets.}
\format{An xts object with 120 rows and 4 columns:
\describe{
  \item{V1}{The asset has full history}
  \item{V2}{This asset has first 12 values missing}
  \item{V3}{This asset has first 12 values missing}
  \item{V4}{This asset has first 24 values missing}
}}
\usage{
data("sim_dat")
}
\description{
A dataset containing the simulated monthly retrun data generated from multivariate normal distribution.
These data have constant correlation of rho = 0.2 and mean return mu = c(0.01, 0.007, 0.008, 0.005) and volatility 
vol = c(0.009, 0.017, 0.020, 0.014).
}
\keyword{data}
