# UnequalReturnHist

## GSoC 2017 Unequal Returns History Project

For the GSoC 2017 Unequal Returns History project I proposed to do the following. 

### Basic Deliverables

- Create an R package "UnequalReturnHist" that implements all three of the methods "Multiple Imputation" (MI), "Combined Backfilling" (CBF), and "Factor Model Monte Carlo" (FMMC) with suitable input and output data structures to support separate or simultaneous use of the three method for: (a) optimization of portfolios whose returns have unequal histories, and (b) risk and performance analysis of portfolios of assets and of optimized portfolios.

- Procedure to implement the MI method via the function uneqhistMI() with a method for computing standard errors for the risk and performance estimators obtained with the MI backfill.

- Procedure to implement the CBF method via the function uneqhistCBF() and standard errors for CBF risk and performance measures to be computed with a bootstrap method like that used by the FMMC method.

### Modification of Deliverables

Shortly after beginning the project it was decided by mutual agreement with my mentors Doug Martin and Yindeng Jiang that there was quite enough to do without including the FMMC method, and so that was dropped.

### Goals Achieved

Successfully created the R package unequalReturnHist that includes the following key functions. 

- uneqhistMI: This function takes an xts object as input and based on function arguments provided returns the risk measures or the Global Minimum Variance (GMV) portfolio statistics such as portfolio weights, portfolio reurn, portfolio standard deviation, portfolio sharpe ratio. It also returns the standard error of Sharpe ratio for different assets in the returns dataset.  Mention the two different methods of computing S.E.'s, i.e., bootstrap and method in Little and Rubin book.  What about supporting functions?

- uneqhistCBF: This function takes an xts object as input and based on function arguments provided returns a matrix of risk measures for the whole combined backfilled dataset or a list containing the Global Minimum Variance portfolio weights and portfolio statistics such as portfolio return, portfolio standard deviation and portfolio sharpe ratio.


Work remaining to be Completed and Targeted for September-October 

- I need to add a method of computing standard errors of risk and performance measure estimators based on the combined backfilled (CBF) method.  This method will be similar to the bootstrap method used in Jiang and Martin (2015).  Paper title and journal here.

-	Add some test cases for the two main functions

-	For each of the MI method, an option will to be provided to fit one of several distributions, including normal, symmetric t, and skewed t-distributions to the regression residuals, and replace the bootstrapped residuals with draws from the fitted distribution.

-	Create a very extensive and detailed vignette for the "UnequalReturnHist" package that includes a substantial number of usage examples.
 


## Installation

To get started, you can install the package from github using devtools.

``` r
# install.packages("devtools")
devtools::install_github("spushpak/UnequalReturnHist")
```

## Example

This is a basic example which shows how to use different functions:

``` r
library(UnequalReturnHist)

## sample data
data("hfunds5_ue_ts")
head(hfunds5_ue_ts)

## Combined backfill
uneqhistCBF(hfunds5_ue_ts, FUN = "riskMeasures")
uneqhistCBF(hfunds5_ue_ts, FUN = "gmvPortfolio")

## Multiple Imputation
uneqhistMI(hfunds5_ue_ts, FUN = "riskMeasures")
uneqhistMI(hfunds5_ue_ts, FUN = "riskMeasures", M = 1000, saveReps = TRUE)
uneqhistMI(hfunds5_ue_ts, FUN = "gmvPortfolio", M = 1000, saveReps = TRUE)


## simulated data
data("sim_dat")
head(sim_dat)

## Combined backfill
uneqhistCBF(sim_dat, FUN = "riskMeasures")
uneqhistCBF(sim_dat, FUN = "gmvPortfolio")

## Multiple Imputation
uneqhistMI(sim_dat, FUN = "riskMeasures")
uneqhistMI(sim_dat, FUN = "riskMeasures", M = 1000, saveReps = TRUE)
uneqhistMI(sim_dat, FUN = "gmvPortfolio", M = 1000, saveReps = TRUE)
```
