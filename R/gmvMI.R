rm(list = ls())

library(xts)
library(MASS)

source("C:/UW/RDev/UnequalReturnHist/R/fitValue.R")
source("C:/UW/RDev/UnequalReturnHist/R/constructRiskStats.R")
source("C:/UW/RDev/UnequalReturnHist/R/constructGMVPortfolio.R")
source("C:/UW/RDev/UnequalReturnHist/R/expectedShortfall.R")


# Read all other functions in the UnequalReturnHist package into R's 
# environment

ccorMat = function(rho,p)
{
  M1 = diag(,nrow = p)
  M2 = rep(1,p)%*%t(rep(1,p))
  (1-rho)*M1 + rho*M2
}

ccorMatInv = function(rho,p)
{
  M1 = diag(,nrow = p)
  M2 = rep(1,p)%*%t(rep(1,p))
  a = rho/(1+(p-1)*rho)
  (M1-a*M2)/(1-rho) 
}

cor2cov = function(R,vol)
{
  diag(vol)%*%R%*%diag(vol)
}



rho <- 0.4
mu <- c(0.01, 0.007, 0.008, 0.015)
vol <- c(0.057, 0.086, 0.072, 0.043)

cor.mat <- ccorMat(rho, 4)
cov.mat <- cor2cov(cor.mat, vol)

set.seed(1)
ret.dat <- mvrnorm(n=60, mu, cov.mat)
colnames(ret.dat) <- c("V1","V2", "V3", "V4")

## date sequence by month
date <- as.character(seq(as.Date("2000/1/1"), by = "month", length.out = 60))
dat.xts <- xts(ret.dat, order.by = as.Date(date))
full.data <- dat.xts

v2.v3.missing <- dat.xts['2000-01-01::2000-12-01', c("V2", "V3")]
v4.missing <- dat.xts['2000-01-01::2001-12-01', "V4"]

dat.xts['2000-01-01::2000-12-01', c("V2", "V3")] <- NA
dat.xts['2000-01-01::2001-12-01', "V4"] <- NA

# Keep a copy of original data - will be updated with fitted values for the
# missing portion
fitted.xts <- dat.xts 

# Find which columns have 'NA' values; so these columns have shorter histtory
miss.hist.var <- colnames(dat.xts)[apply(dat.xts, 2, anyNA)]

# Find which column has full history
long.hist.var <- setdiff(colnames(dat.xts), miss.hist.var)

# Find how many observatins are missing for each short history columns
na_count <- apply(dat.xts, 2, function(x) sum(is.na(x)))
na_count <- as.matrix(na_count)
colnames(na_count) <- "num.miss"

# Alternative way to find the full history column
#long.hist.var <- rownames(na_count)[na_count[, "num.miss"]==0]

# Drop full history group from 'na_count' as 'NA' count for this group is zero
na_count <- na_count[miss.hist.var, "num.miss", drop=F]

# Sort count of missing obs from smallest to largest
miss.itr <- sort(unique(na_count[,"num.miss"]))

# Full length of the dataset
full.length <- nrow(dat.xts)

resid.mat <- matrix(NA, nrow = full.length, ncol = length(miss.hist.var))
colnames(resid.mat) <- miss.hist.var

# Start with the short history columns which have the smallest number of 
# observations missing i.e. this group is the longest short history group. Then 
# move on to the next shorter history group and so on.
for (i in miss.itr) {
  short.hist.var <- rownames(na_count)[na_count[, "num.miss"]==i]
  temp.dat <- dat.xts[, c(long.hist.var, short.hist.var), drop=F]
  num.miss <- i
  reg.dat <- fitValue(short.hist.var, long.hist.var, temp.dat, num.miss)
  
  fitted.xts[reg.dat$miss.val, short.hist.var] <- reg.dat$fitted.dat[reg.dat$miss.val, short.hist.var]
  resid.mat[(length(reg.dat$miss.val)+1):nrow(resid.mat), short.hist.var] <- reg.dat$err.mat[, short.hist.var]
} ############# end of for loop (for each short history group)

# No. of bootstrap samples
M <- 1000

new.dat <- fitted.xts

# Create a list to hold the risk and performance measures
risk_GMVportfolio <- vector("list", M)

# Create a list to hold the covariance of the assets for all replicates 
cov_gmvport <- vector("list", M)

portfolio.stats <- vector("list", M)

# Bootstrap sampling - draw a random sample of size k (no. of missing obs) from 
# (n-k) residuals with replacement. These k residuals are added to the k fitted 
# values.
for(i in 1:M){
  for (j in miss.itr) {
    short.hist.var <- rownames(na_count)[na_count[, "num.miss"]==j]
    miss.val <- which(is.na(resid.mat[, short.hist.var[1]]))
    num.miss <- length(miss.val)
    sample.indx <- which(!is.na(resid.mat[, short.hist.var[1]]))
    
    r.sample <- sample(sample.indx, size=num.miss, replace=TRUE)
    boot.resid <- resid.mat[r.sample, short.hist.var, drop=F]
    colnames(boot.resid) <- short.hist.var
    boot.resid <- as.xts(boot.resid, order.by = index(fitted.xts)[miss.val])
    new.dat[miss.val, short.hist.var] <- fitted.xts[miss.val, short.hist.var, drop=F] + boot.resid[, short.hist.var, drop=F]
  }
  
  #block.dat <- new.dat[, miss.hist.var, drop=F]
  #risk.metrics[[i]] <- constructRiskStats(block.dat, miss.hist.var)
  
  risk.metrics <- constructRiskStats(new.dat)
  gmvport.wt <- constructGMVPortfolio(new.dat)
  colnames(gmvport.wt) <- "GMV_Portfolio_Wt"
  
  risk_GMVportfolio[[i]] <- rbind(risk.metrics, t(gmvport.wt))
  cov_gmvport[[i]] <- cov(new.dat)
  port.ret <-  as.numeric(risk.metrics["Mean", ]%*%gmvport.wt)
  port.sd <- sqrt(as.numeric(t(gmvport.wt)%*%cov(new.dat)%*%gmvport.wt))
  port.sharpeRatio <- port.ret / port.sd
  portfolio.stats[[i]] <-  rbind(port.ret, port.sd, port.sharpeRatio)
}

# # Create a matrix to hold the aggregate risk and performance measures
# risk.vals <- matrix(0, nrow=6, ncol=length(miss.hist.var),
#                     dimnames=list(c("Skewness", "Kurtosis", "Mean", "Volatility", "Sharpe Ratio", "Expected Shortfall"), miss.hist.var))
# 
# risk.vals["Skewness", ] <- colMeans(do.call("rbind", lapply(risk.metrics, "[", "Skewness", TRUE)))
# risk.vals["Kurtosis", ] <- colMeans(do.call("rbind", lapply(risk.metrics, "[", "Kurtosis", TRUE)))
# risk.vals["Mean", ] <- colMeans(do.call("rbind", lapply(risk.metrics, "[", "Mean", TRUE)))
# risk.vals["Volatility", ] <- colMeans(do.call("rbind", lapply(risk.metrics, "[", "Volatility", TRUE)))
# risk.vals["Sharpe Ratio", ] <- colMeans(do.call("rbind", lapply(risk.metrics, "[", "Sharpe Ratio", TRUE)))
# risk.vals["Expected Shortfall", ] <- colMeans(do.call("rbind", lapply(risk.metrics, "[", "Expected Shortfall", TRUE)))
# 
# 
# # Create an empty matrix to hold the GMV portfolio weigths
# # aggregating over all the replicates
# portfolio.vals <- matrix(0, nrow=ncol(new.dat), ncol=1,
#                          dimnames=list(colnames(new.dat), "Portfolio weights"))
# 
# portfolio.vals <- colMeans(do.call("rbind", lapply(gmvport.wtlist, "[", TRUE, TRUE)))


risk.GMVport <-  Reduce("+", risk_GMVportfolio) / length(risk_GMVportfolio)
avg.cov <- Reduce("+", cov_gmvport) / length(cov_gmvport)
avg.portfolio.stats <- Reduce("+", portfolio.stats ) / length(portfolio.stats)

risk.GMVport
avg.cov
avg.portfolio.stats

# GMV Portfolio weight and portfolio return
# Portfolio weights are found for each replicate (Method 1)
w_gmv1 <- as.matrix(risk.GMVport["GMV_Portfolio_Wt", ])

# Portfolio return averaging over M portfolio returns in the list
portfolio.ret1 <- as.numeric(avg.portfolio.stats["port.ret", ])
portfolio.ret1

portfolio.sd1 <- as.numeric(avg.portfolio.stats["port.sd", ])
portfolio.sd1

portfolio.SR1 <- portfolio.ret1 / portfolio.sd1
portfolio.SR1

# Using average of the portfolio vector and average of the mean return vector
#avg.portfolio.ret <- as.numeric(risk.GMVport["GMV_Portfolio_Wt", , drop=F]%*%t(risk.GMVport["Mean", , drop=F]))
#avg.portfolio.ret

# Method 2
inv.avg.cov <- solve(avg.cov)

ones <- matrix(c(1, 1), nrow=ncol(dat.xts), ncol=1)
numerator <- inv.avg.cov%*%ones
denominator <- t(ones)%*%inv.avg.cov%*%ones
w_gmv2 <- numerator / as.numeric(denominator)
w_gmv2

# Using the average of the mean return vector and the portfolio vector
portfolio.ret2 <- as.numeric(risk.GMVport["Mean", , drop=F]%*%w_gmv2)
portfolio.ret2

portfolio.sd2 <- sqrt(as.numeric(t(w_gmv2)%*%avg.cov%*%w_gmv2))
portfolio.sd2

portfolio.SR2 <- portfolio.ret2 / portfolio.sd2
portfolio.SR2

# Method 3
# Invert each cov matrix in the list
inv.cov_gmvport <- lapply(cov_gmvport, solve)
avg.inv.cov <- Reduce("+", inv.cov_gmvport) / length(inv.cov_gmvport)

ones <- matrix(c(1, 1), nrow=ncol(dat.xts), ncol=1)
numerator <- avg.inv.cov%*%ones
denominator <- t(ones)%*%avg.inv.cov%*%ones
w_gmv3 <- numerator / as.numeric(denominator)
w_gmv3

# Using the average of the mean return vector and the portfolio vector
portfolio.ret3 <- as.numeric(risk.GMVport["Mean", , drop=F]%*%w_gmv3)
portfolio.ret3

portfolio.sd3 <- sqrt(as.numeric(t(w_gmv3)%*%solve(avg.inv.cov)%*%w_gmv3))
portfolio.sd3

portfolio.SR3 <- portfolio.ret3 / portfolio.sd3
portfolio.SR3

# For the full dataset
w_gmv4 <- constructGMVPortfolio(full.data)
w_gmv4

portfolio.ret4 <- as.numeric(colMeans(full.data)%*%w_gmv4)
portfolio.ret4

portfolio.sd4 <- sqrt(as.numeric(t(w_gmv4)%*%cov(full.data)%*%w_gmv4))
portfolio.sd4

portfolio.SR4 <- portfolio.ret4 / portfolio.sd4
portfolio.SR4


# Three GMV weights from three methods
cbind(w_gmv1, w_gmv2, w_gmv3, w_gmv4)
port.returns <- rbind(portfolio.ret1, portfolio.ret2, portfolio.ret3, portfolio.ret4)
row.names(port.returns) <- c("Method 1", "Method 2", "Method 3", "Original Data")

port.sds <- rbind(portfolio.sd1, portfolio.sd2, portfolio.sd3, portfolio.sd4)
row.names(port.returns) <- c("Method 1", "Method 2", "Method 3", "Original Data")

port.SRs <- rbind(portfolio.SR1, portfolio.SR2, portfolio.SR3, portfolio.SR4)
row.names(port.returns) <- c("Method 1", "Method 2", "Method 3", "Original Data")

comparison.mat <- cbind(port.returns, port.sds, port.SRs)
colnames(comparison.mat) <- c("Return", "Std Dev", "Sharpe Ratio")
comparison.mat

