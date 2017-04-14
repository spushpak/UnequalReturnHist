rm(list=ls())
setwd("C:/UW/RDev/CBF")

source("cbf_functions.R")
library(xts)
library(moments)
library(PerformanceAnalytics)

# Read data
#dat.mat <- read.csv("Wilshire and MSCI EM.csv", header = T, stringsAsFactors = F)
#dat.mat <- read.csv("multiple_grps.csv", header = T, stringsAsFactors = F)
dat.mat <- read.csv("hfunds5_ue_ts.csv", header = T, stringsAsFactors = F)

# Convert the data to xts object (probably it's not necessary)
dat.xts <- xts(dat.mat[, -1], order.by = as.Date(dat.mat[, 1], "%m/%d/%Y"))

# Keep a copy of original data
old.xts <- dat.xts 

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

# Create rgression expression for 'lm()' function
regressor.list <- paste(long.hist.var, collapse = "+")

# Full length of the dataset
full.length <- nrow(dat.xts)

new.dat <- dat.xts[, long.hist.var]

# Start with the short history columns which have the smallest number of 
# observations missing i.e. this group is the longest short history group. Then 
# move on to the next shorter history group and so on.
for (i in miss.itr) {
  short.hist.var <- rownames(na_count)[na_count[, "num.miss"]==i]
  temp.dat <- dat.xts[, c(long.hist.var, short.hist.var), drop=F]
  num.miss <- i
  reg.dat <- fitValue(short.hist.var, long.hist.var, temp.dat, num.miss)
  new.dat <- uneqhistCBF(reg.dat, new.dat)  
} ############# end of for loop (for each short history group)


################################################################################
# # Create the start index of each repeating block in the new dataset. Each repeating
# # block has length equal to the row numbers of initial dataset.
# row.grp <- seq(1, nrow(new.dat), by=full.length)
# 
# # Summary statistics of the combined backfilled dataset
# #desc.stats(new.dat)
# 
# # Summary statistics of the original dataset
# #desc.stats(dat.xts)
# 
# # Total length of the newly constructed data
# prod(full.length-miss.itr)*full.length
# 
# # No. of repeatng blocks
# num.block <- prod(full.length-miss.itr)
# 
# # Create an empty list
# risk.metrics <- vector("list", num.block)
# 
#   # Repeat the above step for each repaeting block 
#   for (j in 1:length(row.grp)) {
#     risk.mat <- matrix(0, nrow=6, ncol=length(miss.hist.var), 
#                        dimnames=list(c("Skewness", "Kurtosis", "Mean", "Volatility", "Sharpe Ratio", "Expected Shortfall"), miss.hist.var))
#     
#     mom.dat <- new.dat[row.grp[j]:(row.grp[j]+full.length-1), miss.hist.var, drop=F]
# 
#     risk.mat["Skewness", ] <- moments::skewness(mom.dat)
#     risk.mat["Kurtosis", ] <- moments::kurtosis(mom.dat)
#     risk.mat["Mean", ] <- colMeans(mom.dat)
#     risk.mat["Volatility", ] <- apply(mom.dat, 2, sd)
#     #risk.mat["Sharpe Ratio", ] <- SharpeRatio(mom.dat, FUN="StdDev")
#     risk.mat["Sharpe Ratio", ] <- risk.mat["Mean", ]/risk.mat["Volatility", ]
#     #risk.mat["Expected Shortfall", ] <- ES(mom.dat, p=.95, method="historical")
#  
#     risk.metrics[[j]] <- risk.mat
#   }
# 
# risk.metrics

####################################################################
# Create an empty matrix to hold the risk measures
risk.vals <- matrix(0, nrow=6, ncol=length(miss.hist.var),
                    dimnames=list(c("Skewness", "Kurtosis", "Mean", "Volatility", "Sharpe Ratio", "Expected Shortfall"), miss.hist.var))

risk.vals["Skewness", ] <- moments::skewness(new.dat[, miss.hist.var])
risk.vals["Kurtosis", ] <- moments::kurtosis(new.dat[, miss.hist.var])
risk.vals["Mean", ] <- colMeans(new.dat[, miss.hist.var])
risk.vals["Volatility", ] <- apply(new.dat[, miss.hist.var], 2, sd)
risk.vals["Sharpe Ratio", ] <- risk.vals["Mean", ]/risk.vals["Volatility", ]
#risk.vals["Expected Shortfall", ] <- ES(new.dat[, miss.hist.var], p=.95, method="historical")

round(risk.vals, digits = 3)



