#' @title Implements Multiple Imputation Method folowing Page(2013)
#' 
#' @description This function implements Multiple Imputation Method folowing 
#' Page(2013)
#' 
#' @importFrom xts xts
#' @importFrom moments skewness kurtosis
#' 
#' 
#' @param dat.mat is the returns data for multiple assets with unequal return history.
#' @param saveReps is a logical flag, to indicate whether the risk measures  
#' need to be computed for all the replicates. Default is FALSE.
#' 
#' @return
#' \item{risk.list}{A list or matrix containing risk measures for either the whole dataset or for each replicate.}
#'    
#' @author Pushpak Sarkar
#' 
#' 
#' 
#' 
#' 
#' 
#' @export
#' 

# library(xts)
# library(moments)
# 
# # Read data
# #dat.mat <- read.csv("Wilshire and MSCI EM.csv", header = T, stringsAsFactors = F)
# #dat.mat <- read.csv("multiple_grps.csv", header = T, stringsAsFactors = F)
# dat.mat <- read.csv("hfunds5_ue_ts.csv", header = T, stringsAsFactors = F)


uneqhistPage <- function(dat.mat, saveReps=FALSE){

# Convert the data to xts object (probably it's not necessary)
dat.xts <- xts(dat.mat[, -1], order.by = as.Date(dat.mat[, 1], "%m/%d/%Y"))

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

risk.list <- constructPageData(fitted.xts, resid.mat, miss.hist.var, miss.itr, 
                               na_count, sReps=saveReps)

return(round(risk.list, digits = 3))

}
