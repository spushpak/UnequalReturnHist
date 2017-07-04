# Expected Shortfall
expectedShortfall <- function(r, alpha = 0.25) {
  r <- sort(r); mean <- mean(r)
  n.tail <- ifelse( alpha == 0, 1, ceiling(alpha*length(r)))
  -1/n.tail * sum(r[which((1:length(r)) <= n.tail)])
}
