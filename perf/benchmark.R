library(KernSmooth)
library(microbenchmark)
x <- 2 * pi * runif(100000)
y <- sin(x) + rnorm(100000) / 10
v <- seq(from = 0, to = 2 * pi, length.out = 1000)
h <- dpill(x, y, gridsize = 1000, range.x = c(0, 2 * pi))
microbenchmark("Local linear" = {locpoly(x, y, bandwidth = h, gridsize = 1000, range.x = c(0, 2 * pi))})