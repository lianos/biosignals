# context("Basic Theory")
# 
# test_that("convolving a step signal is kosher", {
#   x <- c(rep(0, 30), rep(1, 30), rep(0, 30), rep(1, 30), rep(0, 30))
# })

library(edge1d)
x <- as.numeric(readLines('/Users/stavros/cBio/edge1d/edge1d-pkg/inst/extdata/fret.txt'))
k <- edge1d:::.gaussianKernel1D(2, bandwidth=20, deriv=0)
xs <- convolve1d(x, k)

plot(x, type='l', lwd=2)
lines(xs, col='blue')
