# context("Basic Theory")
#
# test_that("convolving a step signal is kosher", {
#   x <- c(rep(0, 30), rep(1, 30), rep(0, 30), rep(1, 30), rep(0, 30))
# })

## library(biosignals)
## fp <- '/Users/stavros/cBio/biosignals/biosignals-pkg/inst/extdata/fret.txt'
## x <- as.numeric(readLines(fp))
## k <- biosignals:::.gaussianKernel1D(2, bandwidth=20, deriv=0)
## ke <- biosignals:::.gaussianKernel1D(8, bandwidth=20, deriv=1)

## xs <- convolve1d(x, k, rescale=TRUE)
## xe <- convolve1d(x, ke, rescale=FALSE)

## plot(x, type='l', lwd=2)
## lines(xs, col='blue')

## ## Test IRanges islands
## quartz()
## starts <- c(1, 256, 771)
## ends <- c(186, 751, 840)

## xs.fp <- convolve1d.fp(x, k, starts, ends, rescale=FALSE)

## plot(x, type='l')
## lines(xs.fp, col='blue')
