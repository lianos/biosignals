context("edge detection")

## This example is taken from the supplement of the FACADE paper:
## http://nar.oxfordjournals.org/content/38/15/e157.full
test_that("idealized step edge detected", {
  k0 <- generateKernel('gaus', bandwidth=10, deriv=0)
  k1 <- generateKernel('gaus', bandwidth=10, deriv=1)
  k2 <- generateKernel('gaus', bandwidth=10, deriv=2)
  
  ## Visualize the kernes
  if (FALSE) {
    kdog <- generateKernel('dog', bandwidth=10, deriv=2)
    ylim=c(min(k0, k1, k2), max(k0, k1, k2))
    plot(k0, type='l', lwd=2, ylim=ylim)
    lines(k1, lwd=2, col='red')
    lines(k2, lwd=2, col='green')
    z0 <- zeroCrossings(k2)
    abline(v=z0, col='grey', lwd=2, lty='dashed')
  }
  
  x <- c(rep(0, 15), rep(20, 15), -1, rep(0, 20))
  x <- x + rnorm(length(x))
  x0 <- convolve1d(x, k0)
  x1 <- convolve1d(x, k1)
  x2 <- convolve1d(x, k2)
  
  quartz()
  
  ylim <- c(min(x, x0, x1, x2), max(x, x0, x1, x2))
  plot(x, type='l', ylim=ylim, col='blue', lwd=2)
  lines(x0, col='grey', lty='dashed', lwd=2)
  lines(x1, col='red', lwd=2)  
  lines(x2, col='green', lwd=2)  
  zz <- zeroCrossings(x2)
  abline(v=zz, col='grey', lwd=2, lty='dashed')
  abline(h=0)
})
