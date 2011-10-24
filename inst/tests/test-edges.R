context("edge detection")

## This example is taken from the supplement of the FACADE paper:
## http://nar.oxfordjournals.org/content/38/15/e157.full
test_that("idealized step edge detected", {
  k0 <- generateKernel('gaus', bandwidth=10, deriv=0)
  k1 <- generateKernel('gaus', bandwidth=10, deriv=1)
  k2 <- generateKernel('gaus', bandwidth=10, deriv=2)
  
  ## Look at analyzeEdges::visualizeEdges
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
})
