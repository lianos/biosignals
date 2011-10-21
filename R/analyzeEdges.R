##' Analyzes 1D data using multiscale analysis.
##' 
##' Finds distinct transition regions and transition ratios.
##' 
##' @param x The 1D data
##' @param scales A vector of scales at which to perform analysis
##' @param threshold A vector of numbers between [0, 1] indicating
##' how "obvious" a step has to be at each scale in order to be
##' considered a transition
##' @param timestep The amount of time one data point represents
##' @param start.time The time at which to start analysis
##' @param end.time The time at which to end analysis
##' @param tran.rad The number of samples to ignore on either side
##' of a transition point when calculating region statistics
##' 
##' @return $dData deriviative scale space,
##' $minmax: a vector that indicates minima and maxima within the data
##' $stats A table of mean and std.dev for each distinct data region
analyzeEdges <- function(x, scales=c(1, 2, 4, 8),
                         tresholds=c(0.1, 0.2, 0.3, 0.4),
                         timestep=1L, start.time=1L, end.time=length(x),
                         tran.rad=1) {
  
  ## Limit the analysis to time region specified?
  if (missing(timestep) || missing(start.time) || missing(end.time)) {
    time <- 1:length(x)
  }
  
  d.data <- createGaussScaleSpace(x, 1, scales)
  
  ## Find the position of local min/max at the most coarse scale
  min.max <- findLocalExtrema(d.data[, ncol(d.data)], tail(thresholds, 1),
                              tail(scales, 1))
}

## To find "real" edges, look for zero-crossings of the gauss,deriv=2 kernel
## and filter this using "clear peaks" found in gauss,deriv=1 kernel
## 
## A "clear peak" has to do with the y-value of the peak's center. This will
## be large -- but what is "large"? This is where you use the "local maximum"
## combined w/ threshold, somehow.
visualizeEdges <- function(x, mu=0, sd=1, bandwidth=10, threshold=0.15) {
  if (missing(x)) {
    ## this is a test
    x <- c(rep(0, 15), rep(20, 15), -1, rep(0, 20))
    x <- x + rnorm(length(x))
  }
  
  k0 <- generateKernel('gaus', bandwidth=bandwidth, mu=mu, sd=sd, deriv=0)
  k1 <- generateKernel('gaus', bandwidth=bandwidth, mu=mu, sd=sd, deriv=1)
  k2 <- generateKernel('gaus', bandwidth=bandwidth, mu=mu, sd=sd, deriv=2)
  
  x0 <- convolve1d(x, k0)
  x1 <- convolve1d(x, k1)
  x2 <- convolve1d(x, k2)
  
  ylim <- c(min(x, x0, x1, x2), max(x, x0, x1, x2))
  plot(x, type='l', ylim=ylim, col='blue', lwd=2)
  lines(x0, col='blue', lty='dashed', lwd=2)
  lines(x1, col='red', lwd=2)  
  lines(x2, col='green', lwd=2)
  
  zz <- zeroCrossings(x2)
  significant <- which(abs(x1) > threshold * max(abs(x1)))
  zz <- intersect(zz, significant)
  abline(v=zz, col='grey', lwd=2, lty='dashed')
  abline(h=0)
  invisible(x)
}
