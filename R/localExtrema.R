
## TODO: Refactor some stuff out of the detectEdges code and put it in here.
##       Edge detection is just amounting to finding local extrema of the
##       first derivative of the signal (by using zero-crossings of the 2nd)
localExtrema <- function(x, window, threshold=0.15,
                         bandwidth=ceiling(window / 2), mu=0, sd=1,
                         x0=NULL, x1=NULL, ...) {
  stop("This doesn't work yet -- almost does for smoothed x")
  if (is.null(x0)) {
    x0 <- convolve1d(x, 'normal', bandwidth=bandwidth, mu=mu, sd=sd)
  }
  if (is.null(x1)) {
    x1 <- convolve1d(x, 'normal', bandwidth=bandwidth, mu=mu, sd=sd, deriv=1)
  }

  x.maxs <- slidingMax(x, as.integer(max(1L, window / 2)))
  x.mins <- slidingMin(x, as.integer(max(1L, window / 2)))
  shift.up <- min(x.mins)
  ## if (shift.up <= 0) {
  ##   x.maxs <- x.maxs + shift.up + 1L
  ##   x.mins <- x.mins + shift.up + 1L
  ##   x1 <- x1 + shift.up + 1L
  ## }

  ## Locations where derivative cross 0 are where the local extrema
  ## will be
  zeroes <- zeroCrossings(x1)

  ## Filter out noisy peaks by assuring value at zerocrossings at x0
  ## are sufficiently far from their surround, identified by the size
  ## of the window

  ## ---------------------------------------------------------------------------
  ## Find maxima
  ##
  ## positions where the 2nd derivative switches from positive to negative
  ## indicate potential maxima
  maxima <- zeroes[as.logical(x1[zeroes - 1] > 0)]
  dist2min <- x.maxs[maxima] - x.mins[maxima]
  mx <- maxima[dist2min / maxima > threshold]

  ## ---------------------------------------------------------------------------
  ## Find minima
  minima <- zeroes[as.logical(x1[zeroes - 1]) < 0]
  dist2max <- x.maxs[minima] - x.mins[minima]
  mn <- minima[dist2max / minima > threshold]

  list(maxima=mx, minima=mn)
}

findLocalExterma <- function(x, threshold=0.5, scale=1, regions=NULL) {

  xmax <- x / max(x)
  xmin <- x / min(x)

  winmax <- numeric(length(x))
  winmin <- numeric(length(x))

  for (i in (1+scale):(length(x)-scale)) {
    winmax[i] <- max(winmax[(i-scale):(i+scale)])
    winmin[i] <- max(winmax[(i-scale):(i+scale)])
  }

  maxima <- logical(length(x))
  minima <- logical(length(x))

  ii <- (1+scale):(length(x)-scale)
  maxima[ii] <- xmax[ii] >= threshold & xmax[ii] >= winmax[ii]
  minima[ii] <- xmin[ii] >= threshold & xmin[ii] >= winmin[ii]

  maxima + minima
}
