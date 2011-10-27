asRle <- function(x, eps=1e-6) {
  stopifnot(is.numeric(x))
  ret <- .Call("Ras_rle", x, eps, PACKAGE="biosignals")
  Rle(ret$values, ret$lengths)
}

expandRle <- function(x) {
  stopifnot(is(x, "Rle"))
  stopifnot(is.numeric(as.vector(x[1])))
  .Call("Rexpand_rle", runLength(x), runValue(x), PACKAGE="biosignals")
}

slidingMax <- function(x, k=5L) {
  ret <- .Call('Rsliding_max', as.numeric(x), as.integer(k),
               PACKAGE="biosignals")
}

slidingMin <- function(x, k=5L) {
  .Call('Rsliding_min', as.numeric(x), as.integer(k), PACKAGE="biosignals")
}

zeroCrossings <- function(x) {
  .Call('Rzero_crossings', as.numeric(x), PACKAGE="biosignals")
}

quantilePositions <- function(signal, bounds=IRanges(1, length(signal)),
                              quantile.breaks=NULL) {
  if (!inherits(signal, "Rle")) {
    ## TODO: Add slice,numeric to IRanges
    signal <- Rle(signal)
  }
  if (is.numeric(bounds) && length(bounds) == 1) {
    bounds <- slice(signal, bounds, includeLower=bounds > 0)
  }
  stopifnot(is(bounds, 'IRanges'))
  if (missing(quantile.breaks) || is.null(quantile.breaks)) {
    quantile.breaks <- c(0.05, 0.10, 0.15, 0.20, 0.80, 0.85, 0.90, 0.95)
  }
  if (min(quantile.breaks) <= 0 || max(quantile.breaks) >= 1) {
    stop("Quantiles must be within (0,1) range.")
  }
  if (is.unsorted(quantile.breaks)) {
    stop("Quantiles need to be in ascending order")
  }
  if (any(start(bounds) < 1L) || any(end(bounds) > length(signal))) {
    stop("Bounding info is out of bounds.")
  }
  signal <- as.numeric(signal)

  quants <- .Call("Rcoverage_quantiles", signal, start(bounds), end(bounds),
                  quantile.breaks, PACKAGE="biosignals")
  colnames(quants) <- paste('quantile', gsub("0\\.", "", quantile.breaks), sep=".")
  as.data.frame(quants)
}

