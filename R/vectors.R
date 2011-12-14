##' Find the maximum values in sliding window
##'
##' @param x The vector of values to search
##' @param k The size of the window to search in. The window is size k*2
##' @return vector of \code{length(x)}, indicating the maximum value found
##' in the window centered at the given position and includes k bins up
##' and downstream.
slidingMax <- function(x, k=5L) {
  .Call('Rsliding_max', as.numeric(x), as.integer(k), PACKAGE="biosignals")
}

##' Find minimum values in sliding window
##'
##' @param x The vector of values to search
##' @param k The size of the window to search in. The window is size k*2
##' @return vector of \code{length(x)}, indicating the minimum value found
##' in the window centered at the given position and includes k bins up
##' and downstream.
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
  colnames(quants) <- paste('quantile', gsub("0\\.", "", quantile.breaks),
                            sep=".")
  as.data.frame(quants)
}

## ----------------------------------------------------------------------------
## These Rle functions are just there to test if the C++ functions are working
## They are not exported
asRle <- function(x, eps=1e-6) {
  stopifnot(is.numeric(x))
  .Call("Ras_rle", x, eps, PACKAGE="biosignals")
}

expandRle <- function(x) {
  stopifnot(is(x, "Rle"))
  stopifnot(is.numeric(as.vector(x[1])))
  .Call("Rexpand_rle", runLength(x), runValue(x), PACKAGE="biosignals")
}

expandRleS4 <- function(x) {
  stopifnot(is(x, "Rle"))
  .Call("Rexpand_rle_S4", x, PACKAGE="biosignals")
}

## ----------------------------------------------------------------------------
