##' Smoothes vector \code{x} using kernel density smoothing
##'
##' @author Anshul Kundaje (original MATLAB implementation)
##' @author Steve Lianoglou (ported to R)
##'
##' @param x input vector
##' @param kernel.type name of kernel
##' @param bandwidth kernel half width
##' @param normalize If set to true, kernel is equivalent to smooth averaging,
##' otherwise it is equivalent of a smooth sum
##' @param trim.length \code{logical(1)}. The vector returned by
##' \code{\link{convolve}} might be longer than the one you gave it. If
##' \code{TRUE}, the returned vector will be as long as your original vector.
##'
##' @return a vector of smoothed values for \code{x}
kdeSmooth <- function(x, kernel, bandwidth=1L, normalize=FALSE,
                      trim.length=TRUE, rescale=TRUE, upper.quantile=0.98,
                      convolve.missing=FALSE, kernel.norm=NULL, ...) {
  if (inherits(x, 'Rle')) {
    return(kdeSmooth.Rle(x, kernel, bandwidth, normalize=normalize,
                         trim.length=trim.length, ...))
  }
  stopifnot(is.numeric(x))
  if (!is.null(dim(x))) {
    ## No convultion of matrices yet
    stop("x can only be 1D")
  }
  
  if (is.numeric(kernel) && missing(bandwidth)) {
    bandwidth <- ceiling(length(kernel) / 2)
  }
  
  if (is.character(kernel)) {
    kernel <- match.arg(kernel, .convolutionKernels())
    kernel.norm <- generateKernel(kernel, bandwidth, normalize=TRUE, ...)
    kernel <- generateKernel(kernel, bandwidth, normalize=normalize, ...)
  }
  stopifnot(is.numeric(kernel))
  
  ## Pad head/tail of vector to sidestep "edge effects"
  x.padded <- c(rep(x[1], bandwidth), x, rep(tail(x, 1), bandwidth))
  not.number <- !is.finite(x.padded)
  has.nan <- any(not.number)

  convolve.missing <- convolve.missing && !is.null(kernel.norm)
  if (convolve.missing) {
    ## Guard against missing values
    if (has.nan) {
      missing.values <- x.padded[not.number]
      x.padded[not.number] <- 0
    }
  }

  xs <- convolve(x.padded, rev(kernel), type='open')

  if (convolve.missing) {
    ## Dividing by normalized result "deals with" missing values in x
    xs <- xs / convolve(as.numeric(!not.number), rev(kernel.norm), type='open')
  }

  ## Shift the signal "back" by to account for extra padding + bandwidth
  xs <- xs[-seq(2*bandwidth)]
  if (trim.length) {
    ## The vector you get might be longer than the one you gave.
    ## xs <- xs[1:min(length(x) + bandwidth, length(xs))]
    xs <- xs[1:min(length(x), length(xs))]
  }

  ## Reset original missing values
  if (convolve.missing && has.nan) {
    xs[not.number] <- missing.values
  }

  ## TODO: Consider looking for + removing outliers at the tails of vector.
  ##       Maybe look at the range of `diff`s you get at the edges and remove
  ##       outliers. Maybe the lenght of tails to look for is a function of
  ##       bandwidth -- explore further
  if (rescale) {
    xs <- (xs / max(xs, na.rm=TRUE)) * max(x, na.rm=TRUE)
  }

  xs
}

kdeSmooth.Rle <- function(x, kernel, bandwidth, normalize=TRUE,
                          trim.length=TRUE, rescale=TRUE, upper.quantile=0.98,
                          ...) {
  stopifnot(inherits(x, 'Rle'))
  if (is.character(kernel)) {
    kernel <- match.arg(kernel, .convolutionKernels())
    kernel <- generateKernel(kernel, bandwidth, normalize=normalize, ...)
  }

  regions <- slice(x, 1, rangesOnly=TRUE)
  v <- Views(x, regions)

  ## Calculate new values in runs
  vs <- viewApply(v, function(vv) {
    vector <- as.numeric(vv)
    ## conv <- convolve(padded, kval, type='open')
    conv <- kdeSmooth(vector, kernel, bandwidth, rescale=rescale,
                      trim.length=TRUE)
    conv
  }, simplify=FALSE)

  .as.Rle(vs, regions, length(x))
}


.fixEdgeEffects <- function(x, bandwidth, is.normed=FALSE) {
  if (is.normed) {
    xn <- x
  } else {
    xn <- xn / max(xn)
  }

  edge.size <- ceiling(length(x) * 0.1)
  ## head
  ## idx.head <- 1:edge.size
  ## diff.head <- diff(x[idx.head])
  ## fix.head <- max(diff.head) > 1.5 * sd(diff.head)
  ## if (any(fix.head)) {
  ##   x[fix.head] <-
  ## }

  ## idx.tail <- (length(x) - edge.size + 1L):length(x)
  ## diff.tail <- diff(x[idx.tail])
  ## fix.tail <- max(diff.tail) > 1.5 * sd(diff.tail)
  ## if (any(fix.tail)) {

  ## }
}

.as.Rle <- function(new.vals, val.ranges, len) {
  if (is.unsorted(start(val.ranges)) || is.unsorted(end(val.ranges))) {
    stop("Illegal ranges")
  }
  o <- findOverlaps(val.ranges, ignoreSelf=TRUE, ignoreRedundant=TRUE,
                    type="any")
  if (length(o) > 0L) {
    stop("Overlapping ranges not handled yet")
  }

  prev.end <- 0L
  runs <- lapply(1:length(new.vals), function(idx) {
    ## Builds 2 column matrix. Col 1 is runValue, col2 is runLength
    val <- new.vals[[idx]]
    .range <- val.ranges[idx]

    ans <- cbind(val, rep(1, length(val)))
    dist <- start(.range) - prev.end -1L
    if (dist > 0) {
      ans <- rbind(c(0, dist), ans)
    }
    prev.end <<- end(.range)
    ans
  })

  end.point <- end(val.ranges[length(val.ranges)])
  if (end.point < len) {
    runs[[length(runs) + 1L]] <- c(0, len - end.point + 1)
  }
  runs <- do.call(rbind, runs)

  smoothed <- Rle(runs[,1], runs[,2])
  smoothed
}
