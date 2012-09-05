setGeneric("detectPeaksByEdges",
function(x, bandwidth=35, mu=0, sd=1, min.height=1L, ...) {
  standardGeneric("detectPeaksByEdges")
})

##' Returns NULL if there was an error
##'
setMethod("detectPeaksByEdges", c(x="numeric"),
function(x, bandwidth, mu, sd, min.height, ignore.from.start=0L,
         ignore.from.end=0L, ...) {
  if (is.numeric(min.height) && min.height > 0) {
    peaks <- detectPeaksByEdges(Rle(x), bandwidth, mu, sd,
                                min.height=min.height,
                                ignore.from.start=ignore.from.start, ...)
    return(peaks)
  }

  edges <- detectEdges(x, bandwidth=bandwidth, mu=mu, sd=sd,
                       threshold=threshold, ...)

  ## --------------------------------------------------------------------------
  ## Ensure that no edges are picked up within the "ignore.from.*" bounds.
  if (is.numeric(ignore.from.start) && ignore.from.start > 0) {
    if (length(edges$start) > 0) {
      edges$start <- edges$start[edges$start >= ignore.from.start]
    }
    if (length(edges$end) > 0) {
      edges$end <- edges$end[edges$end >= ignore.from.start]
    }
  }


  if (is.numeric(ignore.from.end) && ignore.from.end > 0) {
    if (length(edges$start) > 0) {
      edges$start <- edges$start[edges$start <= length(x)-ignore.from.end + 1L]
    }
    if (length(edges$end) > 0) {
      edges$end <- edges$end[edges$end <= length(x) - ignore.from.end + 1L]
    }
  }

  ## ---------------------------------------------------------------------------
  ## Cleaning vagaries in edge calling that shouldn't happen when looking for
  ## peaks

  ## Look for and remove edges that start after the last end, or ones that end
  ## before the first start
  if (length(edges$start) > 0 && length(edges$end) > 0) {
    bad.start <- edges$start >= max(edges$end)
    if (any(bad.start)) {
      edges$start <- edges$start[!bad.start]
    }
    bad.end <- edges$end <= min(edges$start)
    if (any(bad.end)) {
      edges$end <- edges$end[!bad.end]
    }
  }

  ## If edges are out of order (insane for peaks), then try to clean it up
  fishy <- .bad.edge.detection(edges$start, edges$end) ||
    !.start.end.inorder(edges$start, edges$end)

  ## TODO: To ignore noise -- add "edge refinement" code here.
  if (fishy) {
    return(NULL)
  }
  
  ## Make edges start at non-zero positions
  shift.starts <- x[edges$start] == 0 & x[edges$start + 1L] > 0L
  if (any(shift.starts)) {
    edges$start[shift.starts] <- edges$start[shift.starts] + 1L
  }

  shift.ends <- x[edges$end] == 0 & x[edges$end - 1L] > 0L
  if (any(shift.ends)) {
    edges$end[shift.ends] <- edges$end[shift.ends] - 1L
  }

  IRanges(edges$start, edges$end)
})


## For smooth.slice code, see commits between [Dec 21, Dec 23]
## smooth.slice tries (tired) to avoid the following situation where we call
## two peaks where one should be:
##
##          /
##      .  /|
## ----/|-/ |------    [ <- min.height line ]
##    . |/  |
##   /      |
##
setMethod("detectPeaksByEdges", c(x="Rle"),
function(x, bandwidth, mu, sd, min.height, trim=c(0.05, 0.95),
         ignore.from.start=0L, ignore.from.end=0L, ...) {
  if (is.numeric(trim)) {
    if (length(trim) != 2 || any(trim <= 0) || any(trim >= 1)) {
      stop("`trim` should be vector of length two between (0,1)")
    }
  }

  bandwidth <- as.integer(bandwidth)[1L]
  min.height <- as.integer(max(1L, min.height))
  F <- getMethod('detectPeaksByEdges', 'numeric')

  islands <- slice(x, lower=min.height, rangesOnly=TRUE)

  edges <- lapply(1:length(islands), function(i) {
    istart <- start(islands[i])
    iend <- end(islands[i])
    iwidth <- iend - istart + 1L
    ix <- as.numeric(x[istart:iend])

    if (is.numeric(trim)) {
      qpos <- quantilePositions(ix, quantile.breaks=trim)
      qstart <- qpos[[1L]][1L]
      qend <- qpos[[2L]][1L]
      qx <- ix[qstart:qend]
    } else {
      qstart <- 1L
      qend <- length(ix)
      qx <- ix
    }

    xx <- c(rep(0L, bandwidth), qx, rep(0L, bandwidth))
    e <- F(xx, bandwidth, mu, sd, min.height=0L,
           ignore.from.start=max(bandwidth - qstart, 0L, ignore.from.start),
           ignore.from.end=max(bandwidth - length(ix) + qend, 0L, ignore.from.end), ...)

    ## NULL is returned on error/out-of-bounds conditions, in this case we set
    ## the "peak" to just be the (trimmed) fenceposts of this "island"
    if (is.null(e) || length(e) == 0L) {
      if (qstart < ignore.from.start || qend > ignore.from.end) {
        ## TODO: This needs to be refactored.
        ## If the inferred peak bounds are out of range given the ignore*
        ## params, we will junk this peak
        return(NULL)
      } else {
        e <- IRanges(qstart, qend)
        values(e) <- DataFrame(fishy=TRUE)
      }
    } else {
      e <- shift(e, -bandwidth + qstart)
      values(e) <- DataFrame(fishy=rep(FALSE, length(e)))
    }

    values(e)$multi.modal <- rep(length(e) > 1, length(e))
    values(e)$island.idx <- i

    ## Shift the edge calls up into the correct region of x
    shift(e, istart - 1L)
  })

  edges <- do.call(c, unname(edges[!sapply(edges, is.null)]))
  edges
})

##' Attempts to remove "noisy" peaks [NOT USED]
##'
##' @param x The numeric vector used to call peaks from
##' @param edges The IRanges object which identifies the start/end edges
##' of a peak.
##' @param peak.threshold The fraction that the minimum y-value at
##' either end of the peak must be less than in order for the peak
##' to be considered "real"
##' @param threshold.maxs Leave this NULL if each peak's height with respect
##' to its minima (foudn at the edges) should be evaluated by themselves.
##' Pass in a number to use a a more "global" maximum for each peak maximum
##' to be compared to.
filterPeaks <- function(x, edges, min.height=1L, peak.threshold=0.75,
                        edge.window=5L, threshold.maxs=NULL) {
  stop("This isn't working yet -- jump to the DEBUG/threshold stuff")

  if (edge.window < 1) {
    stop("edge.window must be non-negative")
  }
  edge.window <- as.integer(ceiling(edge.window))

  if (!(is.integer(as.vector(x[1L])))) {
    stop("An integer-type of vector is required (vector or Rle)")
  }
  if (any(start(ranges) < 1 | end(ranges) > length(x))) {
    stop("Illegal start/end values in `edges` IRanges object")
  }
  if (peak.threshold < 0 || peak.threshold > 1) {
    stop("peak.threshold must be numeric in [0,1] range")
  }
  if (peak.threshold == 0) {
    return(IRanges())
  }
  if (peak.threshold == 1) {
    return(edges)
  }

  views <- Views(x, edges)
  maxs <- viewMaxs(views)
  keep.height <- maxs >= min.height

  edges <- edges[keep.height]
  maxs <- maxs[keep.height]

  ## DEBUG: This thresholding code is not correct
  if (is.null(threshold.maxs)) {
    start.window <- IRanges(start(edges) - edge.window + 1L, width=edge.window)
    end(start.window) <- pmax(1L, end(start.window))
    start(start.window) <- pmax(1L, start(start.window))
    val.start <- viewMins(Views(x, start.window))

    end.window <- IRanges(end(edges) + edge.window - 1L, width=edge.window)
    start(end.window) <- pmin(length(x), start(end.window))
    end(end.window) <- pmin(length(x), end(end.window))
    val.end <- viewMins(Views(x, end.window))

    edge.min <- pmin(val.start, val.end)
    keep <- (edge.min + 0.5 / maxs) <= peak.threshold & maxs > 0
  } else {
    ## ??
    keep <- maxs / threshold.maxs <= peak.threshold & maxs > 0
  }


  edges[keep]
}

###############################################################################
## Helper functions to check that peak calling works
## ----------------------------------------------------------------------------
.start.end.inorder <- function(starts, ends) {
  is.start <- c(rep(TRUE, length(starts)), rep(FALSE, length(ends)))
  check <- is.start[order(c(starts, ends))]
  all(rle(check)$lengths == 1)
}

.bad.edge.detection <- function(starts, ends) {
  length(starts) == 0L || length(starts) != length(ends) ||
  min(ends) <= min(starts) || max(starts) >= max(ends)
}

###############################################################################
## Sidlined utility functions that were used for an elementary effort of doing
## adaptive bandwidth selection in local "coverage islands". There is some
## real literature out there that talks about this, so use some of those
## methods instead.

## This is a utility function used for "adaptive bandwidth" hunting.
##
## This functionality is temporarily removed until I can smoke it out better
validateBandwidthAndWindows <- function(bandwidth, window.size) {
  if (length(bandwidth) < 1) {
    stop("bandwidth required")
  }

  if (length(bandwidth) > 1L) {
    o <- order(bandwidth, decreasing=TRUE)
    bandwidth <- bandwidth[o]
    if (length(window.size) == 1L) {
      window.size <- rep(window.size, length(bandwidth))
    } else {
      if (length(window.size) != length(bandwidth)) {
        stop("length of edge.window vector != length of bandwidth vector")
      }
      window.size <- window.size[o]
    }
  }

  if (length(window.size) != length(bandwidth)) {
    stop("length of window size must equal bandwidth length")
  }

  bandwidth <- as.integer(bandwidth)
  window.size <- as.integer(window.size)

  if (any(bandwidth) <= 0) {
    stop("bandwidth(s) must be > 0")
  }

  list(bandwidth=bandwidth, window.size=window.size)
}

##' Refines the edge calls in pre-specified regions of x.
##'
##' This function is used during "adaptive bandwidth" hunting in local
##' coverage islands -- this functionality is temporarily removed until
##' I can implement it better.
##'
##' @param x The coverage vector.
##' @param bandwidth Bandwidth of the kernel
##' @param mu The man of the kernel
##' @param sd
refineEdges <- function(x, bandwidth, mu, sd, starts, ends,
                        half.window=ceiling(bandwidth * 2),
                        .strand='+', max.iter=10,
                        detangle.by=c('max', 'min'), ...) {
  detangle.by <- match.arg(detangle.by)

  bw <- validateBandwidthAndWindows(bandwidth, half.window)
  bandwidth <- bw$bandwidth
  half.window <- bw$half.window

  if (length(starts) == 0 || length(ends) == 0) {
    obs.quantiles <- quantilePositions(x, IRanges(1, length(x)))
  } else {
    obs.quantiles <- data.frame(start=1, end=length(x))
  }

  if (length(starts) == 0) {
    starts <- obs.quantiles[[1L]][1L]
  }
  if (length(ends) == 0) {
    ends <- obs.quantiles[[ncol(obs.quantiles)]][[1L]]
  }

  ## Are the starts and ends order out of whack?

  if (!.start.end.inorder(starts, ends)) {
    if (detangle.by == 'max') {
      starts <- min(starts)
      ends <- max(ends)
    } else {
      starts <- max(starts)
      ends <- min(ends)
    }
  }

  list(start=starts, end=ends)
}

