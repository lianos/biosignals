setGeneric("detectPeaksByEdges",
function(x, bandwidth=40, mu=0, sd=1, threshold=0.15,
         half.window=ceiling(bandwidth / 2), min.width=10, do.clean=TRUE,
         ...) {
  standardGeneric("detectPeaksByEdges")
})

.start.end.inorder <- function(starts, ends) {
  is.start <- c(rep(TRUE, length(starts)), rep(FALSE, length(ends)))
  check <- is.start[order(c(starts, ends))]
  all(rle(check)$lengths == 1)
}

.bad.edge.detection <- function(starts, ends) {
  length(starts) == 0L || length(starts) != length(ends) ||
  min(ends) <= min(starts) || max(starts) >= max(ends)
}

##' Returns NULL if there was an error
setMethod("detectPeaksByEdges", c(x="numeric"),
function(x, bandwidth, mu, sd, threshold, half.window, min.width, do.clean,
         .strand='+', min.height=0,
         ignore.from.start=NULL, ignore.from.end=NULL, ...) {
  if (min.height > 0) {
    peaks <- detectPeaksByEdges(Rle(x), bandwidth, mu, sd, threshold,
                                half.window, min.width, do.clean, bandwidths,
                                .strand=.strand, min.height=min.height, ...)
    return(peaks)
  }

  edges <- detectEdges(x, bandwidth=bandwidth[1L], mu=mu, sd=sd,
                       threshold=threshold, half.window=half.window[1L], ...)

  ## ---------------------------------------------------------------------------
  ## Ensure that no edges are picked up within the "ignore.from.*" bounds.
  if (is.numeric(ignore.from.start)) {
    if (length(edges$start) > 0) {
      edges$start <- edges$start[edges$start >= ignore.from.start]
    }
    if (length(edges$end) > 0) {
      edges$end <- edges$end[edges$end >= ignore.from.start]
    }
  }

  if (is.numeric(ignore.from.end)) {
    if (length(edges$start) > 0) {
      edges$start <- edges$start[edges$start <= length(x) - ignore.from.end + 1L]
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

  ## if (do.clean && fishy) {
  ##   edges <- refineEdges(x, bandwidth, edges$start, edges$end,
  ##                        half.window=half.window, .strand=.strand, ...)
  ## }

  if (.bad.edge.detection(edges$start, edges$end) ||
      !.start.end.inorder(edges$start, edges$end)) {
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

  ir <- IRanges(edges$start, edges$end)
  values(ir) <- DataFrame(fishy=rep(fishy, length(ir)))
  if (min.width > 0) {
    ir <- ir[width(ir) > min.width]
  }
  ir
})

setMethod("detectPeaksByEdges", c(x="Rle"),
function(x, bandwidth, mu, sd, threshold, half.window, min.width, do.clean,
         bandwidths, pad.by=1, min.height=0L, .strand='+',
         smooth.x=TRUE, ...) {
  # if (smooth.x) {
  #   xs <- convolve1d(x, kernel='normal', bandwidth=bandwidth, mu=mu, sd=sd)
  # } else {
  #   xs <- x
  # }

  xs <- x

  f <- getMethod('detectPeaksByEdges', 'numeric')
  islands <- slice(xs, lower=min.height, rangesOnly=TRUE,
                   includeLower=min.height != 0)
  islands <- islands[width(islands) >= min.width]
  nbad <- 0L

  edges <- lapply(1:length(islands), function(i) {
    istart <- start(islands[i])
    iend <- end(islands[i])
    istart.pad <- max((istart - pad.by), 1L)
    iend.pad <- min(iend + pad.by, length(x))
    ix <- as.numeric(x[istart.pad:iend.pad])

    e <- f(ix, bandwidth, mu, sd, threshold, half.window, do.clean=do.clean,
           bandwidths=bandwidths, .strand=.strand, min.height=0,
           ignore.from.start=istart - istart.pad,
           ignore.from.end=iend.pad - iend, ...)

    ## Return a NULL on error/out-of-bounds conditions
    if (is.null(e)) {
      nbad <<- nbad + 1L
      return(NULL)
    }
    if (length(e) == 0) {
      return(NULL)
    }

    ## Shift the edge calls back as far as we padded ix from its start
    ## and up into the correct region of x
    e <- shift(e, istart.pad - istart + (istart - 1L))
    values(e)$multi.modal <- rep(length(e) > 1, length(e))
    e
  })

  if (nbad > 0) {
    message("couldn't locate peaks in ", nbad, " islands\n")
  }
  browser()
  edges <- do.call(c, unname(edges[!sapply(edges, is.null)]))
  if (!is.null(edges) && length(edges) > 0) {
    edges <- edges[width(edges) >= min.width]
  }
  edges
})

## Checks to see that bandwidth and window.size vectors are kosher and
## returuns them in decreasing order
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

## TODO: We are currently assuming the edge calling in local "ilsands"
##       so if the start,edge order looks wonky, we just take min/max
##       Make this a smarter "refining" of mixed start/edge pairs!
##       Make take the steeper of the two edges (or something)

##' Refines the edge calls in pre-specified regions of x.
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

