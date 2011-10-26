setGeneric("detectPeaksByEdges",
function(x, bandwidth=40, mu=0, sd=1, threshold=0.15, edge.window=2*bandwidth,
         min.width=10, do.clean=TRUE, bandwidths=ceiling(c(.8, .5) * bandwidth),
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
function(x, bandwidth, mu, sd, threshold, edge.window, min.width, do.clean,
         bandwidths, .strand='+', min.height=0,
         ignore.from.start=NULL, ignore.from.end=NULL, ...) {
  if (min.height > 0) {
    peaks <- detectPeaksByEdges(Rle(x), bandwidth, mu, sd, threshold,
                                edge.window, min.width, do.clean, bandwidths,
                                .strand=.strand, min.height=min.height, ...)
    return(peaks)
  }

  edges <- detectEdges(x, bandwidth=bandwidth, mu=mu, sd=sd,
                       threshold=threshold, edge.window=edge.window, ...)

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
      edges$start <- edges$start[edges$start <= length(x) - ignore.from.end]
    }
    if (length(edges$end) > 0) {
      edges$end <- edges$end[edges$end <= length(x) - ignore.from.end]
    }
  }

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

  fishy <- .bad.edge.detection(edges$start, edges$end) ||
    !.start.end.inorder(edges$start, edges$end)

  if (do.clean && fishy) {
    edges <- refineEdges(x, bandwidth, edges$start, edges$end,
                         bandwidths=bandwidths, .strand=.strand, ...)
  }

  if (.bad.edge.detection(edges$start, edges$end) ||
      !.start.end.inorder(edges$start, edges$end)) {
    return(NULL)
  }

  ir <- IRanges(edges$start, edges$end)
  values(ir) <- DataFrame(fishy=rep(fishy, length(ir)))
  if (min.width > 0) {
    ir <- ir[width(ir) > min.width]
  }
  ir
})

setMethod("detectPeaksByEdges", c(x="Rle"),
function(x, bandwidth, mu, sd, threshold, edge.window, min.width, do.clean,
         bandwidths, pad.by=10, min.height=10L, .strand='+', ...) {
  f <- getMethod('detectPeaksByEdges', 'numeric')
  islands <- slice(x, lower=min.height, rangesOnly=TRUE,
                   includeLower=min.height != 0)
  islands <- islands[width(islands) >= min.width]

  nbad <- 0L

  edges <- lapply(1:length(islands), function(i) {
    istart <- start(islands[i])
    iend <- end(islands[i])
    istart.pad <- max((istart-pad.by), 1L)
    iend.pad <- min(iend+pad.by, length(x))
    ix <- as.numeric(x[istart.pad:iend.pad])

    e <- f(ix, bandwidth, mu, sd, threshold, edge.window, do.clean=do.clean,
           bandwidths=bandwidths, .strand=.strand, min.height=0,
           ignore.from.start=istart - istart.pad,
           ignore.from.end=iend.pad - iend, ...)

    if (is.null(e)) {
      nbad <<- nbad + 1L
      return(NULL)
    }
    if (length(e) == 0) {
      return(NULL)
    }

    values(e)$multi.modal <- length(e) > 1
    starts <- pmax(1L, start(e) - (istart - istart.pad))
    start(e) <- starts

    ends <- pmin(length(ix) - (iend.pad - iend), end(e) - pad.by)
    ## ends <- pmax(ends, starts + 1) ## edge was found in front padding!
    end(e) <- ends
    shift(e, istart)
  })

  if (nbad > 0) {
    message("couldn't locate peaks in ", nbad, " islands\n")
  }

  edges <- do.call(c, unname(edges))
  if (!is.null(edges) && length(edges) > 0) {
    edges <- edges[width(edges) >= min.width]
  }
  edges
})

refineEdges <- function(x, bandwidth, starts, ends,
                        bandwidths=ceiling(c(.8, .5) * bandwidth), .strand='+',
                        max.iter=10, detangle.by=c('max', 'min'), ...) {
  bandwidths <- sort(bandwidths, decreasing=TRUE)
  detangle.by <- match.arg(detangle.by)

  ## If there are no fenceposts, try to shrink bandwith
  while ((length(starts) == 0 || length(ends) == 0) && length(bandwidths)) {
    bandwidth <- ceiling(bandwidths[1])
    ## cat("trimming fat w/ bandwidth:", bandwidth, "\n")
    de <- detectEdges(x, bandwidth=bandwidth)
    ends <- de$end
    starts <- de$start
    bandwidths <- bandwidths[-1]
  }

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

  return(list(start=starts, end=ends))
  ## TODO: We are currently assuming the edge calling in local "lisands"
  ##       so if the start,edge order looks wonky, we just take min/max
  ##       Make this a smarter "refining" of mixed start/edge pairs!
  ##       Make take the steeper of the two edges (or something)
  ## df <- data.frame(pos=all.bounds,
  ##                  start=c(rep(TRUE, length(starts)), rep(FALSE, length(ends))))
  ## df <- df[order(df$pos, decreasing=.strand == '+'), ]
  ## r <- rle(df$start)
  ## i <- 1
  ## axe <- which(r$lengths > 1)

  ## ## If there are two starts in a row, only keep the first one
  ## ## If there are two ends in a row, only keep the first one also!
  ## while(length(axe) > 0 && i <= max.iter) {
  ##   if (r[axe]) {
  ##     ## This is the start of a peak

  ##   } else {
  ##     ## This is the end of a peak
  ##   }
  ##   rm.idx <- max(1, axe[1L] - 1)
  ##   df <- df[-sum(r$lengths[1:rm.idx]),]
  ##   r <- rle(df$start)
  ##   i <- i + 1
  ##   r <- rle(df$start)
  ##   axe <- which(r$lengths > 1)
  ## }


  ## if (.strand == '+') {
  ##   df <- df[rev(seq(nrow(df))),]
  ## }

  ## ret <- list(start=df$pos[df$start], end=df$pos[!df$start])
  ## if (length(ret$end) == 0 || length(ret$start) == 0) {
  ##   stop("0-length edges should not have happend w/ the quantile trick!")
  ## }

  ## ret
}

