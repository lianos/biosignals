setGeneric("detectPeaksByEdges",
function(x, bandwidth=40, mu=0, sd=1, threshold=0.15, edge.window=2*bandwidth,
         min.width=10, do.clean=TRUE, bandwidths=ceiling(c(.8, .5) * bandwidth),
         ...) {
  standardGeneric("detectPeaksByEdges")
})

.is.fishy.edge.detection <- function(starts, ends) {
  length(starts) == 0L || length(starts) != length(ends)
}

##' Returns NULL if there was an error
setMethod("detectPeaksByEdges", c(x="numeric"),
function(x, bandwidth, mu, sd, threshold, edge.window, min.width, do.clean,
         bandwidths, .strand='+', min.height=10, ...) {
  if (min.height > 0) {
    peaks <- detectPeaksByEdges(Rle(x), bandwidth, mu, sd, threshold,
                                edge.window, min.width, do.clean, bandwidths,
                                .strand=.strand, min.height=min.height, ...)
    return(peaks)
  }

  edges <- detectEdges(x, bandwidth=bandwidth, mu=mu, sd=sd,
                       threshold=threshold, edge.window=edge.window, ...)

  looks.fishy <- .is.fishy.edge.detection(edges$start, edges$end)

  if (do.clean && looks.fishy) {
    edges <- refineEdges(x, bandwidth, edges$start, edges$end,
                         bandwidths=bandwidths, .strand=.strand, ...)
  }

  if (.is.fishy.edge.detection(edges$start, edges$end)) {
    return(NULL)
  }

  ir <- IRanges(edges$start, edges$end)
  values(ir) <- DataFrame(fishy=rep(looks.fishy, length(ir)))
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
    ix <- as.numeric(cvr[(istart-pad.by):(iend+pad.by)])

    e <- f(ix, bandwidth, mu, sd, threshold, edge.window, do.clean=do.clean,
           bandwidths=bandwidths, .strand=.strand, min.height=0, ...)

    if (is.null(e)) {
      nbad <<- nbad + 1L
      return(NULL)
    }
    if (length(e) == 0) {
      return(NULL)
    }

    values(e)$multi.modal <- length(e) > 1
    starts <- pmax(1L, start(e) - pad.by)
    start(e) <- starts

    ends <- pmin(length(ix) - pad.by, end(e) - pad.by)
    ends <- pmax(ends, starts + 1) ## edge was found in front padding!
    end(e) <- ends
    shift(e, istart)
  })

  if (nbad > 0) {
    message("couldn't locate peaks in ", nbad, " islands\n")
  }

  edges <- do.call(c, unname(edges))
  edges <- edges[width(edges) >= min.width]
  edges
})

refineEdges <- function(x, bandwidth, starts, ends,
                        bandwidths=ceiling(c(.8, .5) * bandwidth), .strand='+',
                        max.iter=10, ...) {
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

  all.bounds <- c(starts, ends)
  o <- order(all.bounds, decreasing=.strand == '+')
  df <- data.frame(pos=all.bounds,
                   start=c(rep(TRUE, length(starts)), rep(FALSE, length(ends))))
  df <- df[o,]
  r <- rle(df$start)
  i <- 1

  while(length(axe <- which(r$lengths > 1)) && i <= max.iter) {
    df <- df[-sum(r$lengths[1:axe]),]
    r <- rle(df$start)
    i <- i + 1
  }


  if (.strand == '+') {
    df <- df[rev(seq(nrow(df))),]
  }

  ret <- list(start=df$pos[df$start], end=df$pos[!df$start])
  if (length(ret$end) == 0 || length(ret$start) == 0) {
    stop("0-length edges should not have happend w/ the quantile trick!")
  }

  ret
}

