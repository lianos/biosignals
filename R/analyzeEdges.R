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
                         thresholds=c(0.1, 0.2, 0.3, 0.4),
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

##' Detects edges across a signal vector
detectEdges <- function(x, bandwidth=10, mu=0, sd=1, threshold=0.15,
                        edge.window=2*bandwidth, thresh.mult=0.3, ...) {
  k1 <- generateKernel('gaussian', bandwidth=bandwidth, mu=mu, sd=sd, deriv=1)
  k2 <- generateKernel('gaussian', bandwidth=bandwidth, mu=mu, sd=sd, deriv=2)
  k3 <- generateKernel('gaussian', bandwidth=bandwidth, mu=mu, sd=sd, deriv=3)

  x1 <- convolve1d(x, k1) ## first derivative of data
  x2 <- convolve1d(x, k2) ## locates peak of first deriv when this crosses 0
  x3 <- convolve1d(x, k3) ## how steep is the edge?

  ## Find where second deriv crosses 0
  zeroes <- zeroCrossings(x2)

  ## Use x1 to filter out zero crossings that represent noise
  ## ax1 <- as.numeric(abs(x1))
  ## mx <- slidingMax(ax1, k=edge.window)
  ## potential <- which(ax1 >= mx * threshold)

  ## ax3 <- abs(x3)
  ## ## mx <- max(abs(x3))
  ## ## potential <- which(abs(x3) >= threshold * mx)
  ## mx <- slidingMax(ax3, k=edge.window)
  ## potential <- which(ax3 >= threshold * mx)

  ## edges <- intersect(zeroes, potential)
  ## edges.start <- edges[as.logical(x2[edges - 1] > 0)]
  ## edges.end <- edges[as.logical(x2[edges - 1] < 0)]

  ## Filter out edges by ensuring that the value at the peaks of x1 (deriv 1)
  ## "clear" the minima of their surround for start of edges, vice versa
  ## for ends of edges
  edges.start <- zeroes[as.logical(x2[zeroes - 1] > 0)]
  bad.starts <- as.logical(x1[edges.start] <= 0)
  if (any(bad.starts)) {
    edges.start <- edges.start[!bad.starts]
  }
  x1s <- x1
  x1s[x1s <= 1] <- 1 ## only look at sliding minima for values north of 0
  min.x1 <- slidingMin(x1s, bandwidth)
  ## min.x1[edges.start] / x1[edges.start]
  start.offsets <- (threshold/2) * as.numeric(x1[edges.start])
  es <- edges.start[(min.x1[edges.start] + start.offsets) /
                    (as.numeric(x1[edges.start]) + start.offsets) < threshold]

  ## browser()

  edges.end <- zeroes[as.logical(x2[zeroes - 1] < 0)]
  bad.ends <- as.logical(x1[edges.end] >= 0)
  if (any(bad.ends)) {
    edges.end <- edges.end[!bad.ends]
  }
  x1e <- x1
  x1e[x1e >= -1] <- -1 ## only look for sliding maxima for vals south of 0
  max.x1 <- slidingMax(x1e, bandwidth)
  ## tmn <- as.numeric(thresh.mult * x1[edges.end]) - 1 ## -1 in case x1[] == 0
  ## ee <- edges.end[(max.x1[edges.end] + tmn) / tmn < threshold]
  end.offsets <- threshold * as.numeric(x1[edges.end])
  ee <- edges.end[(max.x1[edges.end] + end.offsets) /
                  (as.numeric(x1[edges.end] + end.offsets)) < threshold]

  list(start=es, end=ee)
}

################################################################################
## Debugging
## -----------------------------------------------------------------------------

## To find "real" edges, look for zero-crossings of the gauss,deriv=2 kernel
## and filter this using "clear peaks" found in gauss,deriv=1 kernel
##
## A "clear peak" has to do with the y-value of the peak's center. This will
## be large -- but what is "large"? This is where you use the "local maximum"
## combined w/ threshold, somehow.
##
## peaks starts when second deriv goes from + to -
## peaks end when second deriv goes from - to +
visualizeEdges <- function(x, mu=0, sd=1, bandwidth=10, threshold=0.15,
                           edge.window=2*bandwidth, nt=NULL,
                           nt.y=0.03 * par('usr')[3:4], .strand='+') {
  if (missing(x)) {
    ## this is a test
    x <- c(rep(0, 15), rep(20, 15), -1, rep(0, 20))
    x <- x + rnorm(length(x))
  }

  x <- as.numeric(x)

  k0 <- generateKernel('gaus', bandwidth=bandwidth, mu=mu, sd=sd, deriv=0)
  k1 <- generateKernel('gaus', bandwidth=bandwidth, mu=mu, sd=sd, deriv=1)
  k2 <- generateKernel('gaus', bandwidth=bandwidth, mu=mu, sd=sd, deriv=2)
  k3 <- generateKernel('gaus', bandwidth=bandwidth, mu=mu, sd=sd, deriv=3)
  ## kd <- generateKernel('dog', bandwidth=bandwidth, mu=mu, sd=sd, dog=1.8)
  ## browser()
  x0 <- convolve1d(x, k0)
  x1 <- convolve1d(x, k1)
  x2 <- convolve1d(x, k2)
  x3 <- convolve1d(x, k3)
  ## xd <- convolve1d(x, kd)

  ylim <- c(min(x, x0, x1, x2), max(x, x0, x1, x2))

  ## browser()

  plot(x, type='l', ylim=ylim, col='blue', lwd=2)
  abline(h=0)
  lines(x0, col='blue', lty='dashed', lwd=2)
  lines(x1, col='red', lwd=2, lty='dashed')
  lines(x2, col='green', lwd=2, lty='dashed')
  lines(x3, col='orange', lwd=1, lty='dashed')


  legend('bottomright', legend=c('dervi1', 'deriv2'), text.col=c('red', 'green'))

  ## zz <- zeroCrossings(x2)
  ## mx <- slidingMax(abs(x1), k=bandwidth * 2)
  ## possibly <- which(abs(x1) >= mx * threshold)
  ## edges <- intersect(zz, possibly)
  edges <- detectEdges(x, bandwidth=bandwidth, mu=mu, sd=sd, threshold=threshold,
                       edge.window=edge.window)
  if (length(edges$start) != length(edges$end)) {
    cat("unbalanced ... cleaning ...\n")
    edges <- cleanBoundaries(edges$start, edges$end, length(x), .strand=.strand)
  }
  abline(v=edges$start, col='green', lty='dashed', lwd=2)
  abline(v=edges$end, col='red', lty='dashed', lwd=2)

  if (is(nt, 'DNAString') || is.character(nt)) {
    acgt <- c(A='green', C='blue', G='darkorange', T='red')
    if (.strand == '-') {
      nt <- reverseComplement(DNAString(nt))
    }
    nt <- as.character(nt)
    nt <- unlist(strsplit(nt, ''))
    cols <- acgt[nt]
    rect(1:length(x) - .45, nt.y[1], 1:length(x) + .45, nt.y[2], border=NA, col=cols)
  }

  invisible(edges)
}

cleanBoundaries <- function(starts, ends, len, .strand='+', max.iter=10) {
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
  if (length(ret$end) == 0) {
    ret$end <- len
  }
  if (length(df$start) == 1) {
    ret$start <- 1
  }
  ret
}

if (FALSE) {
  library(biosignals)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg19)
  cvr <- readRDS("/Users/stavros/cBio/projects/biosignals/biosignals-pkg/inst/extdata/coverage.rda")
  ## cvr <- load.it('')
  ## islands <- slice(x, lower=1, rangesOnly=TRUE)
  islands <- slice(cvr, lower=10, rangesOnly=TRUE)
  maxs <- viewMaxs(Views(cvr, islands))
  check <- which(maxs > 100)
  bchr <- unmasked(Hsapiens$chr21)

  edges <- lapply(1:length(islands), function(i) {
    istart <- start(islands[i]) - 10
    iend <- end(islands[i]) + 10
    tryCatch({
      e <- detectEdges(as.numeric(cvr[istart:iend]), bandwidth=40)
    }, error=function(e) cat("error:", i, "\n"))

    if (length(e$start) != length(e$end)) {
      cat("unmatched:", i, " ... decreasing bandwidth\n")
      e2 <- detectEdges(as.numeric(cvr[istart:iend]), bandwidth=20)
      if (length(e2$start) != length(e2$end)) {
        cat('... still screwed, just boxing it')
        e <- cleanBoundaries(e$start, e$end, len=iend - istart + 1)
      } else {
        e <- e2
      }


      if (length(e$start) != length(e$end)) {
        cat('... BAD clean\n')
      } else {
        cat('... cleaned\n')
      }
    }
    if (length(e$start) + length(e$end) != 2) {
      cat("more than 1:", i, "\n")
    }
    e
  })

  for (i in 1:30) {
    istart <- start(islands[i]) - 10
    iend <- end(islands[i]) + 10
    cat(sprintf("%d : chr21:%d-%d\n", i + 1, istart, iend))
    quartz()
    de <- visualizeEdges(cvr[istart:iend], bandwidth=50, nt=bchr[(istart):(iend)])
  }

}
