##' Detects edges across a raw signal vector
##'
##' The gaussian smoothed first and second derivative of the data are used to
##' locate the edges in the non-negative vector \code{x}. Potential edges occur
##' at the peaks (and troughs) of the first derivative. These are located where
##' the second derivative crosses zero.
##'
##' The \code{edge.window} and code{threshold} arguments are used to filter out
##' noisy peaks.
##'
##' @param x A numeric (coverage) vector. This method assumes x is non-negative.
##' @param bandwidth The bandwidth for the kernel.
##' @param mu The mean of the kernel
##' @param sd The std deviation of the kernel
##' @param threshold A number between 0 and 1 used to filter out noise. It
##' indicates how obvious a step has to be in order to be considered a
##' transition. (not used now)
##' @param half.window The region around each point in \code{x} to use when
##' trying to filter out noise by step-size (not used now)
detectEdges <- function(x, bandwidth=10, mu=0, sd=1, threshold=0.15,
                        half.window=ceiling(bandwidth * 2), ...) {
  ## NOTE: The original noise-cleaning code was moved out and stored in
  ##       the `with-original-noise-fiter` branch. cf. `localExtrema.R` to get
  ##       an idea about the direction to go in.
  ## k1 <- generateKernel('gaussian', bandwidth=bandwidth, mu=mu, sd=sd, deriv=1)
  k2 <- generateKernel('gaussian', bandwidth=bandwidth, mu=mu, sd=sd, deriv=2)
  
  ## First derivatve not necessary if no noise-filtering is done
  ## x1 <- convolve1d(x, k1) ## first derivative of data
  x2 <- convolve1d(x, k2) ## locates peak of first deriv when this crosses 0

  ## Locations where second derivative cross 0 are the potential regions
  ## where an edge will be.
  zeroes <- zeroCrossings(x2)

  ## ---------------------------------------------------------------------------
  ## Find edges which might start a peak -- peaks of first derivative
  ##
  ## positions where the 2nd derivative switches from positive to negative
  ## indicate the potential start of an "upswing" (peak)
  edges.start <- zeroes[as.logical(x2[zeroes - 1] > 0)]
  ## Code to filter noise was removed from here.
  
  ## ---------------------------------------------------------------------------
  ## Find edges which might end a peak
  ##
  ## positions where 2nd derivative switches from negative to positive
  ## indicate potential end of a peak -- this is where an edge is
  edges.end <- zeroes[as.logical(x2[zeroes - 1] < 0)]
  ## Code to filter noise was removd here.
  
  list(start=edges.start, end=edges.end)
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
                           nt.y=0.03 * par('usr')[3:4], .strand='+',
                           draw.derivs=c(1, 2), ...) {
  args <- list(...)

  if (missing(x)) {
    ## this is a test
    x <- c(rep(0, 15), rep(20, 15), -1, rep(0, 20))
    x <- x + rnorm(length(x))
  }

  x <- as.numeric(x)

  k0 <- generateKernel('gaus', bandwidth=bandwidth, mu=mu, sd=sd, deriv=0)
  k1 <- generateKernel('gaus', bandwidth=bandwidth, mu=mu, sd=sd, deriv=1)
  k2 <- generateKernel('gaus', bandwidth=bandwidth, mu=mu, sd=sd, deriv=2)

  x0 <- convolve1d(x, k0)
  x1 <- convolve1d(x, k1)
  x2 <- convolve1d(x, k2)

  if (is.null(args$ylim)) {
    ylim <- c(min(0, x, x0, x1, x2), max(x, x0, x1, x2))
  } else {
    ylim <- args$ylim
  }

  plot(x, type='l', ylim=ylim, col='black', lwd=2)
  abline(h=0, lwd=3)

  draw.smooth <- args$draw.smooth
  if (is.null(draw.smooth)) {
    draw.smooth <- TRUE
  }
  if (draw.smooth) {
    lines(x0, col='black', lwd=1.5, lty='dotted')
  }

  if (is.numeric(draw.derivs)) {
    legend.text <- character()
    legend.col <- character()
    if (1 %in% draw.derivs) {
      lines(x1, col='red', lwd=1.5, lty='dotted')
      legend.text <- '1st deriv'
      legend.col <- 'red'
    }
    if (2 %in% draw.derivs) {
      lines(x2, col='blue', lwd=1.5, lty='dotted')
      legend.text <- c(legend.text, '2nd deriv')
      legend.col <- c(legend.col, 'blue')
    }
    legend('bottomright', legend=legend.text, text.col=legend.col)
  }

  edges <- detectPeaksByEdges(x, bandwidth=bandwidth, mu=mu, sd=sd, ...)

  draw.edges <- args$draw.edges
  if (is.null(draw.edges)) {
    draw.edges <- TRUE
  }
  if (draw.edges) {
    abline(v=start(edges), col='green', lty='dashed', lwd=2)
    abline(v=end(edges), col='red', lty='dashed', lwd=2)
    tryCatch({
      rect(start(edges), ylim[1], end(edges), ylim[2], col='#6f6f6f66',
           border=NA, density=-1)
    }, error=function(e) NULL)
  }

  if (is(nt, 'DNAString') || is.character(nt)) {
    nt.y[1] <- min(-2, nt.y[1])
    nt.y[2] <- max(2, nt.y[2])
    if (nt.y[2] > 0) {
      nt.y <- nt.y - (nt.y[2] + .25)
    }
    acgt <- c(A='green', C='blue', G='darkorange', T='red')
    if (.strand == '-') {
      nt <- reverseComplement(DNAString(nt))
    }
    nt <- as.character(nt)
    nt <- unlist(strsplit(nt, ''))
    cols <- acgt[nt]

    rect(1:length(x) - .45, nt.y[1], 1:length(x) + .45, nt.y[2],
         border=NA, col=cols)
  }

  invisible(edges)
}


if (FALSE) {
  library(biosignals)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg19)
  cvr <- readRDS("/Users/stavros/cBio/projects/biosignals/biosignals-pkg/inst/extdata/coverage.rds")
  bchr <- unmasked(Hsapiens$chr21)
  ## cvr <- load.it('')
  ## islands <- slice(x, lower=1, rangesOnly=TRUE)

  ##############################################################################
  ## Detect peaks by smoothing individual coverage islands
  islands <- slice(cvr, lower=10, rangesOnly=TRUE)
  maxs <- viewMaxs(Views(cvr, islands))
  check <- which(maxs > 100)
  bchr <- unmasked(Hsapiens$chr21)

  edges <- lapply(1:length(islands), function(i) {
    istart <- start(islands[i])
    iend <- end(islands[i])
    e <- detectPeaksByEdges(as.numeric(cvr[istart:iend]), bandwidth=40,
                            .strand='+')
    shift(e, istart)
  })
  edges <- do.call(c, unname(edges))
  edges <- cbind(data.frame(start=start(edges), end=end(edges)),
                 as.data.frame(values(edges)))
  for (i in 1:30) {
    istart <- start(islands[i]) - 10
    iend <- end(islands[i]) + 10
    cat(sprintf("%d : chr21:%d-%d\n", i + 1, istart, iend))
    quartz()
    de <- visualizeEdges(cvr[istart:iend], bandwidth=40,
                         nt=bchr[(istart):(iend)])
  }

  ##############################################################################
  ## As a function now
  e <- detectPeaksByEdges(cvr, bandwidth=40, min.height=10, min.width=21)
  edges <- cbind(data.frame(start=start(e), end=end(e)),
                 as.data.frame(values(e)))
  edges$width <- edges$end - edges$start
  edges$seqnames <- 'chr21'
  edges$strand <- '+'
  eg <- as(edges, 'GRanges')
  ## compare w/ it the old one
  info <- specializedSignalHumps(cvr.smooth, '+', eps=1e-6, lower=10,
                                 min.width=21, quantile.breaks=quantile.breaks)
  info <- data.frame(start=start(info), end=end(info),
                     is.distal=values(info)$is.distal,
                     fence.start=values(info)$quantile.02,
                     fence.end=values(info)$quantile.98)
  info$seqnames <- 'chr21'
  info$strand <- '+'
  ig <- as(info, 'GRanges')

  in.edges <- eg[!eg %in% ig]
  in.peaks <- ig[!ig %in% eg]

  for (i in sample(1:length(in.peaks), 30)) {
    istart <- start(in.peaks)[i] - 10
    iend <- end(in.peaks[i]) + 10
    quartz()
    visualizeEdges(cvr[istart:iend], bandwidth=40, nt=bchr[istart:iend])
    title(sprintf('chr21:%d-%d', istart, iend))
  }

  wtf <- GRanges('chr21', IRanges(38887641, 38887691), '+')
  ##############################################################################
  ## Find edges from the entire signal, then look for info within islands?
  cvr.x <- as.numeric(cvr)
  all.edges <- detectEdges(cvr.x, bandwidth=40, .strand='+')


}
