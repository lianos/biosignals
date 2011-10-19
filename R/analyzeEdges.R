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
                         tresholds=c(0.1, 0.2, 0.3, 0.4),
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
