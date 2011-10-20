setGeneric("convolve1d", function(x, ...) standardGeneric("convolve1d"))

setMethod("convolve1d", c(x="numeric"),
function(x, kernel='normal', rescale=TRUE, bandwidth=20, deriv=0,
         starts=NULL, ends=NULL, ...) {
  if (is.character(kernel)) {
    kernel <- generateKernel(kernel, bandwidth=bandwidth, ...)
  }
  stopifnot(is.numeric(kernel))
  stopifnot(is.logical(rescale) && length(rescale) == 1L)
  
  kernel <- rev(kernel)
  
  if (!is.null(starts) || !is.null(ends)) {
    if (is.null(starts) || is.null(ends)) {
      stop("both starts and ends need to be defined")
    }
    if (length(starts) != length(ends)) {
      stop("length of starts and ends must be equal")
    }
    if (any(starts >= ends)) {
      stop("all starts must be less than ends")
    }
    if (min(starts) < 1) {
      stop("start positions must be >= 1")
    }
    if (max(ends) > length(x)) {
      stop("end positions must be <= length(x)")
    }
    starts <- as.integer(starts)
    ends <- as.integer(ends)
    
    ret <- .Call("Rfencepost_convolve_1d", x, kernel, starts, ends, rescale,
                 PACKAGE="biosignals")
  } else {
    ret <- .Call("Rconvolve_1d", x, kernel, rescale, PACKAGE="biosignals")
  }
  
  ret
})

setMethod("convolve1d", c(x="Rle"),
function(x, kernel='normal', rescale=TRUE, bandwidth=20, deriv=0, lower=0,
         ...) {
  islands <- slice(x, lower=lower, includeLower=lower != 0, rangesOnly=TRUE)
  ret <- convolve1d(as.numeric(x), kernel, rescale=rescale,
                    bandwidth=bandwidth, starts=start(islands),
                    ends=end(islands), ...)
  Rle(ret)
})
