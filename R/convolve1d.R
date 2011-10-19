convolve1d <- function(x, kernel, bandwidth=20, rescale=TRUE) {
  # x <- c(rep(x[1], length(kernel)), x, rep(tail(x, 1), length(kernel)))
  ret <- .Call("Rconvolve_1d", as.numeric(x), as.numeric(kernel), rescale,
               PACKAGE="biosignals")
  ret
}

convolve1d.Rle <- function(x, kernel, bandwidth=20, rescale=TRUE,
                           lower=0) {
  xn <- as.numeric(x)
  islands <- slice(x, lower=lower, includeLower=lower != 0, rangesOnly=TRUE)
  starts <- start(islands)
  ends <- ends(islands)
  ret <- convolve1d.fp(x, kernel, starts, ends, bandwidth=bandwidth,
                       rescale=rescale)
}

convolve1d.fp <- function(x, kernel, starts, ends, bandwidth=20,
                          rescale=TRUE) {
  ret <- .Call("Rfencepost_convolve_1d", as.numeric(x), as.numeric(kernel),
               as.integer(starts), as.integer(ends), rescale)
}
