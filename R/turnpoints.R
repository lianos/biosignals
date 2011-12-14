turnpoints <- function(x, threshold=0.5, win.length=1L,
                       starts=NULL, ends=NULL) {
  x.rmax <- x / max(x)
  x.rmin <- x / min(x) ## this only works when min(x) != negative

  winmax <- numeric(length(x))
  winmin <- numeric(length(x))

  global.idxs <- (1L + win.length):(length(x) - win.length)
  for (i in global.idxs) {
    idxs <- (i-win.length):(i+win.length)
    winmax[i] <- max(x.rmax[idxs])
    winmin[i] <- min(x.rmin[idxs])
  }

  maxima <- numeric(length(x))
  minima <- numeric(length(x))

  ii <- global.idxs
  maxima[ii] <- x.rmax[ii] >= threshold & x.rmax[ii] >= winmax[ii]
  minima[ii] <- x.rmin[ii] >= threshold & x.rmin[ii] >= winmin[ii]

  list(maxima=maxima, minima=minima)
}
