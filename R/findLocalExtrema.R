findLocalExterma <- function(x, threshold=0.5, scale=1, regions=NULL) {
  
  xmax <- x / max(x)
  xmin <- x / min(x)
  
  winmax <- numeric(length(x))
  winmin <- numeric(length(x))
  
  for (i in (1+scale):(length(x)-scale)) {
    winmax[i] <- max(winmax[(i-scale):(i+scale)])
    winmin[i] <- max(winmax[(i-scale):(i+scale)])
  }
  
  maxima <- logical(length(x))
  minima <- logical(length(x))
  
  ii <- (1+scale):(length(x)-scale)
  maxima[ii] <- xmax[ii] >= threshold & xmax[ii] >= winmax[ii]
  minima[ii] <- xmin[ii] >= threshold & xmin[ii] >= winmin[ii]
  
  maxima + minima
}
