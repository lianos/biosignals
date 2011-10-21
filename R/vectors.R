slidingMax <- function(x, k=5L) {
  ret <- .Call('Rsliding_max', as.numeric(x), as.integer(k),
               PACKAGE="biosignals")
}

slidingMin <- function(x, k=5L) {
  .Call('Rsliding_min', as.numeric(x), as.integer(k), PACKAGE="biosignals")
}

zeroCrossings <- function(x) {
  .Call('Rzero_crossings', as.numeric(x), PACKAGE="biosignals")
}