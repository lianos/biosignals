##' Computes the gaussian scale space of a 1D dataset
##' 
##' Scale parameters are in data-spacing units
##' 
##' @param x A 1d data vector
##' @param deriv The order of gaussian scale space to compute.
##' 0 is a smoothing scale space (no derivatives);
##' 1 is an edge detecting space
##' 
##' @scales A vector of scales (sigma) to compute
##' 
##' @return A scale space representation of the data
createGaussScaleSpace <- function(x, deriv=0, scales=1:20) {
  stopifnot(deriv <= 5)
  m <- matrix(nrow=length(x), ncol=length(scales))
  width <- 3
  
  for (i in 1:length(scales)) {
    kernel <- .gaussianKernel1D(i, deriv)
    m[, i] <- kdeSmooth(x, kernel)
  }
  
  m
}

.gaussianKernel1D <- function(sigma, deriv=0, width=3) {
  stopifnot(deriv <= 5)
  ## generateKernel creates a window length using this function:
  ## win.len <- 2L * ceiling(bandwidth) + 1L
  ## 
  ## The MATLAB GuassianKernel1D calculates window width to be
  ## 2 * width(=3) * ceiling(scale) + 1 
  ## 
  ## We need bandwidth = width * ceiling(scale)
  ##sbandwidth <- width * ceiling(sigma)
  ##k <- generateKernel('normal', bandwidth=sbandwidth, sd=sigma)
  range <- 1:(2*width*ceiling(sigma) + 1)
  center <- ceiling(length(range) / 2)
  k <- (1 / (sigma * sqrt(2*pi))) * exp(-((range - center)^2) / (2*sigma^2))
  if (deriv > 0) {
    ## Calculate upto 5 derivs
    center <- ceiling(length(k) / 2)
    krange <- seq_along(k)
    derivs <- matrix(1, nrow=5, ncol=length(k))
    derivs[2,] <- -((krange - center) / (sigma^2))
    derivs[3,] <- ((krange - center)^2 - (sigma^2)) / sigma^4
    derivs[4,] <- -((krange - center)^3 - 3*sigma^2*(krange - center)) / sigma^6
    derivs[5,] <- ((krange - center)^4 - 6*sigma^2*(krange-center)^2 + 3*sigma^4) / sigma^8
    k <- k * derivs[deriv + 1L,]
  }
  k
}
