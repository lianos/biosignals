##' @nord
##' The types of kernel we know how to use
.convolutionKernels <- function() {
  c('uniform', 'rectangular', 'gaussian', 'normal', 'triangular',
    'epanechnikov', 'quartic', 'biweight', 'tricube', 'triweight',
    'cosine', 'tukey', 'laplacian', 'lofg', 'dog')
}

##' Generates the kernel window function
##'
##' @author Anshul Kundaje (original MATLAB implementation)
##' @author Steve Lianoglou (ported to R)
##'
##' @param kernel.type Name of kernel to use
##' @param bandwidth kernel half width
##' @param normalize \code{logical(1)}: normalizes kernel if (\code{TRUE}).
##'
##' @return a vector with kernel values
generateKernel <- function(kernel.type='normal', bandwidth=1L, normalize=TRUE,
                           ...) {
  kernel.type <- match.arg(kernel.type, .convolutionKernels())
  if (bandwidth <= 0) {
    stop("kernel bandwidth must be > 0")
  }

  args <- list(...)
  mu <- if (!is.null(args$mu)) args$mu else 0
  sd <- if (!is.null(args$sd)) args$sd else 1

  ## lofg = laplacian of gaussian == first deriv of gaus
  ## dog = difference of gaussians
  is.gauss.like <- c('normal', 'gaussian', 'lofg', 'dog')
  if (kernel.type %in% is.gauss.like) {
    ## build a kernel values to calculate that are centered around the mean
    ## of the distribution
    if (missing(bandwidth)) {
      bandwidth <- ceiling(3 * sd)
    }
    win.len <- ceiling(2L * bandwidth + 1L)
    # center <- ceiling(win.len / 2)
    kin <- seq(-3*sd, 3*sd, length.out=win.len) + mu
  } else {
    win.len <- 2L * ceiling(bandwidth) + 1L
    kin <- seq(-1, 1, length.out=win.len)
  }

  if (kernel.type %in% c('normal', 'gaussian')) {
    kout <- dnorm(kin, mu, sd)
    if (is.numeric(args$deriv) && args$deriv > 0) {
      deriv <- args$deriv
      stopifnot(deriv < 5)
      # center <- ceiling(length(kval) / 2) + mu
      # krange <- seq_along(kval)
      derivs <- matrix(1, nrow=4, ncol=length(kin))
      derivs[1,] <- -(kin / (sd^2))
      derivs[2,] <- (kin^2 - (sd^2)) / sd^4
      ## derivs[2,] <- kin^2/sd^4 - sd^2
      derivs[3,] <- -(kin^3 - 3*sd^2*kin) / sd^6
      derivs[4,] <- (kin^4 - 6*sd^2*kin^2 + 3*sd^4) / sd^8
      kout <- kout * derivs[deriv,]
    }
  } else if (kernel.type %in% c('unform', 'rectangular')) {
    kout <- rep(1, win.len)
  } else if (kernel.type == 'triangular') {
    kout <- 1 - abs(kin)
  } else if (kernel.type == 'epanechnikov') {
    kout <- (3/4) * (1 - kin^2)
  } else if (kernel.type %in% c('quartic', 'biweight')) {
    kout <- (15 / 16) * (1 - kin^2)^2
  } else if (kernel.type %in% c('triweight', 'tricube')) {
    kout <- (35/32) * (1 - kin^2)^3
  } else if (kernel.type == 'cosine') {
    kout <- (pi / 4) * cos(pi * kin / 2)
  } else if (kernel.type == 'laplacian') {
    kout <- c(1, -2, 1)
    ## return(kval)
  } else if (kernel.type == 'lofg') {
    ## laplacian of gaussian (first deriv)
    ## http://www.ce.rit.edu/~cockburn/courses/cvf05/Hw1_corrected.pdf
    ## \frac{d}{dx}G = - \frac{2}{\sqrt{2\pi}\sigma^3} \times
    ## \exp (- \frac{x^2}{\sigma^2} )
    kout <- (2 / (sqrt(2*pi) * sd^3)) * kin * exp(-1 * kin^2 / sd^2)
  } else if (kernel.type == 'dog') {
    ## diference of gaussians
    ## http://en.wikipedia.org/wiki/Difference_of_Gaussians
    dog <- if (is.numeric(args$dog)) args$dog else 1.6
    kout <- dnorm(kin, mu, sd) - dnorm(kin, mu, sd * dog)
  } else if (kernel.type == 'tukey') {
    stop("wut's a tukey?")
  } else {
    stop("wut kernel is this: ", kernel.type)
  }

  ## if (normalize) {
  ##   kval <- kval / sum(kval)
  ## } else {
  ##   kval <- kval / max(kval)
  ## }

  kout
}

# dog <- function(x=seq(-4, 4, length.out=20), mu=0, sd1=1, sd2=2) {
#   dnorm(x, mu, sd1) - dnorm(x, mu, sd2)
# }

g1 <- function(mu=0, sd=1, bandwidth=5) {
  if (missing(bandwidth)) {
    bandwidth <- ceiling(3 * sd)
  }
  win.len <- ceiling(2L * bandwidth + 1L)
  center <- ceiling(win.len / 2)
  x <- seq(win.len) - (center + mu)
  (-x / (sd^3 * sqrt(2*pi))) * exp(-1 * (x^2 / (2 * sd^2)))
}

g <- function(mu=0, sd=1, bandwidth=5) {
  if (missing(bandwidth)) {
    bandwidth <- ceiling(3 * sd)
  }
  win.len <- ceiling(2L * bandwidth + 1L)
  center <- ceiling(win.len / 2)
  x <- seq(win.len) - (center + mu)
  (1 / (sqrt(2*pi) * sd)) * exp(-(x^2/(2*sd^2)))
}
