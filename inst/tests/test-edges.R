context("edge detection")

## This example is taken from the supplement of the FACADE paper:
## http://nar.oxfordjournals.org/content/38/15/e157.full
test_that("idealized step edge detected", {
  ## Idealized step
  x <- c(rep(0, 15), rep(20, 15), -1, rep(0, 20))
  x <- x + rnorm(length(x)) ## with some noise

  edges <- detectEdges(x)
  
  ## Give some slack on what the start/end calls are given the randomness
  ## in the noise added to x, the real start/end should be 16 and 31
  expect_true(edges$start %in% c(15:17))
  expect_true(edges$end %in% c(30:32))
})

