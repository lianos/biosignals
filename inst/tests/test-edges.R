context("edge detection")

## This example is taken from the supplement of the FACADE paper:
## http://nar.oxfordjournals.org/content/38/15/e157.full
test_that("idealized step edge detected", {
  ## Idealized step
  x <- c(rep(0, 15), rep(20, 15), -1, rep(0, 20))
  x <- x + rnorm(length(x)) ## with some noise

  edges <- detectEdges(x)
  expect_equal(edges$start, 16)
  expect_equal(edges$end, 31)
})

