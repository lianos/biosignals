context("turnpoints")

test_that("sliding max/mins use +/- windows", {
  vals  <- c(-1, -1,  1,  1, 3, 5, 5, 5, 3, 3,  3,  10)
  xmax  <- c( 1,  1,  3,  5, 5, 5, 5, 5, 5, 10, 10, 10)
  xmin  <- c(-1, -1, -1, -1, 1, 1, 3, 3, 3, 3,  3,  3)
  
  expect_equal(slidingMax(vals, 2), xmax, info="sliding max")
  expect_equal(slidingMin(vals, 2), xmin, info="sliding max")
})