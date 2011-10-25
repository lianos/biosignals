context("vector operations")

test_that("terminal quantile calculation is sane", {
  x <- c(6, 6, 6, 6, 6, 6, 7, 8, 9, 9, 10, 10, 9, 9, 10, 13, 13,
         18, 18, 18, 19, 25)
  qtiles <- quantilePositions(x)

  ## TODO: FIX THIS BUG
  ## 2011-10-25: Was getting a weird edge effect, see quantile.95:
  ##
  ##      q.05 q.1 q.15 q.2 q.8 q.85 q.9 q.9 q.95
  ##         1   3    5   7   8   20  21  22    1
  expect_false(is.unsorted(as.integer(qtiles)))
})
