context("vector operations")

test_that("Rle construction works", {
  x1 <- c(1, 1, 1, 1, 2, 2, 3, 4, 4, 4, 5, 6, 8, 8, 8.01)
  expect_true(all(biosignals:::asRle(x1) == Rle(x1)))
})

test_that("Rle expansion works", {
  x <- c(0, 0, 1, 1, 1, 0.1, 0.1, 0.2, 3, 3, 4)
  r <- Rle(x)
  expect_equal(biosignals:::expandRle(r), x)
})

test_that("terminal quantile calculation is sane", {
  quantile.breaks <- c(0.05, 0.10, 0.15, 0.20, 0.80, 0.85, 0.90, 0.95)
  x <- c(6, 6, 6, 6, 6, 6, 7, 8, 9, 9, 10, 10, 9, 9, 10, 13, 13,
         18, 18, 18, 19, 25)
  qtiles <- quantilePositions(x, quantile.breaks=quantile.breaks)


  ## 2011-10-25: Was getting a weird edge effect, see quantile.95:
  ##
  ##      q.05 q.1 q.15 q.2 q.8 q.85 q.9 q.9 q.95
  ##         1   3    5   7   8   20  21  22    1
  expect_false(is.unsorted(as.integer(qtiles)))

  ## This isn't exactly checking quantiles are correct, but
  ## close enough for now:
  ##   all values in x at qtile positions should be valid
  cuts <- cut(cumsum(x) / sum(x), quantile.breaks)
  qq <- as.integer(cuts)
  q.at <- which(c(TRUE, diff(qq) != 0) & !is.na(qq))
  expect_true(all(q.at %in% qtiles))
})
