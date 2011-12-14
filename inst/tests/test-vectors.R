context("vector operations")

slidingExtrema <- function(x, k, extreme=c('min', 'max')) {
  extreme <- match.fun(match.arg(extreme))
  ans <- numeric(length(x))
  
  for (i in 1:length(x)) {
    idxs <- (i - k):(i + k)
    idxs <- idxs[idxs > 0 & idxs <= length(x)]
    ans[i] <- extreme(x[idxs])
  }
  
  ans
}

test_that("sliding min/max functions work", {
  x <- sample(1:50)
  expect_equal(slidingMax(x, 5), slidingExtrema(x, 5, 'max'))
  expect_equal(slidingMin(x, 5), slidingExtrema(x, 5, 'min'))
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
