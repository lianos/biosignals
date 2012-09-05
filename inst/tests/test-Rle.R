context("Rle")

test_that("Rle construction works", {
  x1 <- c(1, 1, 1, 1, 2, 2, 3, 4, 4, 4, 5, 6, 8, 8, 8.01)
  expect_true(all(biosignals:::asRle(x1) == Rle(x1)))
})

test_that("Rle constructor handles weird input", {
  ## An Rle with uniform input killed my entire sunday
  r1 <- biosignals:::asRle(rep(10, 5))
  ir1 <- Rle(rep(10, 5))
  expect_true(all(r1 == ir1))
  
  r2 <- biosignals:::asRle(numeric())
  ir2 <- Rle(numeric())
  expect_true(all(r2 == ir2))
})

test_that("Rle expansion works", {
  x <- c(0, 0, 1, 1, 1, 0.1, 0.1, 0.2, 3, 3, 4)
  r <- Rle(x)
  expect_equal(biosignals:::expandRle(r), x)
  expect_equal(biosignals:::expandRleS4(r), x)
})

test_that("Convolution over sparse Rle is bueno", {
  cvr <- readRDS(system.file('extdata', 'coverage.rds', package="biosignals"))
  all.islands <- slice(cvr, lower=0, includeLower=FALSE, rangesOnly=TRUE)
  
  ## these look like "normal" signals
  smooth.idx <- c(801, 9021, 9022)
  normal.islands <- all.islands[smooth.idx]
  normal.smooth <- convolve1d(Views(cvr, normal.islands))
  cvrs <- subject(normal.smooth)
  
  expect_equal(length(cvrs), length(cvr))
  expect_equal(length(normal.smooth), length(normal.islands))
  
  ## expect that a range that wasn't convolved isn't different,
  not.smooth <- setdiff(1:length(all.islands), smooth.idx)
  ns.idx <- sample(not.smooth)[1]
  expect_equal(cvrs[all.islands[100]], cvr[all.islands[100]])
  
  ## data in range that was convolved should be different
  s.idx <- sample(smooth.idx)[1]
  expect_false(all(cvrs[all.islands[801]] == cvr[all.islands[801]]))
})