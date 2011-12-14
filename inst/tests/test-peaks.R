context("peak calling")

test_that("double starts (ends) are refined reasonably in a local peak", {
  x <- c(   1,    0,    0,    0,    0,    0,    0,    0,    0,    9,  778,  933,
          955,  984,  995,  998, 1002, 1003, 1003, 1004, 1004, 1005, 1006, 1008,
         1009, 1009, 1012, 1013, 1013, 1013, 1013, 1013, 1014, 1014, 1015, 1015,
         1015, 1015, 1020, 1044, 1059, 1061, 1094, 1097, 1102, 1136,  793,  975,
         1181, 1206, 1245, 1256, 1257, 1256, 1256, 1257, 1258, 1257, 1260, 1264,
         1264, 1264, 1263, 1262, 1263, 1264, 1264, 1264, 1263, 1263, 1263, 1263,
         1263, 1263, 1258, 1235, 1220, 1218, 1186, 1184, 1180, 1138,  728,  397,
          248,  570,  544,  547,  547,  550,  553,  552,  552,  552,  549,  546,
          547,  549,  548,  548,  548,  547,  547,  547,  547,  548,  547,  547,
          547,  547,  547,  546,  546,  546,  545,  544,  543,  542,  526,  520,
          441,   65,   41,   24,   19,   16,   13,   12,   11,   11,   10,    7,
            5,    3,    2,    2,    1,    1,    1,    1,    1)
  de <- detectPeaksByEdges(x, bandwidth=36, ignore.from.start=5)
  expect_is(de, 'IRanges')
  
  ## This test was for noise filtering, which has been removed. The 'no-filter'
  ## code will find 2 (UGLY) peaks here
  expect_equal(length(de), 2, info="check `ignore.from.start`")
})

test_that("signals do not cause IRanges with negative widths", {
  x <- c(18, 21, 20, 20, 19, 12, 13, 13,  9,  9, 10, 10, 10,
         10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11,
         11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 54, 52, 62,
         62, 63, 63, 63, 62, 62, 62, 62, 61, 61, 61, 61,
         61, 61, 61, 61, 61, 61, 61, 60, 60, 60, 60, 60, 60,
         60, 60, 60, 60, 60, 59, 59, 59, 15, 14,  1,  1,
         0 , 0 , 0 , 0 , 0 , 0 , 0 , 0)
  de <- detectPeaksByEdges(x, bandwidth=36, bandwidths=c(28, 21),
                           ignore.from.start=10, ignore.from.end=10)
  expect_is(de, 'IRanges')
  expect_equal(length(de), 1)
  expect_true(width(de)[1] > 35 && width(de)[1] < 42)
})

