source("adam.R")
# These are tests. You can run them with testthat::test_file("test-adam_series.R").
# These tests are in another file because they may be slow,
# and the test harness lets you run one file at a time, not one test at a time.

relerr <- function(obs, expected) {(obs - expected) / expected}

test_that("adam series is inverse of pr2ar", {
  params <- makeParameters_Adam()
  time_periods <- 5
  ar <- rep(0.4, time_periods)
  pr <- ar2pr_series_Adam(ar, adam_parameters)
  ar_trajectory <- pr2ar_Adam(pr, adam_parameters)
  infer_ar <- ar_trajectory$alpha
  expect_equal(length(infer_ar) + 1, length(ar))
  expect_lt(abs(relerr(infer_ar[length(infer_ar)], ar[length(ar)])), 0.02)
})
