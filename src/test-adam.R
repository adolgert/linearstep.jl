source("adam.R")
# These are tests. You can run them with testthat::test_file("test-adam.R").

abserr <- function(observed, expected) {
  abs((observed - expected) / expected)
}


test_that("step-up distributions have lengths 8 and 9", {
  r1 <- xiF.x(0.1)
  expect_equal(length(r1), 8)
  expect_lt(abserr(sum(r1), 1), 1e-7)

  r2 <- xiF_0(0.1)
  expect_equal(length(r2), 8)
  expect_lt(abserr(sum(r2), 1), 1e-7)

  r3 <- xiF_1(0.1)
  expect_equal(length(r3), 9)
  expect_lt(abserr(sum(r3), 1), 1e-7)

  r4 <- xiF.2(c(0.1, 0.2))
  expect_equal(length(r4), 8)
  expect_lt(abserr(sum(r4), 1), 1e-7)

  r5 <- xiF.3(0.1)
  expect_equal(length(r5), 9)
  expect_lt(abserr(sum(r5), 1), 1e-7)

  r6 <- xiF_h(0.1)
  expect_equal(length(r6), 8)
  expect_lt(abserr(sum(r6), 1), 1e-7)

  r7 <- xiF_1(0.9)
  expect_lt(abserr(sum(r3), 1), 1e-7)

})


test_that("stage up superinfection for adam model", {
  state_cnt <- 9
  sus_matrix <- SuSblock_Adam(list(N = state_cnt))
  expect_equal(dim(sus_matrix)[1], state_cnt + 1)
  expect_equal(dim(sus_matrix)[2], state_cnt + 1)
})


test_that("pfpr by age runs", {
  skip("this makes plots, so skip it.")
  params <- makeParameters_Adam()
  pfprXage(.1, 0, params) -> cXX0
  pfprXage(.1, .3, params) -> cXX1
  pfprXage(.5, .1, params) -> cXX2
  pfprXage(.5, .3, params) -> cXX3
})



test_that("pr can become ar for a series", {
  params <- makeParameters_Adam()
  ar <- c(0.3, 0.31, 0.34, 0.35, 0.32)
  pr <- ar2pr_series_Adam(ar, params)
  expect_equal(length(ar), length(pr))
  expect_true(all(is.finite(pr)))
})


test_that("cohort creates stable wave", {
  params <- makeParameters_Adam()
  cXX = cohortXX_Adam(.1,params)
  XX = cohort2ages_Adam(cXX,params)
  XX = ar2stableWave_Adam(.1, params)
})


test_that("burden runs on attack rates", {
  params <- makeParameters_Adam()
  ar <- c(0.1, 0.12, 0.13, 0.12, 0.11, 0.13)
  burden <- ar2burden_Adam(ar, params)
  expect_equal(nrow(burden), length(ar))
})


test_that("integrated foi works for constant", {
    h <- 0.001
    step_days <- 10
    agg_days <- 365
    step_cnt <- 37 * 2  # 2+ years
    ar <- rep(1 - exp(-h * step_days), step_cnt)
    ar <- c(ar, rep(1 - exp(-0.015 * step_days), step_cnt))
    foi <- aggregate_foi(ar, agg_days, step_days)
    expect_equal(length(foi), 4)
    expect_lt(abs(foi[1] - h * 365), 1e-7)
    expect_lt(abs(foi[length(foi)] - .015 * 365), 1e-7)
})


test_that("can do 20 steps", {
  par <- makeParameters_Adam()
  tt = 1:20
  alphat = c(0.1, .1 + .01*(1 + sin(2*pi*tt/36.5)))
  XX0 = ar2stableWave_Adam(alphat[1], par)
  XX = XX0
  xx = XX2pr29_Adam(XX0,par)
  for (i in 2:length(alphat)) {
    XX = PtimesX_Adam(alphat[i], XX, par)
    xx = c(xx,XX2pr29_Adam(XX,par))
  }
  alpha1 = pr2ar_Adam(xx, par, alphat[1],XX0)
  alpha2 = pr2ar_Adam(xx, par)
})
