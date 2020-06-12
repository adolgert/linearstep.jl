source("scale_ar.R")


test_that("can get work from districts", {
  dist_dict <- data.table(
    name = c("KOBOKO", "MARACHA", "OYAM", "GULU"),
    id = c(1, 2, 3, 4)
  )
  work <- task_work(dist_dict, 3)
  expect_equal(work$district_name, "OYAM")
  expect_equal(work$draw_id, 1)
  expect_equal(work$district, 3)

  work <- task_work(dist_dict, 3, median_only = TRUE)
  expect_equal(work$district_name, "OYAM")
  expect_equal(work$draw_id, 0)
  expect_equal(work$district, 3)

  work <- task_work(dist_dict, 6)
  expect_equal(work$district_name, "MARACHA")
  expect_equal(work$draw_id, 2)
  expect_equal(work$district, 2)
})


test_that("expanding case dates is normal", {
  year_range <- 2000:2019
  case_dt <- data.table(
    date = as.Date(c("2015-07-01", "2015-08-01", "2015-09-01",
             "2019-10-01", "2019-11-01", "2019-12-01")),
    cases = c(3391, 3146, 2632, 3521, 958, 5655)
  )
  date_dt <- prediction_dates(year_range)
  expanded_list <- expand_case_dates(case_dt, date_dt)
  expect_true("cases" %in% names(expanded_list))
  expect_true("index_range" %in% names(expanded_list))
  # One list is all ten-day periods. The other is monthly. So it's a factor of
  # 3 plus endpoints.
  expect_equal((length(expanded_list$index_range) + 2) / length(expanded_list$cases), 3)
})


constant_cases <- function(year_range) {
  case_dt <- data.table(expand.grid(y = year_range, m = 1:12, d = 1))
  case_dt[, m := as.character(m)]
  case_dt[, m := ifelse(nchar(m) == 1, paste0("0", m), m)]
  case_dt[, d := as.character(d)]
  case_dt[, d := ifelse(nchar(d) == 1, paste0("0", d), d)]
  case_dt[, date := as.Date(paste0(y, m, d), format = "%Y%m%d")]
  case_dt <- case_dt[order(date)]
  case_dt$cases <- rep(2000, nrow(case_dt))
  case_dt
}


test_that("adjusting to nothing does nothing", {
  # Make very complete case data where all values are the same.
  year_range <- 2001:2019
  case_dt <- constant_cases(year_range)
  date_dt <- prediction_dates(2000:2019)
  expanded_list <- expand_case_dates(case_dt, date_dt)
  constant_ar <- 0.4
  ar <- rep(constant_ar, nrow(date_dt))
  adjusted <- adjust_ar(ar, case_dt, date_dt)
  adj_ar <- adjusted$scaled_ar
  difference <- max((ar - constant_ar)^2)
  expect_lt(difference, 1e-4)
})
