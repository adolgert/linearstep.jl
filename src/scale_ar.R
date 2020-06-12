
#' Tells you about the environment, only reading from it.
read_cluster_environment <- function() {
  sge_setup <- Sys.getenv("SGE_CLUSTER_NAME") != ""
  mmc_data_dir = fs::path("", "ihme", "malaria_modeling")
  mmc_data_dir_exists = dir.exists(mmc_data_dir) # You could have this locally.
  if (mmc_data_dir_exists) {
    base_dir <- mmc_data_dir
  } else {
    # You're working from local data.
    base_dir <- fs::path_expand(fs::path("~", "data"))
  }
  # We ask the scheduler for one thread, so give data.table a clue to use 1 thread.
  thread_cnt <- as.integer(Sys.getenv("SGE_HGR_fthread"))
  if (!is.finite(thread_cnt) || thread_cnt < 1) {
    thread_cnt <- 0L  # 0 means use all available.
  }
  list(
    sge_task_job = Sys.getenv("SGE_TASK_ID") != "",
    sge_task_id = as.integer(Sys.getenv("SGE_TASK_ID")),
    sysname = factor(Sys.info()[1][["sysname"]], levels = c("Darwin", "Linux", "Windows")),
    base_data_dir = base_dir,
    compute_env = ifelse(mmc_data_dir_exists && sge_setup, "cluster", "local"),
    user = Sys.info()["effective_user"],
    thread_cnt = thread_cnt
  )
}


#' Sets values in the environment for running on the cluster.
setup_for_cluster_environment <- function(cluster_env) {
  # This script will write files with privileges so others can modify them.
  Sys.umask(002)

  data.table::setDTthreads(cluster_env$thread_cnt)
  # Turn off numerous MKL messages.
  invisible(Sys.setenv(MKL_VERBOSE = "0"))
}


### Paths
script_paths <- function(base_dir) {
  ug_dir <- fs::path(base_dir, "projects", "uganda2020")
  paths <- list(
    pfpr = fs::path(ug_dir, "outputs", "district_pfpr_draws1", "district_pfpr_draws.csv"),
    case = fs::path(ug_dir, "inputs", "DHIS2", "district_indicators_cleaned.csv"),
    shapefile = fs::path(
      ug_dir, "inputs", "uganda_subcounties_2019-wgs84", "uganda_subcounties_2019-wgs84.shp"),
    out_dir = fs::path(ug_dir, "outputs", "adam_test_draws")
  )
  paths$err <- fs::path(paths[["out_dir"]], "err.txt")
  if (!dir.exists(paths[["out_dir"]])) {
    dir.create(paths[["out_dir"]], showWarnings = FALSE, recursive = TRUE)
    # Don't tell me the warning that it already exists, but die if
    # it's not there after we tried to make it. The sleep is for
    # shared filesystems being slow.
    Sys.sleep(2)
  }  # else no need to create
  file.create(paths$err)
  for (check_idx in 1:length(paths)) {
    if (!(file.exists(paths[[check_idx]]) || !dir.exists(paths[[check_idx]]))) {
      not_found <- paste("Cannot find file", paths[[check_idx]])
      message(not_found)
      stop(not_found)
    }
  }
  paths
}


#' Ensures a file can be read or written by another group member.
group_writeable <- function(filename, user) {
  system(paste0("chown ", user, ":ihme-malaria ", filename))
}


canonical_district_names <- function(shp_path) {
  # Canonical district names are in the subcounties shape file, not the districts shape file.
  ug_dt <- as.data.table(sf::st_set_geometry(sf::st_read(shp_path), NULL))
  # Use unique because some districts have more than one polygon in the shapefile.
  dist_dict <- unique(ug_dt[, "District"])
  dist_dict <- dist_dict[, .(name = District, id = .I)]
  dist_dict
}


#' Makes a function that plots in a directory
#'
#' @param directory Where the plots go.
#' @param task_idx The index of this task, to go in the plot name.
#' @param where Whether to show the plot or save it. Values are screen, file, none.
#' @return a function(name, plot_commands). The name shouldn't have an extension.
#'   The plot commands is a block of commands to make the plot.
#' @export
plot_to_directory <- function(directory, task_idx, where) {
  stopifnot(where %in% c("screen", "file", "none"))
  if (where == "screen") {
    message("plotting to screen")
    f <- function(name, plot_commands) {
      eval(substitute(plot_commands), parent.frame())
    }
  } else if (where == "file") {
    message(paste("plotting to", directory))
    f <- function(name, plot_commands) {
      pdf(
        file = fs::path(directory, paste0(task_idx, "_", name), ext = "pdf"),
        width = 6,
        height = 6
        )
      eval(substitute(plot_commands), parent.frame())
      invisible(dev.off())
    }
  } else {
    message("plotting disabled")
    f <- function(name, plot_commands) {}
  }
  f
}


task_work <- function(district_table, dist_draw_idx, median_only = FALSE) {
  stopifnot(is.numeric(dist_draw_idx))
  district_cnt <- nrow(district_table)
  if (median_only) {
    dist_draw_table <- data.table(expand.grid(dist = seq(district_cnt), draw = 0L))
  } else {
    dist_draw_table <- data.table(expand.grid(dist = seq(district_cnt), draw = seq(100L)))
  }
  ## Identify our district and draw to work on, by name and ID.
  district = dist_draw_table[dist_draw_idx]$dist
  list(
    district = district,
    district_name = district_table[id == district]$name,
    draw_id = dist_draw_table[dist_draw_idx]$draw
  )
}


### Epi Functions


#' The original peekahead to give attack rate from PfPR.
#'
#' @param pr a series of PfPR as a vector
#' @param parameters Not really any parameters
pr2ar_peek <- function(pr, parameters = NULL) {
  list(alpha = pmin((1 - exp(-10 / 200)) * pr / (1 - pr), .99))
}


#' A peekahead to give attack rate from PfPR, using more accurate formula.
#'
#' @param pr a series of PfPR as a vector
#' @param parameters Not really any parameters
pr2ar_strang <- function(pr, parameters = NULL) {
  Q <- exp(-10 / 200)  # 10 day step, 200 day recovery rate.
  Qrt <- exp(-5 / 200)
  list(alpha = (1 - Q) * pr / (Qrt * (2 - pr - Qrt)))
}


#' Returns a smoothed daily time-series
smoothX = function(x, in_tstep, out_tstep, b) {
  ksmooth(1:length(x), x, kernel = "normal",
          bandwidth = (b / in_tstep),
          n.points = in_tstep*(length(x) - 1)/out_tstep + 1)
}


# Prep PfPR so that it is for this district, draw, and ordered by date.
prep_pfpr <- function(pfpr_path, district, draw_id, limit = Inf) {
  input_pr_dt <- fread(pfpr_path)[district_id == district]
  input_pr_dt$date <- as.Date(paste0(input_pr_dt$year, "-", input_pr_dt$month, "-01"))
  if (median_only) {
    # Compute the median of the draws for this district.
    pr_dt <- input_pr_dt[, .(value = median(value)), by = .(date)]
  } else {
    pr_dt <- input_pr_dt[draw == draw_id, .(date, value)]
  }
  pr_dt[, date := as.Date(date)]
  pr_dt <- pr_dt[order(date)]
  if (is.finite(limit)) {
    return(pr_dt[1:limit])
  } else {
    pr_dt
  }
}


#' Read case data for this district with minor elimination of unusable records.
#'
#' @param case_path filesystem path to the case data.
#' @param district_name Capitalized district name, as a string.
#' @return A data.dable(date, cases), where the date column is YYY-MM-DD and
#'     casess are integers.
case_data <- function(case_path, district_name) {
  case_dt <- fread(case_path)[organisationunitname == district_name]
  setnames(case_dt, "Total_cases", "cases")
  case_dt[, date := as.Date(periodid)]
  case_dt <- case_dt[!(is.na(date) | is.na(cases)), .(date, cases)]
  case_dt
}


#' At what exact dates do we make predictions? Make a data.table of those.
prediction_dates <- function(year_range) {
  # Generate full set of dates
  date_dt <- data.table(expand.grid(y = year_range, m = 1:12, d = c(1, 11, 21)))
  date_dt[, m := as.character(m)]
  date_dt[, m := ifelse(nchar(m) == 1, paste0("0", m), m)]
  date_dt[, d := as.character(d)]
  date_dt[, d := ifelse(nchar(d) == 1, paste0("0", d), d)]
  date_dt[, date := as.Date(paste0(y, m, d), format = "%Y%m%d")]
  date_dt[order(date)]
}


#' The cases are recorded for a subset of dates, with missingness, so interpolate.
expand_case_dates <- function(case_dt, date_dt) {
  # Interpolate missing dates
  l_idx <- which(date_dt$date == min(case_dt$date))
  u_idx <- which(date_dt$date == max(case_dt$date))
  comp_dt <- date_dt[l_idx:u_idx][d == "01"]
  cases <- approx(x = case_dt$date, y = case_dt$cases, xout = comp_dt$date)$y
  list(cases = cases, index_range = l_idx:u_idx)
}


#' Look for short-term trends as a ratio over long-term trends.
case_ratio <- function(cases, date_dt, index_range) {
  # Smooth
  smooth_cases_short <- smoothX(cases, 30, 10, 30)$y
  smooth_cases_long <- smoothX(cases, 30, 10, 200)$y

  # Case ratio
  case_r <- smooth_cases_short/smooth_cases_long
  case_r_dt <- cbind(date_dt[index_range], case_r)
  list(
    smooth_cases_short = smooth_cases_short,
    smooth_cases_long = smooth_cases_long,
    case_ratio_dt = case_r_dt
  )
}


#' Compare case ratio to median in order to assign seasonal signal.
seasonality_from_case_ratio <- function(case_ratio_dt, date_dt) {
  season_dt <- case_ratio_dt[, .(case_r = median(case_r)), by = .(m, d)]
  merge_dt <- merge(
    date_dt[, .(date, y, m, d)],
    case_ratio_dt[, .(date, case_r)],
    by = "date",
    all.x = T
  )
  fill_dt <- merge(merge_dt[is.na(case_r), .(date, y, m, d)], season_dt, by = c("m", "d"))
  full_case_ratio_dt <- rbind(merge_dt[!is.na(case_r)], fill_dt)[order(date)]
  full_case_ratio_dt
}


#' Use case data to adjust attack rate over a set of observation times.
adjust_ar <- function(smooth_ar, case_dt, date_dt) {
  # This creates an average daily force of infection for each time period.
  # It's not the integrated FOI of the ten-day step. It's per day.
  smooth_foi <- ar2foi(smooth_ar)
  smooth_foi_na_cnt <- sum(is.na(smooth_foi))
  if (smooth_foi_na_cnt > 0L) {
    stop("na in smooth_foi")
  }

  case_dates_list <- expand_case_dates(case_dt, date_dt)
  case_ratio_list <- with(case_dates_list, case_ratio(cases, date_dt, index_range))
  full_case_ratio_dt <- seasonality_from_case_ratio(
    case_ratio_list$case_ratio, date_dt)

  ## Scale FoI and calculate resulting PR
  ## note that smooth_foi is longer and getting truncated
  scaled_foi <- smooth_foi*full_case_ratio_dt$case_r[1:length(smooth_foi)]

  list(
    smooth_ar = smooth_ar,
    smooth_foi = smooth_foi,
    smooth_cases_short = case_ratio_list$smooth_cases_short,
    smooth_cases_long = case_ratio_list$smooth_cases_long,
    case_ratio_dt = case_ratio_list$case_ratio_dt,
    index_range = case_dates_list$index_range,
    full_case_ratio_dt = full_case_ratio_dt,
    scaled_foi = scaled_foi,
    scaled_ar = foi2ar(scaled_foi)
  )
}
