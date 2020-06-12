################################################################################
## Purpose: Read in district-level PfPR, smooth it, calculate the associated
##          attack-rate/FoI, then read in disctrict cases and modify the FoI
##          with the temporal signal
## Author: Austin Carter, aucarter@uw.edu
################################################################################

### Setup
rm(list = ls())

## Packages
library(data.table)
library(ggplot2)
library(pr2ar)
library(sf)

source("adam.R")
source("scale_ar.R")


cluster_env <- read_cluster_environment()
paths <- script_paths(cluster_env$base_data_dir)

## Arguments
if (cluster_env$compute_env == "cluster") {
  setup_for_cluster_environment(cluster_env)
}

if (cluster_env$sge_task_job) {
  dist_draw_idx <- cluster_env$sge_task_id
  pr_cutoff <- Inf
  median_only <- TRUE
  use_adam <- TRUE
  plot_to <- "file"
} else {
  # Only process a few PfPR values so it's faster.
  pr_cutoff <- 25L
  median_only <- TRUE
  if (median_only) {
    dist_draw_idx <- 115
  } else {
    dist_draw_idx <- 731
  }
  use_adam <- FALSE
  plot_to <- "screen"
}
save_plot <- plot_to_directory(paths[["out_dir"]], dist_draw_idx, plot_to)

### Tables
dist_dict <- canonical_district_names(paths[["shapefile"]])
work <- task_work(dist_dict, dist_draw_idx, median_only)
# district = 94
# draw_id = 1

case_dt <- case_data(paths$case, work$district_name)
pr_dt <- prep_pfpr(paths$pfpr, work$district, work$draw_id, limit = pr_cutoff)
save_plot("pr_raw", {
  plot(pr_dt$date, pr_dt$value, type = "l", ylim = c(0, 0.6))
  title(main = "Raw PR")
})
# Smooth
smooth_pr <- smoothX(pr_dt$value, 30, 10, 360)$y

# Convert to AR
# This trajectory is a list with several components, one of which is attack rate.
adam_parameters <- makeParameters_Adam()
if (use_adam) {
  ar_trajectory <- pr2ar_Adam(smooth_pr, adam_parameters)
  peek_ar <- NULL
} else {
  ar_trajectory <- pr2ar_strang(smooth_pr, adam_parameters)
  peek_ar <- ar_trajectory$alpha
}
smooth_ar <- ar_trajectory$alpha
trajectory_na_cnt <- sum(is.na(smooth_ar))
if (trajectory_na_cnt > 0L) {
  write(paste(dist_draw_idx, "na trajectory", trajectory_na_cnt, "\n"), paths$err)
  write(smooth_ar, paths$err)
  stop("na in trajectory")
}

date_dt <- prediction_dates(2000:2019)
adjusted <- adjust_ar(smooth_ar, case_dt, date_dt)
scaled_ar <- adjusted$scaled_ar
with(adjusted, {
  save_plot("pr_in", {
    if (is.null(peek_ar)) {
      peek_ar <- pr2ar_strang(smooth_pr, adam_parameters)$alpha
    }
    plot(smooth_pr, type = "l", ylim = c(0, 1))
    lines(smooth_ar, col = "blue")
    lines(smooth_foi, col = "red")
    lines(peek_ar, col = "green")
    title(main = "pr_ar_in")
    legend(2, 0.9, c("pr", "ar", "foi", "peek_ar"),
           col = c("black", "blue", "red", "green"), lwd = 2)
  })
  save_plot("smooth", {
    plot(smooth_cases_short, type = "l", col = "blue")
    lines(smooth_cases_long)
    title(main = "smoothing")
  })
  save_plot("caseratio", {
    plot(case_ratio_dt$date, case_ratio_dt$case_r, type = "l")
    title(main = "case ratio")
  })
})

scaled_pr <- ar2pr_series_Adam(scaled_ar, adam_parameters)
save_plot("adj_ar", {
  plot(smooth_ar, type = "l", ylim = 0:1)
  title(main = "scaled ar and pr")
  legend(
    2,
    0.9,
    c("ar in", "ar out", "pr out", "pr in"),
    col = c("black", "blue", "red", "green"),
    lwd = 2)
  lines(scaled_ar, col = "blue")
  lines(scaled_pr, col = "red")
  lines(smooth_pr, col = "green")
})

## Construct table
out_dt <- with(
  work,
  data.table(
    dist_id = district,
    draw = draw_id,
    date = date_dt[1:length(scaled_pr)]$date,
    pr = scaled_pr,
    ar = scaled_ar)
)
pr_out_name <- fs::path(paths$out_dir, paste0(work$district, "_", work$draw_id, ".csv"))
write.csv(out_dt, pr_out_name, row.names = F)
if (cluster_env$compute_env == "cluster") {
  group_writeable(pr_out_name, cluster_env$user)
}
