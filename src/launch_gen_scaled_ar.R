################################################################################
## Purpose: Launch jobs to produce district-draws of scaled AR and PR
## Author: Austin Carter, aucarter@uw.edu
################################################################################

library(data.table)

### Setup
rm(list=ls())
windows <- Sys.info()[1][["sysname"]]=="Windows"
root <- ifelse(windows,"J:/","/home/j/")
user <- ifelse(windows, Sys.getenv("USERNAME"), Sys.getenv("USER"))
code.dir <- paste0(ifelse(windows, "H:", paste0("/homes/", user, "/")))

### Arguments
cluster_project <- "proj_mmc"

### Code
dist_draw_table <- data.table(expand.grid(dist = seq(135), draw = seq(100)))

run_cnt <- nrow(dist_draw_table)
# run_cnt <- 2
job_string <- paste0("qsub -l m_mem_free=2.0G -l fthread=1 -l h_rt=24:00:00 -l archive=True -q all.q ",
                     "-cwd -P ",cluster_project," ",
                     "-e /share/temp/sgeoutput/", user, "/errors ",
                     "-o /share/temp/sgeoutput/", user, "/output ",
                     "-N gen_scaled_ar ",
                     "-t 1:", run_cnt, " ",
                     code.dir, "bin/execRscript.sh ",
                     code.dir, "/dev/analytics-pipeline/gen_scaled_ar/gen_scaled_ar.R")
print(job_string)
system(job_string)  




## Combine individual files
dir <- "/ihme/malaria_modeling/projects/uganda2020/outputs/adjusted_district_draws/district_draws/"
file_list <- list.files(dir)

# dt <- rbindlist(lapply(file_list, function(file) {
#   split <- strsplit(gsub(".csv", "", file), "_")[[1]]
#   dist <- as.integer(split[1])
#   draw <- as.integer(split[2])
#   data.table(dist, draw)
# }))
# 
# dplyr::setdiff(dist_draw_table, dt)

dt <- rbindlist(lapply(file_list, function(file) {
  dt <- fread(paste0(dir, file))
}))

write.csv(dt, "/ihme/malaria_modeling/projects/uganda2020/outputs/adjusted_district_draws/adjusted_district_draws.csv", row.names = F)

## Generate annual PR summary by district
dt[, y := as.integer(substring(date,1,4))]
sum_dt <- dt[, .(pr_median = median(pr), pr_lower = quantile(pr, 0.025), pr_upper = quantile(pr, 0.975)), by = .(dist_id, y)]
write.csv(sum_dt, "/ihme/malaria_modeling/projects/uganda2020/outputs/adjusted_district_draws/annual_adjusted_district_pr_summary.csv", row.names = F)

sum_dt[dist_id == 5]
range(sum_dt$pr_median)

## Generate annual AR summary by district
hold_dt <- copy(dt)
hold_dt[, c("date", "pr") := NULL]
hold_dt[, id := seq_len(.N), by = .(y, draw, dist_id)]
cast_dt <- dcast(hold_dt, y + draw + dist_id ~ id, value.var = "ar")
mat <- 1 - as.matrix(cast_dt[, 4:ncol(cast_dt)])
cast_dt[, annual_ar := 1 - apply(t(mat), 2, prod, na.rm = T)]
sum_ar <- cast_dt[, .(ar_median = median(annual_ar), ar_lower = quantile(annual_ar, 0.025), ar_upper = quantile(annual_ar, 0.975)), by = .(dist_id, y)]
write.csv(sum_ar, "/ihme/malaria_modeling/projects/uganda2020/outputs/adjusted_district_draws/annual_adjusted_district_ar_summary.csv", row.names = F)
