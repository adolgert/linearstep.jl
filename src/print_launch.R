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
# Just the median surfaces
dist_draw_table <- data.table(expand.grid(dist = seq(135), draw = 0))

run_cnt <- nrow(dist_draw_table)
# run_cnt <- 2
job_string <- paste0("qsub -l m_mem_free=2.0G -l fthread=1 -l h_rt=4:00:00 -l archive=True -q all.q ",
                     "-cwd -P ",cluster_project," ",
                     "-e /share/temp/sgeoutput/", user, "/errors ",
                     "-o /share/temp/sgeoutput/", user, "/output ",
                     "-N gen_scaled_ar ",
                     "-t 1:", run_cnt, " ",
                     code.dir, "bin/execRscript.sh ",
                     code.dir, "/dev/analytics-pipeline/gen_scaled_ar/gen_scaled_ar.R")
print(job_string)
