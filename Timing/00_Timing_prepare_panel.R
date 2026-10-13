#!/usr/bin/env Rscript
# Input: original panel, hourly ERA5 and PM2.5; output: complete adjusted panel and match audit.
# nohup Rscript 00_Timing_prepare_panel.R > 00_Timing_prepare_panel.log 2>&1 &
# ps -ef | grep 00_Timing_prepare_panel.R
# pkill -f 00_Timing_prepare_panel.R

.code_dir <- Sys.getenv("TIMING_REVISION_DIR", unset = "")
if (!nzchar(.code_dir)) {
  .script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  .code_dir <- if (length(.script_arg)) dirname(normalizePath(sub("^--file=", "", .script_arg[1]), winslash = "/")) else getwd()
}
source(file.path(.code_dir, "timing_helpers.R"))
source(file.path(.code_dir, "hourly_weather_controls.R"))

tr_weather_prepare_panel()
