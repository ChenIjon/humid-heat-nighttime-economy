#!/usr/bin/env Rscript
# 用途：使用完全模拟的数据演示原主模型、EWME 和少量日期块 bootstrap。
# 输入：examples/main_panel_synthetic.csv；输出：examples/output_demo/。
# 运行：Rscript examples/run_demo.R；不读取正式 data/main_panel.csv。
# nohup Rscript examples/run_demo.R > demo.log 2>&1 &
# ps -ef | grep examples/run_demo.R
# pkill -f examples/run_demo.R

DEMO_N_BOOT <- 4L
DEMO_N_CORES <- 1L
DEMO_SEED <- 20261013L
DEMO_INPUT <- "main_panel_synthetic.csv"
DEMO_OUTPUT <- "output_demo"

run_demo <- function(example_dir) {
  source(file.path(dirname(example_dir), "Main_model.R"), local = globalenv(), encoding = "UTF-8")
  input <- file.path(example_dir, DEMO_INPUT)
  output <- file.path(example_dir, DEMO_OUTPUT)
  if (!file.exists(input)) stop("Synthetic input missing; run examples/generate_demo_data.R first.")
  panel <- read.csv(input, stringsAsFactors = FALSE)
  config <- main_model_defaults()
  config$n_boot <- DEMO_N_BOOT; config$n_cores <- DEMO_N_CORES
  config$seed <- DEMO_SEED; config$force_refit <- TRUE
  message("DEMONSTRATION ONLY: fully synthetic data; estimates are not research results.",
          "\nInput: ", input, "\nInput files: 1\nOutput: ", output,
          "\nBootstrap replicates per outcome: ", config$n_boot,
          "\nOverwrite: same-name demo outputs; no cache reuse")
  t0 <- Sys.time()
  result <- main_model_run(panel, level_model_outcomes(), output, config)
  stopifnot(nrow(result$data) == nrow(panel), nrow(result$diagnostics) == 3L,
            all(result$diagnostics$Converged),
            all(result$ewme$N_boot_success == DEMO_N_BOOT),
            all(is.finite(result$ewme$CI_lower)), all(is.finite(result$ewme$CI_upper)))
  main_model_write_report(result, output, title = "SYNTHETIC DEMO - NOT RESEARCH RESULTS",
                          source_line = "Fully synthetic input; no real transactions")
  writeLines(c("This records the demonstration runtime, not the server research environment.",
               capture.output(sessionInfo())), file.path(output, "demo_session_info.txt"))
  message("DEMO PASS: ", nrow(result$data), " rows; 3 outcomes; ", DEMO_N_BOOT,
          " successful bootstrap fits per outcome.\nElapsed seconds: ",
          round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
          "\nSaved to: ", output)
  invisible(result)
}

if (sys.nframe() == 0L) {
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  example_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/"))
  run_demo(example_dir)
}
