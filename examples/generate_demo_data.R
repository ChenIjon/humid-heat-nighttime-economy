#!/usr/bin/env Rscript
# 用途：生成完全模拟的城市日面板；不读取、拟合或复制真实交易数据。
# 输入：下方模拟参数；输出：main_panel_synthetic.csv，UTF-8。
# 运行：Rscript examples/generate_demo_data.R

DEMO_SEED <- 20261013L
DEMO_CITIES <- c("Beijing", "Shanghai", "Guangzhou")
DEMO_START_DATE <- as.Date("2024-05-01")
DEMO_N_DAYS <- 120L
DEMO_CSV_NAME <- "main_panel_synthetic.csv"

generate_demo_panel <- function() {
  set.seed(DEMO_SEED)
  dates <- DEMO_START_DATE + seq_len(DEMO_N_DAYS) - 1L
  d <- expand.grid(city = DEMO_CITIES, DATE = dates, stringsAsFactors = FALSE)
  d <- d[order(d$city, d$DATE), ]; rownames(d) <- NULL
  n <- nrow(d); city_id <- match(d$city, DEMO_CITIES)
  day <- as.integer(d$DATE - DEMO_START_DATE)
  weekday <- as.POSIXlt(d$DATE)$wday
  d$calendar_type <- c("Sunday", "Monday", "Tuesday", "Wednesday", "Thursday",
                       "Friday", "Saturday")[weekday + 1L]
  d$calendar_type[day == 0L] <- "Statutory holiday"
  d$calendar_type[day == 10L] <- "Make-up workday"
  d$TEMP <- 25 + 4 * sin(2 * pi * day / 120) + city_id + rnorm(n, 0, 2)
  d$RH <- pmin(95, pmax(25, 65 + 10 * cos(2 * pi * day / 80) + rnorm(n, 0, 8)))
  d$WS <- 0.5 + rgamma(n, shape = 3, rate = 1.2)
  d$PRECIP <- ifelse(runif(n) < 0.65, 0, rgamma(n, shape = 1.5, rate = 0.4))
  d$PM25 <- 5 + rgamma(n, shape = 4, rate = 0.16)
  # 以下系数只是演示设定，不代表论文估计或现实关系。
  errors <- unlist(lapply(DEMO_CITIES, function(city)
    as.numeric(stats::arima.sim(list(ar = 0.3), n = DEMO_N_DAYS, sd = 0.08))))
  log_mu <- log(2500) + 0.1 * city_id + 0.008 * (d$TEMP - 25) -
    0.0004 * (d$RH - 65)^2 - 0.005 * d$PRECIP +
    0.05 * (weekday %in% c(0L, 6L)) + errors
  d$model_cnt <- pmax(1L, rnbinom(n, mu = exp(log_mu), size = 80))
  average_spending <- rgamma(n, shape = 120,
    scale = exp(log(80) + 0.03 * city_id + 0.003 * (d$TEMP - 25)) / 120)
  d$model_amt <- d$model_cnt * average_spending
  d$model_amt_per_cnt <- d$model_amt / d$model_cnt
  d[, c("city", "DATE", "TEMP", "RH", "WS", "PRECIP", "PM25",
        "calendar_type", "model_cnt", "model_amt", "model_amt_per_cnt")]
}

write_demo_panel <- function(output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  d <- generate_demo_panel()
  path <- file.path(output_dir, DEMO_CSV_NAME)
  write.csv(d, path, row.names = FALSE, fileEncoding = "UTF-8")
  message("Synthetic data only; no real input files.\nRows: ", nrow(d),
          "; cities: ", length(unique(d$city)), "; dates: ", min(d$DATE),
          " to ", max(d$DATE), "\nOutput: ", path,
          "\nOverwrite: same-name synthetic CSV")
  invisible(d)
}

if (sys.nframe() == 0L) {
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/"))
  write_demo_panel(script_dir)
}
