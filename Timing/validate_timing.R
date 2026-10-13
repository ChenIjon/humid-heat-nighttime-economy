#!/usr/bin/env Rscript
# Synthetic validation; outputs only in a temporary directory. Run: Rscript validate_timing.R
.code_dir <- Sys.getenv("TIMING_REVISION_DIR", unset = "")
if (!nzchar(.code_dir)) {
  .script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  .code_dir <- if (length(.script_arg)) dirname(normalizePath(sub("^--file=", "", .script_arg[1]), winslash = "/")) else getwd()
}
source(file.path(.code_dir, "timing_helpers.R"))
source(file.path(.code_dir, "hourly_weather_controls.R"))

LOCAL_BASE <- .code_dir
TEST_DIR <- tempfile("timing_v2_smoke_"); dir.create(TEST_DIR)
TIMING_ROOT <- TEST_DIR
SOURCE_PANEL <- file.path(TEST_DIR, "source.csv")
ERA5_FILE <- file.path(TEST_DIR, "era.csv")
PM25_FILE <- file.path(TEST_DIR, "pm.csv")
PANEL_FILE <- file.path(TEST_DIR, "data", "hourly_profile_panel_adjusted.csv")
PANEL_META_FILE <- file.path(TEST_DIR, "data", "panel_metadata_adjusted.rds")

# 北京/拉萨/乌鲁木齐夜间起点与跨午夜匹配；故意让一夜缺PM小时。
src <- CJ(city = c("Beijing", "Lasa", "Wulumuqi"), DATE = as.Date(c("2024-06-01", "2024-06-02")), hour_position = 0:11)
src[, night_start_hour := fifelse(city == "Lasa", 19L, fifelse(city == "Wulumuqi", 20L, 18L))]
src[, night_id := paste(city, DATE, sep = "__")]
src[, actual_clock_hour := (night_start_hour + hour_position) %% 24L]
src[, city_hour_id := paste(city, hour_position, sep = ".")]
src[, c("TEMP", "RH", "hour_cnt", "night_total_cnt", "hour_amt", "night_total_amt", "calendar_type", "sunset_minutes_CST") :=
  list(25, 70, 1, 12, 10, 120, "Saturday", 1100)]
fwrite(src, SOURCE_PANEL)
era <- src[, .(city, date = DATE + (night_start_hour + hour_position) %/% 24L,
  hour = actual_clock_hour, precipitation_mm = as.numeric(hour_position) / 10,
  wind_speed_mps = 1 + as.numeric(hour_position) / 10)]
pm <- era[, .(city, date, hour, PM25 = 10 + hour)]
pm <- pm[!(city == "Beijing" & date == as.Date("2024-06-02") & hour == 18L)]
fwrite(era, ERA5_FILE); fwrite(pm, PM25_FILE)
d <- tr_weather_prepare_panel()
stopifnot(nrow(d) == 60L, uniqueN(d$night_id) == 5L,
  all(d$PRECIP_hourly == d$hour_position / 10),
  all(d$WS_hourly == 1 + d$hour_position / 10),
  all(d$PM25_hourly == 10 + d$actual_clock_hour),
  all(abs(d$PRECIP_night_mean - 0.55) < 1e-12),
  all(abs(d[, sum(PRECIP_hour_dev), by = night_id]$V1) < 1e-12))
stopifnot(nrow(tr_weather_read_panel()) == 60L)
OC <- tr_outcome_spec(); valid <- rep(TRUE, nrow(d)); valid[1L] <- FALSE
filtered <- tr_weather_complete_nights(copy(d), valid, OC)
stopifnot(nrow(filtered) == 48L, uniqueN(filtered$night_id) == 4L)
cat("PASS: 城市夜间窗口、跨午夜、整夜缺失剔除、夜均值分解\n")

# 随机小时控制的条件设计秩；零支持降水列需被安全剔除。
set.seed(42)
r <- CJ(night_id = 1:100, hour = 0:11)
for (v in WEATHER_VARIABLES) {
  r[, (paste0(v, "_hourly")) := rgamma(.N, shape = 2, rate = 1)]
  r[, (paste0(v, "_night_mean")) := mean(get(paste0(v, "_hourly"))), by = night_id]
  r[, (paste0(v, "_hour_dev")) := get(paste0(v, "_hourly")) - get(paste0(v, "_night_mean"))]
}
WX <- tr_weather_design(r, TEST_DIR)
stopifnot(length(WX) == 69L)
A <- as.matrix(r[, ..WX]); first <- match(r$night_id, r$night_id)
stopifnot(qr(A - A[first, ])$rank == 69L)
CONTROL_MODE <- "hourly_hour"
WX36 <- tr_weather_design(copy(r), TEST_DIR)
stopifnot(length(WX36) == 36L)
CONTROL_MODE <- "mean_deviation_hour"
X <- as(cbind("(Intercept)" = 1, A, "s(city_part_id):TEMP" = rnorm(nrow(r))), "dgCMatrix")
X[, "weather_PRECIP_hour_dev_hour_0"] <- 0
k <- tr_reduce_frozen_design(X, colnames(X), file.path(TEST_DIR, "columns.rds"),
  city_hour_baselines = "OnlyCity.0", night_groups = r$night_id)
stopifnot(!"weather_PRECIP_hour_dev_hour_0" %in% k, !"(Intercept)" %in% k,
  sum(grepl("^weather_", k)) == 68L)
cat("PASS: 69/36控制列与条件秩、无支持降水列剔除\n")

# 冻结Newton求解器只从模型脚本抽取定义，在100个合成夜晚上验证。
expr <- parse(file.path(LOCAL_BASE, "02_Timing_model.R"), encoding = "UTF-8")
for (e in expr) if (is.call(e) && identical(e[[1L]], as.name("<-")) &&
  is.symbol(e[[2L]]) && as.character(e[[2L]]) %in% c("solve_cond", "gain_tol_at")) eval(e)
RIDGE_REL <- 1e-10; GRAD_TOL <- 1e-6; NEWTON_GAIN_TOL <- 1e-6; INNER_MAXIT <- 25L
XX <- as(A, "dgCMatrix"); bb <- rnorm(ncol(A), sd = 0.02)
eta <- as.numeric(XX %*% bb); hh <- matrix(eta, nrow = 12L)
pr <- exp(hh); pr <- sweep(pr, 2, colSums(pr), "/")
y <- as.numeric(pr) * 100
fit <- solve_cond(XX, y, rep(100, nrow(r)), r$night_id, r$hour, rep(0, ncol(A)))
stopifnot(is.finite(fit$nll), max(abs(fit$beta - bb)) < 1e-5)
weather_test <- tr_weather_tests(fit$beta, solve(fit$H), colnames(XX), "cnt", "hessian", TEST_DIR)
stopifnot(nrow(weather_test$coefficients) == 69L, nrow(weather_test$joint) == 9L,
  setequal(weather_test$coefficients$variable, WEATHER_VARIABLES),
  setequal(weather_test$coefficients$component, c("night_mean", "hour_dev")))
cat("PASS: 合成条件模型恢复参数与小时/联合Wald输出\n")

# 四阶段真实mgcv设计构建必须一致；反事实必须冻结所有气象控制。
fixture <- copy(r)
fixture[, hour_position := hour]
fixture[, city := "Beijing"]
fixture[, DATE := as.Date("2024-06-01") + night_id - 1L]
fixture[, city_hour_id := paste(city, hour, sep = ".")]
fixture[, TEMP := rnorm(1L, 25, 4), by = night_id]
fixture[, RH := runif(1L, 40, 80), by = night_id]
fixture[, sunset_minutes_CST := runif(1L, 1080, 1200), by = night_id]
fixture[, calendar_type := sample(TR_CALENDAR_LEVELS, 1L), by = night_id]
fixture[, hour_cnt := rpois(.N, 15) + 1L]
fixture[, night_total_cnt := sum(hour_cnt), by = night_id]
fixture[, hour_amt := hour_cnt * 10]
fixture[, night_total_amt := sum(hour_amt), by = night_id]
stage_cn <- list(); model_env <- NULL
for (outcome in c("cnt", "amt")) for (stage in c("01_Timing_penalty.R", "02_Timing_model.R", "03_Timing_inference.R", "04_Timing_grid_ci.R")) {
  ee <- new.env(parent = globalenv()); ee$OUT <- file.path(TEST_DIR, paste(outcome, stage, sep = "_"))
  dir.create(ee$OUT); ee$OC <- tr_outcome_spec(outcome); ee$SPLINE_DF <- 3L; ee$N_CORES <- 1L; ee$N_BOOT <- 4L
  steps <- parse(file.path(LOCAL_BASE, stage), encoding = "UTF-8")
  active <- FALSE
  for (e in steps) {
    assignment <- is.call(e) && identical(e[[1L]], as.name("<-")) && is.symbol(e[[2L]])
    lhs <- if (assignment) as.character(e[[2L]]) else ""
    if (!active && lhs == "d" && identical(e[[3L]][[1L]], as.name("tr_weather_read_panel"))) {
      ee$d <- copy(fixture); active <- TRUE; next
    }
    if (!active) next
    if (lhs %in% c("iT", "re_idx")) break
    # 上游文件一致性另有独立测试，此处仅验证共同设计块。
    if (grepl("(lamL|model)\\$cn", paste(deparse(e), collapse = ""))) next
    eval(e, envir = ee)
  }
  stage_cn[[paste(outcome, stage)]] <- ee$cn
  stopifnot(sum(grepl("^weather_", ee$cn)) == 69L)
  if (stage == "02_Timing_model.R") model_env <- ee
}
stopifnot(all(vapply(stage_cn, identical, logical(1), stage_cn[[1L]])))
for (e in parse(file.path(LOCAL_BASE, "02_Timing_model.R"), encoding = "UTF-8")) {
  if (is.call(e) && identical(e[[1L]], as.name("<-")) && is.symbol(e[[2L]]) &&
      as.character(e[[2L]]) %in% c("block2", "B2", "make_nd", "frozen_X")) eval(e, envir = model_env)
}
ctrl <- which(model_env$B2 == "ctrl")
for (exposure in c("TEMP", "RH")) {
  nd <- model_env$make_nd(exposure, rep(mean(model_env$d[[exposure]]) + 1, nrow(model_env$d)))
  cf <- model_env$frozen_X(nd, exposure)
  stopifnot(max(abs(cf[, ctrl, drop = FALSE] - model_env$X[, ctrl, drop = FALSE])) == 0)
}
cat("PASS: 四阶段mgcv设计一致、反事实冻结小时气象控制\n")

# 重复小时键必须拒绝；原始小时数据变化必须拒绝旧面板。
fwrite(rbind(pm, pm[1L]), PM25_FILE)
err <- try(tr_weather_hour_table(PM25_FILE,
  c(city = "city", date = "date", hour = "hour", PM25_hourly = "PM25"), "PM25_hourly", "重复测试"), silent = TRUE)
stopifnot(inherits(err, "try-error"))
stopifnot(inherits(try(tr_weather_read_panel(), silent = TRUE), "try-error"))
cat("PASS: 重复键和过期数据保护\nALL PASS. 临时输出：", TEST_DIR, "\n")
