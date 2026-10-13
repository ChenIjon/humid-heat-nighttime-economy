# -*- coding: UTF-8 -*-
# 小时协变量配置与工具；ERA5/PM2.5 路径在下方设置。
# 输入：旧版小时结局面板、北京时间 ERA5 小时表、PM2.5 城市小时表。
# 输出：完整夜调整面板、匹配审计、标准化记录和协变量检验。
suppressPackageStartupMessages(library(data.table))

SOURCE_PANEL <- file.path(TIMING_DATA_DIR, "hourly_profile_panel.csv")
ERA5_FILE <- "/data/ccpeng/NEI/ERA5_population_weighted_hourly_long.csv"
PM25_FILE <- "/data/ccpeng/NEI/PM25_city_hourly_2024.csv"
PM25_COLUMNS <- c(city = "city", date = "date", hour = "hour", value = "PM25")
CONTROL_MODE <- "mean_deviation_hour"  # 69列；"hourly_hour" 为36列
WEATHER_VARIABLES <- c("PRECIP", "WS", "PM25")
WEATHER_HOURS <- 0:11
PANEL_FILE <- file.path(TIMING_DATA_DIR, "hourly_profile_panel_adjusted.csv")
PANEL_META_FILE <- file.path(TIMING_DATA_DIR, "panel_metadata_adjusted.rds")
TIMING_MODEL_VERSION <- paste0("conditional_city_quarter_weather_v2_", CONTROL_MODE)

tr_weather_settings <- function() list(mode = CONTROL_MODE, weather = WEATHER_VARIABLES,
  hours = WEATHER_HOURS, source = SOURCE_PANEL, era5 = ERA5_FILE,
  pm25 = PM25_FILE, pm25_columns = PM25_COLUMNS)

tr_weather_hour_table <- function(path, columns, values, label) {
  if (!file.exists(path)) stop(label, "文件不存在：", path,
    "。请提供真实小时数据，不可用整夜或日均值代替。")
  x <- fread(path)
  tr_require_columns(x, unname(columns), label)
  x <- x[, unname(columns), with = FALSE]
  setnames(x, unname(columns), names(columns))
  x[, city := as.character(city)]
  x[, date := tr_parse_date(date)]
  hour_num <- suppressWarnings(as.numeric(x$hour))
  if (anyNA(x$date) || anyNA(x$city) || any(!nzchar(x$city)) ||
      any(!is.finite(hour_num) | hour_num != floor(hour_num) | hour_num < 0 | hour_num > 23))
    stop(label, "城市、日期或小时键不合法；小时应为北京时间0--23。")
  x[, hour := as.integer(hour_num)]
  for (v in values) set(x, j = v, value = suppressWarnings(as.numeric(x[[v]])))
  tr_assert_unique(x, c("city", "date", "hour"), label)
  x
}

tr_weather_prepare_panel <- function() {
  t0 <- Sys.time()
  if (identical(normalizePath(SOURCE_PANEL, winslash = "/", mustWork = FALSE),
                normalizePath(PANEL_FILE, winslash = "/", mustWork = FALSE)))
    stop("原始面板与调整面板必须使用不同路径。")
  if (!CONTROL_MODE %in% c("mean_deviation_hour", "hourly_hour")) stop("未知气象控制模式")
  paths <- c(SOURCE_PANEL, ERA5_FILE, PM25_FILE)
  cat("输入文件数：", length(paths), "\n", paste(paths, collapse = "\n"),
      "\n输出：", PANEL_FILE, "\n覆盖策略：覆盖同名调整面板和审计文件，原始面板只读\n")
  if (any(!file.exists(paths))) stop("输入文件缺失：", paste(paths[!file.exists(paths)], collapse = "; "))
  d <- fread(SOURCE_PANEL)
  tr_require_columns(d, c("city", "DATE", "night_id", "night_start_hour", "hour_position",
    "actual_clock_hour", "TEMP", "RH", "hour_cnt", "night_total_cnt", "hour_amt",
    "night_total_amt", "city_hour_id", "calendar_type", "sunset_minutes_CST"), "原版面板")
  d[, DATE := tr_parse_date(DATE)]
  if (anyNA(d$DATE) || any(!is.finite(d$hour_position)) ||
      any(!d$hour_position %in% WEATHER_HOURS) || any(!d$night_start_hour %in% 0:23)) stop("原版夜间日期/小时不合法")
  tr_assert_unique(d, c("night_id", "hour_position"), "原版面板")
  night_check <- d[, .(n = .N, n_hour = uniqueN(hour_position), n_date = uniqueN(DATE),
    n_city = uniqueN(city), n_start = uniqueN(night_start_hour)), by = night_id]
  if (any(night_check$n != 12L | night_check$n_hour != 12L | night_check$n_date != 1L |
          night_check$n_city != 1L | night_check$n_start != 1L)) stop("原版面板包含不完整或键不一致的城市夜晚")
  d[, weather_date := DATE + as.integer((night_start_hour + hour_position) %/% 24L)]
  d[, weather_hour := as.integer((night_start_hour + hour_position) %% 24L)]
  if (any(d$weather_hour != d$actual_clock_hour)) stop("实际时钟小时与夜间起点不一致")
  # 三个小时气象字段独立保存；原面板TEMP/RH与原PRECIP/WS/PM25均不替换。
  era <- tr_weather_hour_table(ERA5_FILE,
    c(city = "city", date = "date", hour = "hour", PRECIP_hourly = "precipitation_mm", WS_hourly = "wind_speed_mps"),
    c("PRECIP_hourly", "WS_hourly"), "ERA5小时表")
  pm_columns <- c(city = PM25_COLUMNS[["city"]], date = PM25_COLUMNS[["date"]],
    hour = PM25_COLUMNS[["hour"]], PM25_hourly = PM25_COLUMNS[["value"]])
  pm <- tr_weather_hour_table(PM25_FILE, pm_columns, "PM25_hourly", "PM2.5小时表")
  cat("ERA5行数：", nrow(era), " PM2.5行数：", nrow(pm), "；开始小时键匹配\n")
  for (v in paste0(WEATHER_VARIABLES, "_hourly")) if (v %in% names(d)) d[, (v) := NULL]
  d[era, on = .(city, weather_date = date, weather_hour = hour),
    c("PRECIP_hourly", "WS_hourly") := .(i.PRECIP_hourly, i.WS_hourly)]
  d[pm, on = .(city, weather_date = date, weather_hour = hour), PM25_hourly := i.PM25_hourly]
  for (v in WEATHER_VARIABLES) d[, (paste0("bad_", v)) := !is.finite(get(paste0(v, "_hourly"))) |
    get(paste0(v, "_hourly")) < 0]
  audit <- d[, .(city = city[1L], DATE = DATE[1L], n_hours = .N,
    bad_PRECIP_hours = sum(bad_PRECIP), bad_WS_hours = sum(bad_WS), bad_PM25_hours = sum(bad_PM25)), by = night_id]
  audit[, retained := bad_PRECIP_hours == 0L & bad_WS_hours == 0L & bad_PM25_hours == 0L]
  dir.create(dirname(PANEL_FILE), recursive = TRUE, showWarnings = FALSE)
  tr_write_csv(audit, file.path(dirname(PANEL_FILE), "weather_match_audit.csv"))
  d <- d[night_id %in% audit[retained == TRUE, night_id]]
  if (!nrow(d)) stop("没有三种小时气象均完整的夜晚；检查匹配审计。")
  d[, c("bad_PRECIP", "bad_WS", "bad_PM25") := NULL]
  for (v in WEATHER_VARIABLES) {
    d[, (paste0(v, "_night_mean")) := mean(get(paste0(v, "_hourly"))), by = night_id]
    d[, (paste0(v, "_hour_dev")) := get(paste0(v, "_hourly")) - get(paste0(v, "_night_mean"))]
  }
  setorder(d, night_id, hour_position)
  tr_write_csv(d, PANEL_FILE)
  meta <- list(settings = tr_weather_settings(), input_md5 = tools::md5sum(paths),
    panel_md5 = unname(tools::md5sum(PANEL_FILE)), created = Sys.time(),
    n_original_nights = nrow(audit), n_retained_nights = uniqueN(d$night_id))
  saveRDS(meta, PANEL_META_FILE)
  tr_write_csv(tr_file_manifest(paths), file.path(dirname(PANEL_FILE), "input_manifest.csv"))
  cat("保留夜晚：", uniqueN(d$night_id), "/", nrow(audit), "；行数：", nrow(d),
    "；城市数：", uniqueN(d$city), "；夜间日期：", as.character(min(d$DATE)), "至", as.character(max(d$DATE)),
    "\n耗时秒：", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2), "\n")
  invisible(d)
}

tr_weather_read_panel <- function() {
  if (!file.exists(PANEL_FILE) || !file.exists(PANEL_META_FILE)) stop("先运行 00_Timing_prepare_panel.R")
  meta <- readRDS(PANEL_META_FILE)
  if (!identical(meta$settings, tr_weather_settings()) ||
      !identical(meta$panel_md5, unname(tools::md5sum(PANEL_FILE))) ||
      !identical(meta$input_md5, tools::md5sum(names(meta$input_md5))))
    stop("输入或配置发生变化，请重新运行 00_Timing_prepare_panel.R，并重跑惩罚/模型/推断。")
  fread(PANEL_FILE)
}

tr_weather_signature <- function() list(version = TIMING_MODEL_VERSION,
  panel_md5 = unname(tools::md5sum(PANEL_FILE)), settings = tr_weather_settings())
tr_weather_check_signature <- function(x) {
  if (!identical(x$weather_signature, tr_weather_signature())) stop("产物与当前小时面板/控制规格不一致，请重跑上游步骤。")
}

tr_weather_complete_nights <- function(d, keep, OC) {
  d[, valid_hour := keep & !is.na(keep)]
  audit <- d[, .(n_hours = .N, valid_hours = sum(valid_hour)), by = .(city, DATE, night_id)]
  audit[, retained := n_hours == 12L & valid_hours == 12L]
  tr_write_csv(audit, file.path(TIMING_ROOT, "data", paste0("analysis_sample_", OC$outcome, ".csv")))
  d <- d[night_id %in% audit[retained == TRUE, night_id]]
  d[, valid_hour := NULL]
  if (!nrow(d)) stop("结局与温湿度筛选后无完整夜晚")
  check <- d[, .(n = .N, nh = uniqueN(hour_position), n_total = uniqueN(get(OC$night)),
    total = get(OC$night)[1L], sum_y = sum(get(OC$hour))), by = night_id]
  if (any(check$n != 12L | check$nh != 12L | check$n_total != 1L |
    abs(check$total - check$sum_y) > 1e-7 * pmax(1, abs(check$total)))) stop("夜间总量与12小时结局之和不一致")
  d[, DATE := tr_parse_date(DATE)]
  d
}

tr_weather_design <- function(d, out) {
  cols <- character(); scaling <- list()
  for (v in WEATHER_VARIABLES) {
    kinds <- if (CONTROL_MODE == "mean_deviation_hour") c("night_mean", "hour_dev") else "hourly"
    for (kind in kinds) {
      field <- paste0(v, "_", kind)
      z <- d[[field]]
      if (any(!is.finite(z))) stop("气象值缺失：", field)
      # mean列按每夜一次计算尺度；偏离列不平移，保留物理含义。
      sample_z <- if (kind == "night_mean") z[!duplicated(d$night_id)] else z
      center <- if (kind == "hour_dev") 0 else mean(sample_z)
      scale <- sd(sample_z); if (!is.finite(scale) || scale == 0) scale <- 1
      zz <- (z - center) / scale
      hours <- if (kind == "night_mean") WEATHER_HOURS[-1L] else WEATHER_HOURS
      for (h in hours) {
        nm <- paste0("weather_", v, "_", kind, "_hour_", h)
        set(d, j = nm, value = zz * as.numeric(d$hour == h)); cols <- c(cols, nm)
      }
      scaling[[length(scaling) + 1L]] <- data.table(variable = v, component = kind,
        center = center, scale = scale, unit = switch(v, PRECIP = "mm_per_hour", WS = "m_per_s", PM25 = "ug_per_m3"))
    }
  }
  tr_write_csv(rbindlist(scaling), file.path(out, "weather_scaling.csv"))
  cat("气象控制规格：", CONTROL_MODE, "；候选列数：", length(cols), "\n")
  cols
}

# 探索性Wald检验：模型Hessian协方差或日期块bootstrap经验协方差。
# 小时偏离列的检验对应夜内小时变化；均值列对应整夜水平的时段分配关系。
tr_weather_tests <- function(beta, covariance, cn, outcome, method, out) {
  ii <- which(grepl("^weather_", cn)); V <- as.matrix(covariance[ii, ii, drop = FALSE])
  if (!length(ii)) stop("设计矩阵没有保留气象控制列")
  se <- sqrt(pmax(0, diag(V))); z <- ifelse(se > 0, beta[ii] / se, NA_real_)
  tab <- data.table(outcome = outcome, term = cn[ii], estimate_standardized = beta[ii],
    SE = se, z = z, p_value = 2 * pnorm(-abs(z)), CI_lo = beta[ii] - 1.96 * se,
    CI_hi = beta[ii] + 1.96 * se, method = method)
  tab[, p_BH := p.adjust(p_value, method = "BH")]
  tab[, variable := sub("^weather_([^_]+)_.*", "\\1", term)]
  tab[, component := sub("^weather_[^_]+_(.*)_hour_[0-9]+$", "\\1", term)]
  tab[, hour_position := as.integer(sub(".*_hour_", "", term))]
  joint <- rbindlist(lapply(WEATHER_VARIABLES, function(v) {
    components <- c("all", unique(tab[variable == v, component]))
    rbindlist(lapply(components, function(comp) {
      jj <- which(tab$variable == v & (comp == "all" | tab$component == comp))
      if (!length(jj)) return(NULL)
      vv <- V[jj, jj, drop = FALSE]; eig <- eigen((vv + t(vv)) / 2, symmetric = TRUE)
      ok <- eig$values > max(abs(eig$values)) * 1e-8 & eig$values > 0
      df <- sum(ok); b <- beta[ii[jj]]
      stat <- if (df == length(jj) && df > 0L) sum(as.numeric(crossprod(eig$vectors[, ok, drop = FALSE], b))^2 / eig$values[ok]) else NA_real_
      data.table(outcome = outcome, variable = v, component = comp, n_terms = length(jj),
        df = df, Wald = stat, p_value = if (is.finite(stat)) pchisq(stat, df, lower.tail = FALSE) else NA_real_,
        status = if (df == length(jj)) "ok" else "singular_covariance_no_joint_test", method = method)
    }))
  }))
  tr_write_csv(tab, file.path(out, paste0("weather_coefficients_", outcome, "_", method, ".csv")))
  tr_write_csv(joint, file.path(out, paste0("weather_joint_tests_", outcome, "_", method, ".csv")))
  invisible(list(coefficients = tab, joint = joint))
}
