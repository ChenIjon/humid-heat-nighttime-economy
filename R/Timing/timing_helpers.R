# Paper: Heat and Humidity Reshape Urban Economic Activity at Night
# Contact: chenchen1107@connect.hku.hk
# Shared paths and model functions.

suppressPackageStartupMessages({library(data.table); library(mgcv); library(splines); library(Matrix); library(parallel)})
TIMING_ROOT <- normalizePath(.code_dir, winslash = "/", mustWork = TRUE)
TIMING_DATA_DIR <- file.path(TIMING_ROOT, "data")
TIMING_OUTPUT_DIR <- file.path(TIMING_ROOT, "output_timing_model")

TR_CALENDAR_LEVELS <- c(
  "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday",
  "Sunday", "Statutory holiday", "Make-up workday"
)

tr_env_int <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (!is.finite(value) || value < 1L) stop(name, " must be a positive integer", call. = FALSE)
  value
}

tr_env_flag <- function(name, default = FALSE) {
  fallback <- if (isTRUE(default)) "1" else "0"
  tolower(Sys.getenv(name, unset = fallback)) %in% c("1", "true", "yes", "y")
}

tr_outcome_spec <- function(outcome = NULL) {
  if (is.null(outcome)) outcome <- Sys.getenv("EXPOSURE_OUTCOME", unset = "cnt")
  outcome <- tolower(trimws(outcome))
  specs <- list(
    cnt = list(outcome = "cnt", label = "cnt", hour = "hour_cnt", night = "night_total_cnt",
               lambda_file = "laplace_lambda.rds", audit_run_mode = "conditional_city_quarter",
               model_file = "conditional_city_quarter.rds"),
    amt = list(outcome = "amt", label = "amt", hour = "hour_amt", night = "night_total_amt",
               lambda_file = "laplace_lambda_amt.rds", audit_run_mode = "conditional_city_quarter_amt",
               model_file = "conditional_city_quarter_amt.rds")
  )
  if (!outcome %in% names(specs)) stop("EXPOSURE_OUTCOME must be 'cnt' or 'amt', got: ", outcome, call. = FALSE)
  s <- specs[[outcome]]
  s$suffix <- if (identical(outcome, "cnt")) "" else paste0("_", s$label)
  s$out <- function(stem, ext) paste0(stem, s$suffix, ".", ext)
  s
}

tr_calendar_block_map <- function(dates, block_length, seed, replicate = NA_integer_) {
  dates <- sort(unique(as.Date(dates)))
  n_days <- length(dates)
  n_blocks <- ceiling(n_days / block_length)
  set.seed(seed)
  starts <- sample.int(n_days, n_blocks, replace = TRUE)
  selected <- unlist(lapply(starts, function(x) {
    ((x - 1L + seq_len(block_length) - 1L) %% n_days) + 1L
  }))[seq_len(n_days)]
  data.table(
    replicate = replicate,
    boot_day = seq_len(n_days),
    source_DATE = dates[selected],
    boot_block = rep(seq_len(n_blocks), each = block_length)[seq_len(n_days)]
  )
}

tr_parallel_lapply <- function(x, fun, n_cores) {
  if (.Platform$OS.type == "unix" && n_cores > 1L)
    parallel::mclapply(x, fun, mc.cores = min(n_cores, length(x)))
  else lapply(x, fun)
}

tr_reduce_frozen_design <- function(X, cn, cache_file, city_hour_baselines = character()) {
  baseline_drop <- integer()
  if (length(city_hour_baselines)) {
    if (anyDuplicated(city_hour_baselines) || !("(Intercept)" %in% cn))
      stop("City-hour baselines require unique levels and an intercept")
    baseline_names <- paste0("factor(city_hour_id)", city_hour_baselines)
    if (sum(!(baseline_names %in% cn)) != 1L)
      stop("Exactly one city-hour baseline must be the factor reference level")
    baseline_drop <- match(c("(Intercept)", baseline_names), cn, nomatch = 0L)
    baseline_drop <- baseline_drop[baseline_drop > 0L]
  }
  active <- setdiff(seq_along(cn), baseline_drop)
  re <- which(grepl("^s\\(city_part_id\\):(TEMP|RH)", cn))
  if (!length(re)) stop("City-part random-slope columns are missing")
  fixed <- setdiff(active, re)
  rk_fixed <- as.integer(Matrix::rankMatrix(X[, fixed, drop = FALSE], method = "qr"))
  cand <- intersect(active, which(grepl("^(calendar_hour_|sunset_dev_hour_)", cn)))
  protected <- cn[setdiff(active, cand)]
  need <- length(fixed) - rk_fixed
  if (file.exists(cache_file)) {
    k <- readRDS(cache_file)
    fixed_k <- match(intersect(k, cn[fixed]), cn)
    if (length(k) == length(active) - need && !anyDuplicated(k) && all(k %in% cn) &&
        !any(cn[baseline_drop] %in% k) && all(protected %in% k) &&
        length(fixed_k) == rk_fixed &&
        as.integer(Matrix::rankMatrix(X[, fixed_k, drop = FALSE], method = "qr")) == rk_fixed) return(k)
    message("reduced design: cached columns do not span the current design; recomputing")
  }
  drop <- integer(0)
  for (c in cand) { if (length(drop) >= need) break
    cols <- setdiff(fixed, c(drop, c))
    if (as.integer(Matrix::rankMatrix(X[, cols, drop = FALSE], method = "qr")) == rk_fixed) drop <- c(drop, c) }
  if (length(drop) < need)
    stop("Cannot obtain full rank in the unpenalized block by dropping control columns only")
  keep <- cn[setdiff(active, drop)]
  kept_fixed <- setdiff(fixed, drop)
  stopifnot(as.integer(Matrix::rankMatrix(X[, kept_fixed, drop = FALSE], method = "qr")) == length(kept_fixed))
  saveRDS(keep, cache_file)
  message("reduced design: dropped ", length(baseline_drop), " conditional-invariant baselines and ",
          length(drop), " redundant controls -> ", length(keep),
          " columns (unpenalized block full rank; penalized RE overlap allowed)")
  keep
}
