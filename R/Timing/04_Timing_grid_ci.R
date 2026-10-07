#!/usr/bin/env Rscript
# Paper: Heat and Humidity Reshape Urban Economic Activity at Night
# Contact: chenchen1107@connect.hku.hk
# Compute grid intervals from saved model and bootstrap outputs; output: grid CI CSV.

.code_dir <- Sys.getenv("TIMING_REVISION_DIR", unset = "")
if (!nzchar(.code_dir)) {
  .script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  .code_dir <- if (length(.script_arg)) dirname(normalizePath(sub("^--file=", "", .script_arg[1]), winslash = "/")) else getwd()
}
source(file.path(.code_dir, "timing_helpers.R"))

OUT <- TIMING_OUTPUT_DIR
SPLINE_DF <- 3L; N_CORES <- tr_env_int("N_CORES", 4L)
GRID_CI_BOOT <- tr_env_int("GRID_CI_BOOT", 1000L)
if (N_CORES < 1L || GRID_CI_BOOT < 20L) stop("N_CORES must be positive and GRID_CI_BOOT must be at least 20")
setDTthreads(N_CORES)
OC <- tr_outcome_spec(); cat("outcome:", OC$outcome, "(hour =", OC$hour, ", night total =", OC$night, ")\n")
MODEL_FILE <- file.path(OUT, "models", OC$model_file)
LAMBDA_FILE <- file.path(OUT, OC$lambda_file)
REPL_FILE <- file.path(OUT, OC$out(sprintf("bootstrap_replicates_B%d", GRID_CI_BOOT), "rds"))
GRID_FILE <- file.path(OUT, OC$out("hour_exposure_grid", "csv"))
STATUS_FILE <- file.path(OUT, OC$out("bootstrap_status", "csv"))
SPEC_FILE <- file.path(OUT, OC$out("bootstrap_specification", "csv"))
CI_FILE <- file.path(OUT, OC$out("hour_exposure_grid_ci", "csv"))
cat("input:", file.path(TIMING_DATA_DIR, "hourly_profile_panel.csv"),
    "\nbootstrap:", REPL_FILE,
    "\noutput:", CI_FILE,
    "\noverwrite: output CSV with the same name\n")

tr_grid_eta <- function(beta, idx, ex, val, eta_obs, d, cn, re_idx, spline_basis) {
  if (length(re_idx) != nlevels(d$city_part_id)) stop("Random-slope columns do not match city-part levels")
  hh <- as.integer(d$hour[idx]); city_part <- as.integer(d$city_part_id[idx])
  co <- matrix(0, nrow = 12L, ncol = SPLINE_DF)
  for (j in seq_len(SPLINE_DF)) for (h in 1:11) {
    k <- match(paste0(ex, "_ns", j, "_hour_", h), cn)
    if (!is.na(k)) co[h + 1L, j] <- beta[k]
  }
  old_shared <- numeric(length(idx))
  new_shared <- numeric(length(idx))
  nb <- as.numeric(predict(spline_basis, val))
  for (j in seq_len(SPLINE_DF)) {
    w <- co[hh + 1L, j]
    old_shared <- old_shared + d[[paste0(ex, "_ns", j)]][idx] * w
    new_shared <- new_shared + nb[j] * w
  }
  dev <- beta[re_idx][city_part]
  eta_obs - old_shared - d[[ex]][idx] * dev + new_shared + val * dev
}
tr_grid_profile <- function(eta) {
  if (length(eta) %% 12L != 0L || any(!is.finite(eta))) stop("Invalid 12-hour prediction vector")
  a <- matrix(eta, nrow = 12L)
  z <- exp(sweep(a, 2L, apply(a, 2L, max), "-"))
  rowMeans(sweep(z, 2L, colSums(z), "/"))
}

model <- readRDS(MODEL_FILE)
lam <- readRDS(LAMBDA_FILE)
if (!identical(lam$objective_version, "conditional_city_quarter_laplace_v5") ||
    !identical(model$cn, lam$cn)) stop("Model and lambda must be from the same v5 design")
if (!is.null(lam$outcome) && !identical(lam$outcome, OC$outcome))
  stop(basename(LAMBDA_FILE), " was tuned for outcome '", lam$outcome, "', not '", OC$outcome, "'")
grid <- fread(GRID_FILE)
if (!all(c("exposure", "percentile", "value", "hour", "dshare_pp") %in% names(grid)) ||
    !setequal(grid$exposure, c("TEMP", "RH")) ||
    anyDuplicated(grid[, .(exposure, percentile, hour)])) stop("Invalid model exposure grid")
points <- lapply(c("TEMP", "RH"), function(ex) {
  p <- unique(grid[exposure == ex, .(percentile, value)]); setorder(p, percentile)
  if (nrow(p) < 2L || anyDuplicated(p$value)) stop("Invalid percentile grid for ", ex)
  p
})
names(points) <- c("TEMP", "RH")
if (nrow(points$TEMP) != nrow(points$RH)) stop("TEMP and RH grids have different lengths")

d <- fread(file.path(TIMING_DATA_DIR, "hourly_profile_panel.csv"))
d[, night_id := factor(night_id)]; d[, city_hour_id := factor(city_hour_id)]
d[, calendar_type := factor(calendar_type, levels = TR_CALENDAR_LEVELS)]; d[, hour := as.numeric(hour_position)]
d[, sunset_dev_min := sunset_minutes_CST - mean(sunset_minutes_CST), by = city]
keep <- is.finite(d[[OC$hour]]) & d[[OC$hour]] >= 0 & is.finite(d$TEMP) & is.finite(d$RH) &
  is.finite(d$sunset_dev_min) & d$RH >= 0 & d$RH <= 100 & d[[OC$night]] > 0
d <- d[keep]; d[, hour_factor := factor(hour, levels = 0:11)]
d[, night_part := factor(cut(hour, breaks = c(-1, 2, 5, 8, 11),
  labels = c("h00_02", "h03_05", "h06_08", "h09_11")))]
d[, city_part_id := interaction(city, night_part, drop = TRUE, sep = ".")]
setorder(d, night_id, hour); d[, grp := .GRP, by = night_id]
if (nrow(d) %% 12L != 0L || any(d$hour != rep(0:11, nrow(d) / 12L)))
  stop("The panel is not ordered as complete 12-hour city-nights")
CH <- setdiff(levels(d$hour_factor), "0"); CC <- levels(d$calendar_type)
SC <- paste0("sunset_dev_hour_", seq_along(CH))
CCal <- unlist(lapply(seq_along(CH), function(i) paste0("calendar_hour_", i, "_", seq_along(CC))))
TS <- unlist(lapply(seq_len(SPLINE_DF), function(j) paste0("TEMP_ns", j, "_hour_", seq_along(CH))))
RS <- unlist(lapply(seq_len(SPLINE_DF), function(j) paste0("RH_ns", j, "_hour_", seq_along(CH))))
hv <- as.character(d$hour_factor); cv <- as.character(d$calendar_type)
for (i in seq_along(CH)) {
  in_h <- hv == CH[i]; d[, (SC[i]) := sunset_dev_min * as.numeric(in_h)]
  for (j in seq_along(CC)) d[, (paste0("calendar_hour_", i, "_", j)) := as.numeric(in_h & cv == CC[j])]
}
bT <- splines::ns(d$TEMP, df = SPLINE_DF, intercept = FALSE)
bR <- splines::ns(d$RH, df = SPLINE_DF, intercept = FALSE)
d[, paste0("TEMP_ns", seq_len(SPLINE_DF)) := as.data.table(unclass(bT))]
d[, paste0("RH_ns", seq_len(SPLINE_DF)) := as.data.table(unclass(bR))]
for (j in seq_len(SPLINE_DF)) for (i in seq_along(CH)) {
  in_h <- as.numeric(hv == CH[i])
  d[, (paste0("TEMP_ns", j, "_hour_", i)) := get(paste0("TEMP_ns", j)) * in_h]
  d[, (paste0("RH_ns", j, "_hour_", i)) := get(paste0("RH_ns", j)) * in_h]
}
nl <- c(TS, RS); ct <- c(SC, CCal); d[, alpha := log(get(OC$night))]
form <- as.formula(paste(OC$hour, "~", paste(c(nl, ct, "factor(city_hour_id)",
  "s(city_part_id, by = TEMP, bs = 're')", "s(city_part_id, by = RH, bs = 're')"), collapse = " + "), "+ offset(alpha)"))
fit0 <- bam(form, data = d, family = poisson(link = "log"), method = "fREML",
            discrete = TRUE, nthreads = N_CORES, control = gam.control(maxit = 200))
X <- as(mgcv::predict.gam(fit0, newdata = as.data.frame(d), type = "lpmatrix"), "dgCMatrix")
full_cn <- colnames(X)
keep_cols <- tr_reduce_frozen_design(X, full_cn, file.path(OUT, "reduced_design_columns.rds"),
  as.character(unique(d[hour == 0, city_hour_id])))
X <- X[, match(keep_cols, full_cn), drop = FALSE]; cn <- keep_cols
if (!identical(cn, model$cn)) stop("Rebuilt grid design differs from the saved model")
re_idx <- list(TEMP = which(grepl("s\\(city_part_id\\):TEMP", cn)),
               RH = which(grepl("s\\(city_part_id\\):RH", cn)))
basis <- list(TEMP = bT, RH = bR)
for (ex in names(points)) {
  expected <- as.numeric(quantile(d[[ex]], points[[ex]]$percentile))
  if (max(abs(expected - points[[ex]]$value)) > 1e-7)
    stop("Exposure grid values differ from the reconstructed panel for ", ex)
}
cat("design p =", ncol(X), " nights =", uniqueN(d$night_id),
    " grid points =", nrow(points$TEMP), "per exposure\n")

grid_eval <- function(beta, idx) {
  if (length(beta) != ncol(X) || length(idx) %% 12L != 0L ||
      any(d$hour[idx] != rep(0:11, length(idx) / 12L)))
    stop("Bootstrap beta or sampled 12-hour ordering differs from the model")
  eta_obs <- as.numeric(X[idx, , drop = FALSE] %*% beta)
  out <- array(NA_real_, dim = c(2L, nrow(points$TEMP), 12L))
  for (ee in seq_along(points)) {
    ex <- names(points)[ee]; pp <- points[[ex]]
    ref_eta <- tr_grid_eta(beta, idx, ex, pp$value[1L], eta_obs, d, cn, re_idx[[ex]], basis[[ex]])
    ref <- tr_grid_profile(ref_eta)
    for (k in seq_len(nrow(pp))) {
      eta <- tr_grid_eta(beta, idx, ex, pp$value[k], eta_obs, d, cn, re_idx[[ex]], basis[[ex]])
      out[ee, k, ] <- 100 * (tr_grid_profile(eta) - ref)
    }
  }
  out
}

point <- grid_eval(model$beta, seq_len(nrow(d)))
for (ee in seq_along(points)) for (k in seq_len(nrow(points[[ee]]))) {
  ex <- names(points)[ee]; pp <- points[[ee]]$percentile[k]
  target <- grid[exposure == ex & percentile == pp][order(hour)]$dshare_pp
  if (length(target) != 12L || max(abs(point[ee, k, ] - target)) > 1e-5)
    stop("Fast grid calculation does not match the model stage at ", ex, ", percentile ", pp)
}
cat("point-estimate grid matches the model stage\n")

seeds <- readRDS(file.path(OUT, "bootstrap_seeds.rds"))
res <- readRDS(REPL_FILE)
status <- fread(STATUS_FILE)
spec <- fread(SPEC_FILE)
if (nrow(spec) != 1L || spec$n_boot != GRID_CI_BOOT ||
    !is.finite(spec$block_length_days) || spec$block_length_days < 1L)
  stop("Bootstrap specification does not match GRID_CI_BOOT")
BLOCK_LENGTH <- as.integer(spec$block_length_days)
if (length(seeds) != GRID_CI_BOOT || length(res) != GRID_CI_BOOT ||
    nrow(status) != GRID_CI_BOOT || !identical(status$replicate, seq_len(GRID_CI_BOOT)))
  stop("Bootstrap seeds, coefficient RDS and status table do not match GRID_CI_BOOT")
good <- which(status$ok & status$converged)
if (length(good) < 20L || any(vapply(good, function(i) {
    r <- res[[i]]
    r$i != i || length(r$beta) != ncol(X) || !all(is.finite(r$beta)) || !isTRUE(r$ok)
  }, logical(1))))
  stop("Too few usable bootstrap coefficients or beta length differs from the model")
cat("usable bootstrap replicates:", length(good), "/", GRID_CI_BOOT, "\n")
d2 <- d[, .(DATE, city, hour)]
bootstrap_index <- function(i) {
  map <- tr_calendar_block_map(d$DATE, BLOCK_LENGTH, seeds[i])
  mrg <- merge(map[, .(boot_day, source_DATE)],
               d2[, .(row = .I, DATE, city, hour)],
               by.x = "source_DATE", by.y = "DATE", allow.cartesian = TRUE, sort = FALSE)
  setorder(mrg, city, boot_day, hour)
  mrg$row
}
rep_grid <- function(i) {
  setDTthreads(1)
  grid_eval(res[[i]]$beta, bootstrap_index(i))
}
first <- rep_grid(good[1L])
last_k <- nrow(points$TEMP)
if (length(res[[good[1L]]]$prof_T) != 12L || length(res[[good[1L]]]$prof_R) != 12L ||
    max(abs(first[1L, last_k, ] - res[[good[1L]]]$prof_T)) > 1e-5 ||
    max(abs(first[2L, last_k, ] - res[[good[1L]]]$prof_R)) > 1e-5)
  stop("Bootstrap endpoint does not match the saved hourly-profile replicate")
cat("bootstrap endpoint matches the saved hourly-profile replicate\n")
surfaces <- vector("list", length(good))
chunk <- 50L
for (s in seq(1L, length(good), by = chunk)) {
  ix <- s:min(s + chunk - 1L, length(good))
  surfaces[ix] <- tr_parallel_lapply(good[ix], rep_grid, N_CORES)
  if (any(vapply(surfaces[ix], inherits, logical(1), "try-error")))
    stop("A grid-prediction worker failed; inspect the log before rerunning")
  cat("grid intervals:", max(ix), "/", length(good), "replicates done\n")
}
lo <- hi <- array(NA_real_, dim = dim(point))
for (ee in seq_len(dim(point)[1L])) for (k in seq_len(dim(point)[2L])) for (h in 1:12) {
  z <- vapply(surfaces, function(a) a[ee, k, h], numeric(1))
  lo[ee, k, h] <- unname(quantile(z, 0.025)); hi[ee, k, h] <- unname(quantile(z, 0.975))
}
ci <- copy(grid)
for (ee in seq_along(points)) for (k in seq_len(nrow(points[[ee]]))) {
  ex <- names(points)[ee]; pp <- points[[ee]]$percentile[k]
  rows <- which(ci$exposure == ex & ci$percentile == pp)
  rows <- rows[order(ci$hour[rows])]
  ci[rows, `:=`(CI_lo = lo[ee, k, ], CI_hi = hi[ee, k, ], n_success = length(good))]
}
ci[, CI_excludes_zero := CI_lo > 0 | CI_hi < 0]
if (any(!is.finite(ci$CI_lo) | !is.finite(ci$CI_hi) | ci$CI_lo > ci$CI_hi))
  stop("Invalid bootstrap grid intervals")
fwrite(ci, CI_FILE)
message("Saved ", basename(CI_FILE), " with pointwise 95% percentile intervals")
