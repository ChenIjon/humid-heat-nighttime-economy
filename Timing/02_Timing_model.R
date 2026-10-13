#!/usr/bin/env Rscript
# Paper: Heat and Humidity Reshape Urban Economic Activity at Night
# The main analyses were conducted using Python (version 3.12.9) and R (version 4.4.1), with the primary statistical models implemented using the mgcv package (version 1.9-3).
# Contact: chenchen1107@connect.hku.hk
# Fit the model using selected penalties; output: model RDS and point-estimate CSVs.
# Input: complete hourly panel with PM2.5, precipitation and wind controls.
# nohup Rscript 02_Timing_model.R > 02_Timing_model.log 2>&1 &
# ps -ef | grep 02_Timing_model.R
# pkill -f 02_Timing_model.R

.code_dir <- Sys.getenv("TIMING_REVISION_DIR", unset = "")
if (!nzchar(.code_dir)) {
  .script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  .code_dir <- if (length(.script_arg)) dirname(normalizePath(sub("^--file=", "", .script_arg[1]), winslash = "/")) else getwd()
}
source(file.path(.code_dir, "timing_helpers.R"))
source(file.path(.code_dir, "hourly_weather_controls.R"))

OUT <- TIMING_OUTPUT_DIR; dir.create(OUT, recursive = TRUE, showWarnings = FALSE); dir.create(file.path(OUT, "models"), recursive = TRUE, showWarnings = FALSE)
SPLINE_DF <- 3L; N_CORES <- tr_env_int("N_CORES", 20L); setDTthreads(N_CORES)
OC <- tr_outcome_spec(); cat("outcome:", OC$outcome, "(hour =", OC$hour, ", night total =", OC$night, ")\n")
LAMBDA_FILE <- file.path(OUT, OC$lambda_file); MODEL_FILE <- file.path(OUT, "models", OC$model_file)
lamL <- readRDS(LAMBDA_FILE); lamT <- lamL$lam[1]; lamR <- lamL$lam[2]
if (!identical(lamL$objective_version, TIMING_MODEL_VERSION))
  stop(basename(LAMBDA_FILE), " does not match the city-part model; rerun 01_Timing_penalty.R on the server first")
if (!is.null(lamL$outcome) && !identical(lamL$outcome, OC$outcome))
  stop(basename(LAMBDA_FILE), " was tuned for outcome '", lamL$outcome, "', not '", OC$outcome, "'")
tr_weather_check_signature(lamL)
cat(sprintf("FROZEN penalties: lam_T = %.6g  lam_R = %.6g\n", lamT, lamR))

d <- tr_weather_read_panel()
d[, night_id := factor(night_id)]; d[, city_hour_id := factor(city_hour_id)]
d[, calendar_type := factor(calendar_type, levels = TR_CALENDAR_LEVELS)]; d[, hour := as.numeric(hour_position)]
d[, sunset_dev_min := sunset_minutes_CST - mean(sunset_minutes_CST), by = city]
keep <- is.finite(d[[OC$hour]]) & d[[OC$hour]] >= 0 & is.finite(d$TEMP) & is.finite(d$RH) &
  is.finite(d$sunset_dev_min) & d$RH >= 0 & d$RH <= 100 & d[[OC$night]] > 0
d <- tr_weather_complete_nights(d, keep, OC); d[, hour_factor := factor(hour, levels = 0:11)]
d[, night_part := factor(cut(hour, breaks = c(-1, 2, 5, 8, 11),
  labels = c("h00_02", "h03_05", "h06_08", "h09_11")))]
d[, city_part_id := interaction(city, night_part, drop = TRUE, sep = ".")]
setorder(d, night_id, hour); d[, grp := .GRP, by = night_id]
cat("input:", PANEL_FILE, "\noutput:", OUT,
    "\nrows:", nrow(d), " nights:", uniqueN(d$night_id), " city-part groups:", uniqueN(d$city_part_id),
    "\noverwrite: existing outputs with the same names\n")
CH <- setdiff(levels(d$hour_factor), "0"); CC <- levels(d$calendar_type)
SC <- paste0("sunset_dev_hour_", seq_along(CH)); CCal <- unlist(lapply(seq_along(CH), function(i) paste0("calendar_hour_", i, "_", seq_along(CC))))
TS <- unlist(lapply(seq_len(SPLINE_DF), function(j) paste0("TEMP_ns", j, "_hour_", seq_along(CH)))); RS <- unlist(lapply(seq_len(SPLINE_DF), function(j) paste0("RH_ns", j, "_hour_", seq_along(CH))))
hv <- as.character(d$hour_factor); cv <- as.character(d$calendar_type)
for (i in seq_along(CH)) { in_h <- hv == CH[i]; d[, (SC[i]) := sunset_dev_min * as.numeric(in_h)]
  for (j in seq_along(CC)) d[, (paste0("calendar_hour_", i, "_", j)) := as.numeric(in_h & cv == CC[j])] }
bT <- splines::ns(d$TEMP, df = SPLINE_DF, intercept = FALSE); bR <- splines::ns(d$RH, df = SPLINE_DF, intercept = FALSE)
d[, paste0("TEMP_ns", seq_len(SPLINE_DF)) := as.data.table(unclass(bT))]; d[, paste0("RH_ns", seq_len(SPLINE_DF)) := as.data.table(unclass(bR))]
for (j in seq_len(SPLINE_DF)) for (i in seq_along(CH)) { in_h <- as.numeric(hv == CH[i])
  d[, (paste0("TEMP_ns", j, "_hour_", i)) := get(paste0("TEMP_ns", j)) * in_h]
  d[, (paste0("RH_ns", j, "_hour_", i)) := get(paste0("RH_ns", j)) * in_h] }
WX <- tr_weather_design(d, OUT)
nl <- c(TS, RS); ct <- c(SC, CCal, WX); d[, alpha := log(get(OC$night))]
form <- as.formula(paste(OC$hour, "~", paste(c(nl, ct, "factor(city_hour_id)",
  "s(city_part_id, by = TEMP, bs = 're')", "s(city_part_id, by = RH, bs = 're')"), collapse = " + "), "+ offset(alpha)"))
fit0 <- bam(form, data = d, family = poisson(link = "log"), method = "fREML", discrete = TRUE, nthreads = N_CORES, control = gam.control(maxit = 200))
X <- as(mgcv::predict.gam(fit0, newdata = as.data.frame(d), type = "lpmatrix"), "dgCMatrix"); FULL_CN <- colnames(X)
keep_cols <- tr_reduce_frozen_design(X, FULL_CN, file.path(OUT, "reduced_design_columns.rds"),
  as.character(unique(d[hour == 0, city_hour_id])), night_groups = d$grp)
X <- X[, match(keep_cols, FULL_CN)]; cn <- keep_cols
if (!identical(cn, lamL$cn)) stop("Lambda file and model design columns differ; rerun 01_Timing_penalty.R")
iT <- which(grepl("s\\(city_part_id\\):TEMP", cn)); iR <- which(grepl("s\\(city_part_id\\):RH", cn)); iU <- c(iT, iR)
if (length(iT) != uniqueN(d$city_part_id) || length(iR) != uniqueN(d$city_part_id))
  stop("City-part slope columns do not match the observed city-part groups")
Sd_fixed <- numeric(ncol(X)); Sd_fixed[iT] <- lamT; Sd_fixed[iR] <- lamR
y <- d[[OC$hour]]; Ntot <- d[[OC$night]]; hour_vec <- d$hour
RIDGE_REL <- 1e-10; GRAD_TOL <- 1e-6; NEWTON_GAIN_TOL <- 1e-6
INNER_MAXIT <- tr_env_int("INNER_MAXIT", 25L)
if (INNER_MAXIT < 1L) stop("INNER_MAXIT must be positive")
cat("design p =", ncol(X), " q_T =", length(iT), " q_R =", length(iR), "\n")

gain_tol_at <- function(nll) max(NEWTON_GAIN_TOL, .Machine$double.eps * max(1, abs(nll)))

solve_cond <- function(Xb, yb, Nb, gb, hb, Sd) { p <- ncol(Xb); beta <- numeric(p)
  gstart <- { dt <- data.table(g = gb, r = seq_along(gb)); dt[, .(s = r[1L]), by = g][order(g)]$s }
  pnf <- function(eta) { dt <- data.table(grp = gb, eta = eta); dt[, mx := max(eta), by = grp]; dt[, lse := mx[1L] + log(sum(exp(eta - mx[1L]))), by = grp]; dt[, p := exp(eta - lse)]; dt[, N := Nb]; dt }
  og <- function(b, wH = TRUE) { eta <- as.numeric(Xb %*% b); pn <- pnf(eta); lse <- pn$lse; pp <- pn$p
    nll <- sum(pn[, .(v = N[1L] * lse[1L]), by = grp]$v) - sum(yb * eta) + 0.5 * sum(Sd * b^2); grad <- as.numeric(crossprod(Xb, Nb * pp - yb)) + Sd * b; H <- NULL
    if (wH) { ii <- rep(seq_along(gb), times = 12L); jj <- integer(length(gb) * 12L); xx <- numeric(length(gb) * 12L); pos <- 0L
      for (bb in 0:11) { colb <- gstart[gb] + bb; val <- Nb * pp * (as.numeric(hb == bb) - pp[colb]); jj[(pos + 1L):(pos + length(gb))] <- colb; xx[(pos + 1L):(pos + length(gb))] <- val; pos <- pos + length(gb) }
      H <- crossprod(Xb, sparseMatrix(i = ii, j = jj, x = xx, dims = c(length(gb), length(gb))) %*% Xb) + Diagonal(x = Sd)
      H <- as(H, "symmetricMatrix"); hdiag <- Matrix::diag(H)
      if (any(!is.finite(hdiag)) || any(hdiag <= 0))
        stop("Conditional Hessian has a nonpositive or nonfinite diagonal")
      diag(H) <- hdiag + RIDGE_REL * hdiag }
    list(nll = nll, grad = grad, H = H) }
  cur <- og(beta); ls_stop <- FALSE; it <- 0L; g0 <- max(abs(cur$grad)); last_rel <- NA_real_
  for (it in seq_len(INNER_MAXIT)) { step <- as.numeric(solve(cur$H, -cur$grad))
    gain <- -0.5 * sum(cur$grad * step)
    if (max(abs(cur$grad)) <= GRAD_TOL || (is.finite(gain) && gain >= 0 && gain <= gain_tol_at(cur$nll))) break
    if (!is.finite(gain) || gain < 0) { ls_stop <- TRUE; break }
    tt <- 1; ok <- FALSE
    while (tt > 1e-6) { cd <- beta + tt * step; if (og(cd, FALSE)$nll < cur$nll - 1e-12) { ok <- TRUE; break }; tt <- tt / 2 }
    if (!ok) { ls_stop <- TRUE; break }
    prev_nll <- cur$nll; beta <- beta + tt * step; cur <- og(beta); last_rel <- (prev_nll - cur$nll) / (1 + abs(prev_nll))
    if (max(abs(cur$grad)) < 1e-6) break }
  rel_g <- max(abs(cur$grad)) / max(g0, .Machine$double.xmin)
  step <- as.numeric(solve(cur$H, -cur$grad)); gain <- -0.5 * sum(cur$grad * step)
  list(beta = beta, nll = cur$nll, ginf = max(abs(cur$grad)), it = it, ls_stop = ls_stop,
       g0 = g0, rel_g = rel_g, newton_gain = gain, last_rel = last_rel, H = cur$H) }
classify_solver <- function(r) { if (!is.finite(r$nll) || !is.finite(r$ginf) || !is.finite(r$newton_gain)) return("failed_nonfinite")
  if (r$ginf <= GRAD_TOL || (r$newton_gain >= 0 && r$newton_gain <= gain_tol_at(r$nll))) return("converged_gradient")
  if (r$ls_stop) return("failed_line_search")
  "failed_maxit" }
is_converged <- function(s) s == "converged_gradient"

block2 <- function(nm) { ifelse(grepl("^TEMP_ns[0-9]+_hour_", nm), "shared_T", ifelse(grepl("^RH_ns[0-9]+_hour_", nm), "shared_R",
  ifelse(grepl("s\\(city_part_id\\):TEMP", nm), "dev_T", ifelse(grepl("s\\(city_part_id\\):RH", nm), "dev_R", "ctrl")))) }
B2 <- block2(cn)
make_nd <- function(ex, val) { nd <- copy(d); bb <- if (ex == "TEMP") bT else bR; nb <- predict(bb, val)
  nd[[ex]] <- val
  for (j in seq_len(SPLINE_DF)) { nd[[paste0(ex, "_ns", j)]] <- nb[, j]
    for (i in seq_along(CH)) nd[[paste0(ex, "_ns", j, "_hour_", i)]] <- nb[, j] * (d$hour == as.numeric(CH[i])) }; nd }
frozen_X <- function(nd, fam) { Xn <- as(mgcv::predict.gam(fit0, newdata = as.data.frame(nd), type = "lpmatrix"), "dgCMatrix")
  Xc <- X; fc <- if (fam == "TEMP") which(B2 %in% c("shared_T", "dev_T")) else which(B2 %in% c("shared_R", "dev_R"))
  Xc[, fc] <- Xn[, match(cn[fc], colnames(Xn))]   # invariant: Xc[-focal] == X[-focal] holds by construction
  Xc }
wshare <- function(eta, g, p, n_part = 4L) { dt <- data.table(g = g, e = eta, p = p)
  dt[, mx := max(e), by = g]; dt[, lse := mx[1L] + log(sum(exp(e - mx[1L]))), by = g]; dt[, pp := exp(e - lse)]
  agg <- dt[, .(s = sum(pp)), by = .(g, p)]; ng <- uniqueN(g)
  M <- matrix(0, ng, n_part); M[cbind(agg$g, agg$p)] <- agg$s; M }
hprof <- function(eta, g, hb) { dt <- data.table(g = g, e = eta, h = hb); dt[, mx := max(e), by = g]
  dt[, lse := mx[1L] + log(sum(exp(e - mx[1L]))), by = g]; dt[, pp := exp(e - lse)]; dt[, .(v = mean(pp)), by = h][order(h)]$v }

full <- solve_cond(X, y, Ntot, d$grp, hour_vec, Sd_fixed); be <- full$beta
full_status <- classify_solver(full)
cat(sprintf("conditional fit: iters=%d ginf=%.3e newton_gain=%.3e status=%s\n", full$it, full$ginf, full$newton_gain, full_status))
if (!is_converged(full_status)) stop("Conditional fit did not converge: ", full_status, "; grad_inf=", signif(full$ginf, 5))
partv <- as.integer(d$night_part); PART_LEVELS <- levels(d$night_part)
qT <- d[, .(lo = as.numeric(quantile(TEMP, .025)), hi = as.numeric(quantile(TEMP, .975))), by = city]
qR <- d[, .(lo = as.numeric(quantile(RH, .025)), hi = as.numeric(quantile(RH, .975))), by = city]
n <- nrow(d); qTp <- quantile(d$TEMP, c(.025, .975)); qRp <- quantile(d$RH, c(.025, .975))
CF <- list(Tp_lo = frozen_X(make_nd("TEMP", rep(qTp[1], n)), "TEMP"), Tp_hi = frozen_X(make_nd("TEMP", rep(qTp[2], n)), "TEMP"),
  Rp_lo = frozen_X(make_nd("RH", rep(qRp[1], n)), "RH"), Rp_hi = frozen_X(make_nd("RH", rep(qRp[2], n)), "RH"),
  Tc_lo = frozen_X(make_nd("TEMP", qT$lo[match(d$city, qT$city)]), "TEMP"), Tc_hi = frozen_X(make_nd("TEMP", qT$hi[match(d$city, qT$city)]), "TEMP"),
  Rc_lo = frozen_X(make_nd("RH", qR$lo[match(d$city, qR$city)]), "RH"), Rc_hi = frozen_X(make_nd("RH", qR$hi[match(d$city, qR$city)]), "RH"))
cityb <- as.character(unique(d[, .(grp, city)])[order(grp)]$city); cidx <- split(seq_along(cityb), cityb)
eta_cf <- function(nm) as.numeric(CF[[nm]] %*% be)
pooled <- rbindlist(lapply(c("TEMP", "RH"), function(ex) { a <- if (ex == "TEMP") "Tp_lo" else "Rp_lo"; z <- if (ex == "TEMP") "Tp_hi" else "Rp_hi"
  dch <- 100 * (colMeans(wshare(eta_cf(z), d$grp, partv)) - colMeans(wshare(eta_cf(a), d$grp, partv)))
  data.table(model = "nonlinear_city_part", exposure = ex, part = PART_LEVELS, change_pp = dch) }))
citychange <- rbindlist(lapply(c("TEMP", "RH"), function(ex) { a <- if (ex == "TEMP") "Tc_lo" else "Rc_lo"; z <- if (ex == "TEMP") "Tc_hi" else "Rc_hi"
  sl <- wshare(eta_cf(a), d$grp, partv); sh <- wshare(eta_cf(z), d$grp, partv)
  rbindlist(lapply(names(cidx), function(ct) { j <- cidx[[ct]]
    data.table(model = "nonlinear_city_part", exposure = ex, city = ct, part = PART_LEVELS,
      share_low = 100 * colMeans(sl[j, , drop = FALSE]),
      share_high = 100 * colMeans(sh[j, , drop = FALSE]),
      change_pp = 100 * (colMeans(sh[j, , drop = FALSE]) - colMeans(sl[j, , drop = FALSE]))) })) }))
hourprof <- rbindlist(lapply(c("TEMP", "RH"), function(ex) { a <- if (ex == "TEMP") "Tp_lo" else "Rp_lo"; z <- if (ex == "TEMP") "Tp_hi" else "Rp_hi"
  data.table(model = "nonlinear_city_part", exposure = ex, hour = 0:11,
    profile_low = 100 * hprof(eta_cf(a), d$grp, hour_vec), profile_high = 100 * hprof(eta_cf(z), d$grp, hour_vec),
    change_pp = 100 * (hprof(eta_cf(z), d$grp, hour_vec) - hprof(eta_cf(a), d$grp, hour_vec))) }))

GRID_N <- tr_env_int("GRID_N", 21L)
if (GRID_N < 2L) stop("GRID_N must be >= 2")
grid_pct <- seq(0.025, 0.975, length.out = GRID_N)
grid_share <- function(ex, val) hprof(as.numeric(frozen_X(make_nd(ex, rep(val, n)), ex) %*% be), d$grp, hour_vec)
hour_surface <- rbindlist(lapply(c("TEMP", "RH"), function(ex) {
  vals <- as.numeric(quantile(d[[ex]], grid_pct))
  ref_share <- grid_share(ex, as.numeric(quantile(d[[ex]], 0.025)))
  rbindlist(lapply(seq_along(grid_pct), function(k) {
    sh <- grid_share(ex, vals[k])
    data.table(model = "nonlinear_city_part", exposure = ex, percentile = grid_pct[k],
      value = vals[k], hour = 0:11, share_pct = 100 * sh, dshare_pp = 100 * (sh - ref_share)) })) }))

Vbeta <- solve(full$H)
tr_weather_tests(be, Vbeta, cn, OC$outcome, "hessian", OUT)
Hinv_diag <- Matrix::diag(Vbeta)
edf_T <- length(iT) - sum(Hinv_diag[iT] * Sd_fixed[iT]); edf_R <- length(iR) - sum(Hinv_diag[iR] * Sd_fixed[iR])
fit_audit <- data.table(run_mode = OC$audit_run_mode, outcome = OC$outcome, lam_T = lamT, lam_R = lamR,
  grad_tol = GRAD_TOL, newton_gain_tol = NEWTON_GAIN_TOL,
  iters = full$it, grad_inf = full$ginf, grad0_inf = full$g0, rel_grad = full$rel_g,
  newton_gain = full$newton_gain, last_step_rel_gain = full$last_rel,
  ls_stop = full$ls_stop, termination_status = full_status, converged = is_converged(full_status), penalized_nll = full$nll,
  weather_mode = CONTROL_MODE, n_weather_candidates = length(WX),
  n_weather_retained = sum(grepl("^weather_", cn)), n_coefficients = length(be), n_nonfinite_coefts = sum(!is.finite(be)),
  edf_TEMP_dev = edf_T, edf_RH_dev = edf_R, ref_df_TEMP_dev = length(iT), ref_df_RH_dev = length(iR),
  beta_rms_TEMP = sqrt(mean(be[iT]^2)), beta_rms_RH = sqrt(mean(be[iR]^2)))
support <- d[, .(n_rows = .N, n_city_nights = uniqueN(night_id)), by = .(city, hour_factor)][order(city, hour_factor)]

fwrite(pooled, file.path(OUT, OC$out("pooled_window_share_change", "csv")))
fwrite(citychange, file.path(OUT, OC$out("city_window_share_change", "csv")))
fwrite(hourprof, file.path(OUT, OC$out("hour_redistribution_profile", "csv")))
fwrite(hour_surface, file.path(OUT, OC$out("hour_exposure_grid", "csv")))
fwrite(fit_audit, file.path(OUT, OC$out("model_fit_audit", "csv")))
fwrite(support, file.path(OUT, OC$out("city_hour_support", "csv")))
saveRDS(list(beta = be, cn = cn, weather_signature = tr_weather_signature(), fit_audit = fit_audit, pooled = pooled, citychange = citychange, hourprof = hourprof),
        MODEL_FILE)
cat("\n=== pooled P2.5->P97.5 per-window share contrast (pp) ===\n"); print(pooled, digits = 5)
cat("\n=== per-city per-window change (pp): summary ===\n")
print(citychange[, .(mean = mean(change_pp), min = min(change_pp), max = max(change_pp), sd = sd(change_pp)), by = .(exposure, part)], digits = 4)
cat("\n=== fit audit ===\n"); print(fit_audit, digits = 5)
message("Saved to ", OUT)
