#!/usr/bin/env Rscript
# Paper: Heat and Humidity Reshape Urban Economic Activity at Night
# Contact: chenchen1107@connect.hku.hk
# Block bootstrap using selected penalties; output: bootstrap RDS and CSVs.

.code_dir <- Sys.getenv("TIMING_REVISION_DIR", unset = "")
if (!nzchar(.code_dir)) {
  .script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  .code_dir <- if (length(.script_arg)) dirname(normalizePath(sub("^--file=", "", .script_arg[1]), winslash = "/")) else getwd()
}
source(file.path(.code_dir, "timing_helpers.R"))

OUT <- TIMING_OUTPUT_DIR; dir.create(OUT, recursive = TRUE, showWarnings = FALSE); SPLINE_DF <- 3L
N_BOOT <- tr_env_int("N_BOOT", 1000L); N_CORES <- tr_env_int("N_CORES", 50L)
BLOCK_LENGTH <- tr_env_int("BLOCK_LENGTH", 7L); SEED0 <- 20260919L; setDTthreads(N_CORES)
`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a
OC <- tr_outcome_spec(); cat("outcome:", OC$outcome, "(hour =", OC$hour, ", night total =", OC$night, ")\n")
LAMBDA_FILE <- file.path(OUT, OC$lambda_file)
lamL <- readRDS(LAMBDA_FILE); lamT <- lamL$lam[1]; lamR <- lamL$lam[2]
if (!identical(lamL$objective_version, "conditional_city_quarter_laplace_v5"))
  stop(basename(LAMBDA_FILE), " does not match the city-part model; rerun the penalty stage first")
if (!is.null(lamL$outcome) && !identical(lamL$outcome, OC$outcome))
  stop(basename(LAMBDA_FILE), " was tuned for outcome '", lamL$outcome, "', not '", OC$outcome, "'")
cat("fixed lambda: lam_T =", signif(lamT, 8), " lam_R =", signif(lamR, 8), "\n")

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
cat("input:", file.path(TIMING_DATA_DIR, "hourly_profile_panel.csv"), "\noutput:", OUT,
    "\nrows:", nrow(d), " nights:", uniqueN(d$night_id), " city-part groups:", uniqueN(d$city_part_id),
    " bootstrap replicates:", N_BOOT, "\noverwrite: existing outputs with the same names\n")
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
nl <- c(TS, RS); ct <- c(SC, CCal); d[, alpha := log(get(OC$night))]
form <- as.formula(paste(OC$hour, "~", paste(c(nl, ct, "factor(city_hour_id)",
  "s(city_part_id, by = TEMP, bs = 're')", "s(city_part_id, by = RH, bs = 're')"), collapse = " + "), "+ offset(alpha)"))
fit0 <- bam(form, data = d, family = poisson(link = "log"), method = "fREML", discrete = TRUE, nthreads = N_CORES, control = gam.control(maxit = 200))
X <- as(mgcv::predict.gam(fit0, newdata = as.data.frame(d), type = "lpmatrix"), "dgCMatrix"); FULL_CN <- colnames(X)
keep_cols <- tr_reduce_frozen_design(X, FULL_CN, file.path(OUT, "reduced_design_columns.rds"),
  as.character(unique(d[hour == 0, city_hour_id])))
X <- X[, match(keep_cols, FULL_CN)]; cn <- keep_cols
if (!identical(cn, lamL$cn)) stop("Lambda file and inference design columns differ; rerun the penalty stage")
iT <- which(grepl("s\\(city_part_id\\):TEMP", cn)); iR <- which(grepl("s\\(city_part_id\\):RH", cn))
if (length(iT) != uniqueN(d$city_part_id) || length(iR) != uniqueN(d$city_part_id))
  stop("City-part slope columns do not match the observed city-part groups")
Sd_fixed <- numeric(ncol(X)); Sd_fixed[iT] <- lamT; Sd_fixed[iR] <- lamR
y <- d[[OC$hour]]; Ntot <- d[[OC$night]]; hour_vec <- d$hour
partv <- as.integer(d$night_part); PART_LEVELS <- levels(d$night_part)
RIDGE_REL <- 1e-10; GRAD_TOL <- 1e-6; NEWTON_GAIN_TOL <- 1e-6
INNER_MAXIT <- tr_env_int("INNER_MAXIT", 25L)
if (INNER_MAXIT < 1L) stop("INNER_MAXIT must be positive")

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
      if (any(!is.finite(hdiag)) || any(hdiag <= 0)) {
        bad <- which(!is.finite(hdiag) | hdiag <= 0)
        stop("Conditional Hessian has a nonpositive or nonfinite diagonal in: ",
             paste(head(colnames(Xb)[bad], 10L), collapse = ", ")) }
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
  list(beta = beta, nll = cur$nll, ginf = max(abs(cur$grad)), newton_gain = gain,
       it = it, ls_stop = ls_stop, g0 = g0, rel_g = rel_g, last_rel = last_rel) }
classify_solver <- function(r) { if (!is.finite(r$nll) || !is.finite(r$ginf) || !is.finite(r$newton_gain)) return("failed_nonfinite")
  if (r$ginf <= GRAD_TOL || (r$newton_gain >= 0 && r$newton_gain <= gain_tol_at(r$nll))) return("converged_gradient")
  if (r$ls_stop) return("failed_line_search")
  "failed_maxit" }
is_converged <- function(s) s == "converged_gradient"

wshare <- function(eta, g, p, n_part = 4L) { dt <- data.table(g = g, e = eta, p = p)
  dt[, mx := max(e), by = g]; dt[, lse := mx[1L] + log(sum(exp(e - mx[1L]))), by = g]; dt[, pp := exp(e - lse)]
  agg <- dt[, .(s = sum(pp)), by = .(g, p)]; ng <- uniqueN(g)
  M <- matrix(0, ng, n_part); M[cbind(agg$g, agg$p)] <- agg$s; M }
hprof <- function(eta, g, hb) { dt <- data.table(g = g, e = eta, h = hb); dt[, mx := max(e), by = g]
  dt[, lse := mx[1L] + log(sum(exp(e - mx[1L]))), by = g]; dt[, pp := exp(e - lse)]; dt[, .(v = mean(pp)), by = h][order(h)]$v }

block2 <- function(nm) { ifelse(grepl("^TEMP_ns[0-9]+_hour_", nm), "shared_T", ifelse(grepl("^RH_ns[0-9]+_hour_", nm), "shared_R",
  ifelse(grepl("s\\(city_part_id\\):TEMP", nm), "dev_T", ifelse(grepl("s\\(city_part_id\\):RH", nm), "dev_R", "ctrl")))) }
B2 <- block2(cn)
make_nd <- function(ex, val) { nd <- copy(d); bb <- if (ex == "TEMP") bT else bR; nb <- predict(bb, val)
  nd[[ex]] <- val
  for (j in seq_len(SPLINE_DF)) { nd[[paste0(ex, "_ns", j)]] <- nb[, j]
    for (i in seq_along(CH)) nd[[paste0(ex, "_ns", j, "_hour_", i)]] <- nb[, j] * (d$hour == as.numeric(CH[i])) }; nd }
frozen_X <- function(nd, fam) { Xn <- as(mgcv::predict.gam(fit0, newdata = as.data.frame(nd), type = "lpmatrix"), "dgCMatrix")
  Xc <- X; fc <- if (fam == "TEMP") which(B2 %in% c("shared_T", "dev_T")) else which(B2 %in% c("shared_R", "dev_R"))
  Xc[, fc] <- Xn[, match(cn[fc], colnames(Xn))]; Xc }   # invariant: Xc[-focal] == X[-focal] by construction
qT <- d[, .(lo = as.numeric(quantile(TEMP, .025)), hi = as.numeric(quantile(TEMP, .975))), by = city]
qR <- d[, .(lo = as.numeric(quantile(RH, .025)), hi = as.numeric(quantile(RH, .975))), by = city]
CF <- list(Tp_lo = frozen_X(make_nd("TEMP", rep(quantile(d$TEMP, .025), nrow(d))), "TEMP"), Tp_hi = frozen_X(make_nd("TEMP", rep(quantile(d$TEMP, .975), nrow(d))), "TEMP"),
  Rp_lo = frozen_X(make_nd("RH", rep(quantile(d$RH, .025), nrow(d))), "RH"), Rp_hi = frozen_X(make_nd("RH", rep(quantile(d$RH, .975), nrow(d))), "RH"),
  Tc_lo = frozen_X(make_nd("TEMP", qT$lo[match(d$city, qT$city)]), "TEMP"), Tc_hi = frozen_X(make_nd("TEMP", qT$hi[match(d$city, qT$city)]), "TEMP"),
  Rc_lo = frozen_X(make_nd("RH", qR$lo[match(d$city, qR$city)]), "RH"), Rc_hi = frozen_X(make_nd("RH", qR$hi[match(d$city, qR$city)]), "RH"))
estimand <- function(beta, idx, gb, cityb) { cidx <- split(seq_along(cityb), cityb)
  pv <- partv[idx]
  pl <- function(a, z) 100 * (colMeans(wshare(as.numeric(CF[[z]][idx, ] %*% beta), gb, pv)) -
                              colMeans(wshare(as.numeric(CF[[a]][idx, ] %*% beta), gb, pv)))
  pc <- function(a, z) { sl <- wshare(as.numeric(CF[[a]][idx, ] %*% beta), gb, pv)
    sh <- wshare(as.numeric(CF[[z]][idx, ] %*% beta), gb, pv)
    out <- t(sapply(cidx, function(j) 100 * (colMeans(sh[j, , drop = FALSE]) - colMeans(sl[j, , drop = FALSE]))))
    dimnames(out) <- list(names(cidx), PART_LEVELS); out }
  list(pT = pl("Tp_lo", "Tp_hi"), pR = pl("Rp_lo", "Rp_hi"), cT = pc("Tc_lo", "Tc_hi"), cR = pc("Rc_lo", "Rc_hi")) }

full <- solve_cond(X, y, Ntot, d$grp, hour_vec, Sd_fixed); cityb0 <- as.character(unique(d[, .(grp, city)])[order(grp)]$city)
if (!is_converged(classify_solver(full))) stop("Full-sample conditional fit did not converge: ", classify_solver(full), "; grad_inf=", signif(full$ginf, 5))
pe <- estimand(full$beta, seq_len(nrow(d)), d$grp, cityb0)
cat("full-sample: pooled_T =", signif(pe$pT, 5), " pooled_R =", signif(pe$pR, 5), " (status:", classify_solver(full), " ginf:", signif(full$ginf, 2), " newton_gain:", signif(full$newton_gain, 2), ")\n")

set.seed(SEED0); seeds <- sample.int(.Machine$integer.max, N_BOOT)
saveRDS(seeds, file.path(OUT, "bootstrap_seeds.rds"))
fwrite(data.table(seed0 = SEED0, block_length_days = BLOCK_LENGTH, n_boot = N_BOOT, ridge_rel = RIDGE_REL,
  grad_tol = GRAD_TOL, newton_gain_tol = NEWTON_GAIN_TOL, outcome = OC$outcome,
  city_slope_groups = "h00_02_h03_05_h06_08_h09_11", n_city_slope_groups = uniqueN(d$city_part_id),
  sampling_unit = "circular_calendar_date_block_shared_across_36_cities", reduced_design = TRUE),
  file.path(OUT, OC$out("bootstrap_specification", "csv")))
fwrite(rbindlist(lapply(seq_len(N_BOOT), function(i) tr_calendar_block_map(d$DATE, BLOCK_LENGTH, seeds[i], i))),
  file.path(OUT, "bootstrap_calendar_date_plan.csv"))
d2 <- d[, .(DATE, city, hour)]
active_bootstrap_columns <- function(Xb, Sd) {
  if (ncol(Xb) != length(Sd)) stop("Bootstrap design and penalty lengths differ")
  supported <- Matrix::colSums(abs(Xb)) > 0
  unsupported <- unname(which(!supported & Sd == 0))
  is_control <- grepl("^(calendar_hour_|sunset_dev_hour_)", colnames(Xb)[unsupported])
  if (any(!is_control))
    stop("Bootstrap has unsupported non-control columns: ",
         paste(colnames(Xb)[unsupported[!is_control]], collapse = ", "))
  dropped <- unsupported
  list(active = setdiff(seq_len(ncol(Xb)), dropped), dropped = dropped)
}
rep_one <- function(i) { setDTthreads(1); t0 <- Sys.time()
  ans <- tryCatch({
    map <- tr_calendar_block_map(d$DATE, BLOCK_LENGTH, seeds[i]); map[, replicate := NULL]
    mrg <- merge(map[, .(boot_day, source_DATE)], d2[, .(row = .I, DATE, city, hour)], by.x = "source_DATE", by.y = "DATE", allow.cartesian = TRUE, sort = FALSE)
    setorder(mrg, city, boot_day, hour); mrg[, gb := .GRP, by = .(city, boot_day)]
    idx <- mrg$row; gb <- mrg$gb; hb <- mrg$hour; cityb <- as.character(mrg[, city[1L], by = gb][order(gb)]$V1)
    Xb <- X[idx, , drop = FALSE]
    cols <- active_bootstrap_columns(Xb, Sd_fixed)
    r <- solve_cond(Xb[, cols$active, drop = FALSE], y[idx], Ntot[idx], gb, hb,
                    Sd_fixed[cols$active])
    beta_full <- numeric(ncol(X)); beta_full[cols$active] <- r$beta; r$beta <- beta_full
    e <- estimand(r$beta, idx, gb, cityb)
    bs <- r$beta; bs[c(iT, iR)] <- 0   # shared-only prediction (RE block zeroed)
    es <- estimand(bs, idx, gb, cityb)
    list(ok = TRUE, i = i, status = classify_solver(r), ginf = r$ginf, newton_gain = r$newton_gain, it = r$it,
      dropped_fixed = length(cols$dropped), dropped_fixed_names = paste(cn[cols$dropped], collapse = ";"),
      secs = as.numeric(difftime(Sys.time(), t0, units = "secs")), pT = e$pT, pR = e$pR,
      pT_s = es$pT, pR_s = es$pR, beta = as.numeric(r$beta),
      cT = e$cT, cR = e$cR, prof_T = 100 * hprof(as.numeric(CF$Tp_hi[idx, ] %*% r$beta), gb, hb) - 100 * hprof(as.numeric(CF$Tp_lo[idx, ] %*% r$beta), gb, hb),
      prof_R = 100 * hprof(as.numeric(CF$Rp_hi[idx, ] %*% r$beta), gb, hb) - 100 * hprof(as.numeric(CF$Rp_lo[idx, ] %*% r$beta), gb, hb), error = NA_character_)
  }, error = function(err) list(ok = FALSE, i = i, error = conditionMessage(err), secs = as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  ans }
cat("running", N_BOOT, "replicates on", N_CORES, "cores ...\n")
res <- vector("list", N_BOOT); chunk <- 100L
for (s in seq(1L, N_BOOT, by = chunk)) { ix <- s:min(s + chunk - 1L, N_BOOT)
  res[ix] <- tr_parallel_lapply(ix, rep_one, N_CORES)
  cat(sprintf("  %d/%d replicates done (%s)\n", max(ix), N_BOOT, format(Sys.time(), "%H:%M:%S"))) }

status <- rbindlist(lapply(res, function(r) data.table(replicate = r$i, ok = isTRUE(r$ok),
  termination_status = if (isTRUE(r$ok)) r$status else NA_character_,
  converged = if (isTRUE(r$ok)) is_converged(r$status) else NA,
  ginf = if (isTRUE(r$ok)) r$ginf else NA_real_, newton_gain = if (isTRUE(r$ok)) r$newton_gain else NA_real_,
  iters = if (isTRUE(r$ok)) r$it else NA_integer_,
  dropped_fixed = if (isTRUE(r$ok)) r$dropped_fixed else NA_integer_,
  dropped_fixed_names = if (isTRUE(r$ok)) r$dropped_fixed_names else NA_character_,
  secs = r$secs, error = r$error %||% NA_character_)))
STATUS_FILE <- file.path(OUT, OC$out("bootstrap_status", "csv"))
fwrite(status, STATUS_FILE)
good <- Filter(function(r) isTRUE(r$ok) && is_converged(r$status), res)
cat("replicates ok:", length(good), "/", N_BOOT, "; median secs:", signif(median(status$secs, na.rm = TRUE), 3), "\n")
if (length(good) < 2L) stop("Fewer than two converged bootstrap replicates; inspect ", basename(STATUS_FILE))

pooled <- rbindlist(lapply(c("TEMP", "RH"), function(ex) {
  M <- t(vapply(good, function(r) if (ex == "TEMP") r$pT else r$pR, numeric(length(PART_LEVELS))))
  pt <- if (ex == "TEMP") pe$pT else pe$pR
  data.table(exposure = ex, part = PART_LEVELS, estimate = pt, boot_SE = apply(M, 2, sd),
    q025 = apply(M, 2, quantile, .025), q975 = apply(M, 2, quantile, .975),
    CI_lo = pt - 1.96 * apply(M, 2, sd), CI_hi = pt + 1.96 * apply(M, 2, sd), n_success = nrow(M)) }))
cities <- sort(unique(as.character(d$city)))
city_tab <- rbindlist(lapply(c("TEMP", "RH"), function(ex) { pt <- if (ex == "TEMP") pe$cT else pe$cR
  rbindlist(lapply(PART_LEVELS, function(pn) { M <- t(vapply(good, function(r) (if (ex == "TEMP") r$cT else r$cR)[, pn], numeric(length(cities))))
    est <- pt[, pn]
    data.table(exposure = ex, city = cities, part = pn, estimate = est, boot_SE = apply(M, 2, sd),
      q025 = apply(M, 2, quantile, .025), q975 = apply(M, 2, quantile, .975),
      CI_lo = est - 1.96 * apply(M, 2, sd), CI_hi = est + 1.96 * apply(M, 2, sd), n_success = nrow(M)) })) }))
fwrite(pooled, file.path(OUT, OC$out("bootstrap_pooled", "csv"))); fwrite(city_tab, file.path(OUT, OC$out("bootstrap_36city", "csv")))
saveRDS(res, file.path(OUT, OC$out(sprintf("bootstrap_replicates_B%d", N_BOOT), "rds")))
cat("\n=== pooled per-window share contrast (pp), fixed-lambda block bootstrap ===\n"); print(pooled, digits = 5)
cat("\n=== 36-city per-window change (pp): estimate & 95% CI ===\n"); print(city_tab, digits = 4)
cat("\n--- failure / convergence diagnostics ---\n")
cat("non-ok (error):", sum(!status$ok), "\n")
print(status[ok == TRUE, .N, by = termination_status])
cat("usable (converged incl. numerical-floor):", sum(status$converged, na.rm = TRUE), "/", N_BOOT, "\n")
message("Saved ", OC$out("bootstrap_pooled", "csv"), " / ", OC$out("bootstrap_36city", "csv"),
        " / ", OC$out("bootstrap_status", "csv"), " / ", OC$out(sprintf("bootstrap_replicates_B%d", N_BOOT), "rds"))
