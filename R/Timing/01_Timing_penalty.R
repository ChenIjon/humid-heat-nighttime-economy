#!/usr/bin/env Rscript
# Paper: Heat and Humidity Reshape Urban Economic Activity at Night
# Contact: chenchen1107@connect.hku.hk
# Select penalties from the local panel; output: laplace_lambda*.rds.

.code_dir <- Sys.getenv("TIMING_REVISION_DIR", unset = "")
if (!nzchar(.code_dir)) {
  .script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  .code_dir <- if (length(.script_arg)) dirname(normalizePath(sub("^--file=", "", .script_arg[1]), winslash = "/")) else getwd()
}
source(file.path(.code_dir, "timing_helpers.R"))

OUT <- TIMING_OUTPUT_DIR; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SPLINE_DF <- 3L; N_CORES <- tr_env_int("N_CORES", 20L); SMOKE_ONLY <- tr_env_flag("SMOKE_ONLY", FALSE); setDTthreads(N_CORES)
OC <- tr_outcome_spec(); cat("outcome:", OC$outcome, "(hour =", OC$hour, ", night total =", OC$night, ")\n")
LAMBDA_FILE <- file.path(OUT, OC$lambda_file)
GRID_FILE <- file.path(OUT, OC$out("laplace_grid", "csv"))
SMOKE_FILE <- file.path(OUT, OC$out("penalty_smoke", "csv"))

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
nl <- c(TS, RS); ct <- c(SC, CCal); d[, alpha := log(get(OC$night))]
form <- as.formula(paste(OC$hour, "~", paste(c(nl, ct, "factor(city_hour_id)",
  "s(city_part_id, by = TEMP, bs = 're')", "s(city_part_id, by = RH, bs = 're')"), collapse = " + "), "+ offset(alpha)"))
fit0 <- bam(form, data = d, family = poisson(link = "log"), method = "fREML", discrete = TRUE, nthreads = N_CORES, control = gam.control(maxit = 200))
X <- as(mgcv::predict.gam(fit0, newdata = as.data.frame(d), type = "lpmatrix"), "dgCMatrix"); cn <- colnames(X)
keep_cols <- tr_reduce_frozen_design(X, cn, file.path(OUT, "reduced_design_columns.rds"),
  as.character(unique(d[hour == 0, city_hour_id])))
X <- X[, match(keep_cols, cn)]; cn <- keep_cols
iT <- which(grepl("s\\(city_part_id\\):TEMP", cn)); iR <- which(grepl("s\\(city_part_id\\):RH", cn)); iU <- c(iT, iR); iF <- setdiff(seq_along(cn), iU)
if (length(iT) != uniqueN(d$city_part_id) || length(iR) != uniqueN(d$city_part_id))
  stop("City-part slope columns do not match the observed city-part groups")
y <- d[[OC$hour]]; Ntot <- d[[OC$night]]; grp <- d$grp; hidx <- d$hour; n <- nrow(d); gstart <- d[, .(start = .I[1L]), by = grp][order(grp)]$start
RIDGE_REL <- 1e-10; GRAD_TOL <- 1e-6; NEWTON_GAIN_TOL <- 1e-6
INNER_MAXIT <- tr_env_int("INNER_MAXIT", 25L)
if (INNER_MAXIT < 1L) stop("INNER_MAXIT must be positive")
cat("design p =", ncol(X), " q_T =", length(iT), " q_R =", length(iR), "\n")

gain_tol_at <- function(nll) max(NEWTON_GAIN_TOL, .Machine$double.eps * max(1, abs(nll)))
gain_pass <- function(nll, ginf, newton_gain) {
  tol <- gain_tol_at(nll)
  if (!is.finite(nll) || !is.finite(ginf) || !is.finite(newton_gain)) return("nonfinite")
  if (ginf <= GRAD_TOL) return("gradient")
  if (newton_gain >= 0 && newton_gain <= tol) return("newton_gain")
  "none"
}
solve_cond <- function(Sd) { beta <- numeric(ncol(X))
  pnf <- function(eta) { dt <- data.table(grp = grp, eta = eta); dt[, mx := max(eta), by = grp]; dt[, lse := mx[1L] + log(sum(exp(eta - mx[1L]))), by = grp]; dt[, p := exp(eta - lse)]; dt[, N := Ntot]; dt }
  og <- function(b, wH = TRUE) { eta <- as.numeric(X %*% b); pn <- pnf(eta); lse <- pn$lse; p <- pn$p
    nll <- sum(pn[, .(v = N[1L] * lse[1L]), by = grp]$v) - sum(y * eta) + 0.5 * sum(Sd * b^2)
    grad <- as.numeric(crossprod(X, Ntot * p - y)) + Sd * b
    H <- NULL; ridge_min_abs <- NA_real_; ridge_max_abs <- NA_real_
    if (wH) { ii <- rep(seq_len(n), times = 12L); jj <- integer(n * 12L); xx <- numeric(n * 12L); pos <- 0L
      for (bb in 0:11) { colb <- gstart[grp] + bb; val <- Ntot * p * (as.numeric(hidx == bb) - p[colb]); jj[(pos + 1L):(pos + n)] <- colb; xx[(pos + 1L):(pos + n)] <- val; pos <- pos + n }
      H <- crossprod(X, sparseMatrix(i = ii, j = jj, x = xx, dims = c(n, n)) %*% X) + Diagonal(x = Sd)
      H <- as(H, "symmetricMatrix"); hdiag <- Matrix::diag(H)
      if (any(!is.finite(hdiag)) || any(hdiag <= 0))
        stop("Conditional Hessian has a nonpositive or nonfinite diagonal")
      ridge <- RIDGE_REL * hdiag
      ridge_min_abs <- min(ridge); ridge_max_abs <- max(ridge)
      diag(H) <- hdiag + ridge }
    list(nll = nll, grad = grad, H = H,
         ridge_min_abs = ridge_min_abs, ridge_max_abs = ridge_max_abs) }
  cur <- og(beta)
  ls_stop <- FALSE; it <- 0L; g0 <- max(abs(cur$grad)); accepted_steps <- 0L
  last_step_scale <- NA_real_; last_rel <- NA_real_
  for (it in seq_len(INNER_MAXIT)) { step <- as.numeric(solve(cur$H, -cur$grad))
    gain <- -0.5 * sum(cur$grad * step)
    if (gain_pass(cur$nll, max(abs(cur$grad)), gain) != "none") break
    if (!is.finite(gain) || gain < 0) { ls_stop <- TRUE; break }
    tt <- 1; ok <- FALSE
    while (tt > 1e-6) { cd <- beta + tt * step; if (og(cd, FALSE)$nll < cur$nll - 1e-12) { ok <- TRUE; break }; tt <- tt / 2 }
    if (!ok) { ls_stop <- TRUE; break }
    prev_nll <- cur$nll; beta <- beta + tt * step; cur <- og(beta)
    accepted_steps <- accepted_steps + 1L; last_step_scale <- tt
    last_rel <- (prev_nll - cur$nll) / (1 + abs(prev_nll)) }
  step <- as.numeric(solve(cur$H, -cur$grad)); gain <- -0.5 * sum(cur$grad * step)
  list(beta = beta, nll = cur$nll, H = cur$H, grad = cur$grad,
       ridge_min_abs = cur$ridge_min_abs, ridge_max_abs = cur$ridge_max_abs,
       ginf = max(abs(cur$grad)),
       newton_gain = gain, it = it, ls_stop = ls_stop, accepted_steps = accepted_steps,
       gain_tol_eff = gain_tol_at(cur$nll), gain_pass = gain_pass(cur$nll, max(abs(cur$grad)), gain),
       g0 = g0, rel_g = max(abs(cur$grad)) / max(g0, .Machine$double.xmin),
       last_step_scale = last_step_scale, last_rel = last_rel) }
inner_ok <- function(r) is.finite(r$nll) && is.finite(r$ginf) && is.finite(r$newton_gain) &&
  (r$ginf <= GRAD_TOL || (r$newton_gain >= 0 && r$newton_gain <= gain_tol_at(r$nll)))
ld <- function(M) { Mc <- as.matrix(M); Mc <- (Mc + t(Mc)) / 2; j <- 0
  repeat { ch <- tryCatch(chol(Mc + diag(j, nrow(Mc))), error = function(e) NULL); if (!is.null(ch)) return(2 * sum(log(diag(ch)))); j <- if (j == 0) 1e-8 * mean(diag(Mc)) else j * 10 } }

Mfun <- function(th) { lamT <- exp(th[1]); lamR <- exp(th[2]); Sd <- numeric(ncol(X)); Sd[iT] <- lamT; Sd[iR] <- lamR
  r <- solve_cond(Sd)
  diagnostics <- list(grad_inf = r$ginf, newton_gain = r$newton_gain,
                      penalized_nll = r$nll, gain_tol_eff = r$gain_tol_eff,
                      gain_pass = r$gain_pass,
                      grad_inf_fixed = max(abs(r$grad[iF])),
                      grad_inf_RE = max(abs(r$grad[iU])),
                      ridge_min_abs = r$ridge_min_abs, ridge_max_abs = r$ridge_max_abs,
                      it = r$it, ls_stop = r$ls_stop, accepted_steps = r$accepted_steps,
                      initial_grad_inf = r$g0, relative_grad = r$rel_g,
                      last_step_scale = r$last_step_scale, last_relative_decrease = r$last_rel)
  if (!inner_ok(r)) return(as.data.table(c(list(M = Inf, status = "failed_inner",
    failure_reason = if (r$ls_stop) "line_search" else "max_iterations"), diagnostics)))
  Hff <- as.matrix(r$H[iF, iF]); Hfu <- as.matrix(r$H[iF, iU]); Huu <- as.matrix(r$H[iU, iU])
  Hprofile <- Huu - crossprod(Hfu, solve(Hff, Hfu))
  as.data.table(c(list(M = as.numeric(r$nll + .5 * ld(Hprofile) -
    .5 * (length(iT) * log(lamT) + length(iR) * log(lamR))),
    status = "ok", failure_reason = NA_character_), diagnostics)) }

sp_names <- names(fit0$sp)
sp_T <- which(grepl("city_part_id.*TEMP", sp_names)); sp_R <- which(grepl("city_part_id.*RH", sp_names))
if (length(sp_T) != 1L || length(sp_R) != 1L) stop("Cannot identify the two bam starting penalties")
start_lam <- as.numeric(fit0$sp[c(sp_T, sp_R)])
if (any(!is.finite(start_lam)) || any(start_lam <= 0)) stop("bam starting penalties must be positive and finite")
base <- log(start_lam)
cat("bam starting penalties: lam_T =", signif(start_lam[1], 7), " lam_R =", signif(start_lam[2], 7), "\n")
if (SMOKE_ONLY) {
  smoke <- tryCatch(Mfun(base), error = function(e)
    data.table(M = NA_real_, grad_inf = NA_real_, newton_gain = NA_real_,
               penalized_nll = NA_real_, gain_tol_eff = NA_real_, gain_pass = NA_character_,
               grad_inf_fixed = NA_real_, grad_inf_RE = NA_real_,
               ridge_min_abs = NA_real_, ridge_max_abs = NA_real_,
               it = NA_integer_, ls_stop = NA, accepted_steps = NA_integer_,
               initial_grad_inf = NA_real_, relative_grad = NA_real_,
               last_step_scale = NA_real_, last_relative_decrease = NA_real_,
               status = "error", failure_reason = conditionMessage(e)))
  smoke[, `:=`(n_rows = nrow(d), n_nights = uniqueN(d$night_id),
               n_cities = uniqueN(d$city), n_city_part = uniqueN(d$city_part_id),
               design_p = ncol(X), q_TEMP = length(iT), q_RH = length(iR),
               start_lam_TEMP = start_lam[1L], start_lam_RH = start_lam[2L],
               grad_tol = GRAD_TOL, newton_gain_tol = NEWTON_GAIN_TOL,
               inner_maxit = INNER_MAXIT, outcome = OC$outcome)]
  smoke_path <- SMOKE_FILE
  fwrite(smoke, smoke_path)
  print(smoke)
  cat("smoke diagnostic:", smoke_path, "\n")
  if (!identical(smoke$status, "ok") || !is.finite(smoke$M))
    stop("Penalty smoke failed; inspect ", basename(smoke_path), " and the log")
  cat("Penalty smoke passed; no lambda was saved.\n")
  quit(save = "no", status = 0L)
}
run_grid <- function(ctr, step, rad) { g <- as.data.table(expand.grid(a = seq(-rad, rad, by = step), b = seq(-rad, rad, by = step)))
  g <- cbind(g, rbindlist(lapply(seq_len(nrow(g)), function(k) Mfun(ctr + c(g$a[k], g$b[k])))))
  if (!any(is.finite(g$M))) stop("No converged inner fits at grid step ", step, "; inspect grad_inf and newton_gain")
  cat("grid step", step, ":", sum(is.finite(g$M)), "/", nrow(g), "valid; best M =", signif(min(g$M), 10), " at offset", unlist(g[which.min(M), .(a, b)]), "\n"); g }
g1 <- run_grid(base, 1, 2)                     # 5x5
c1 <- base + c(g1[which.min(M), a], g1[which.min(M), b])
g2 <- run_grid(c1, 0.25, 0.5)                  # 5x5 refine
c2 <- c1 + c(g2[which.min(M), a], g2[which.min(M), b])
lam <- exp(c2); cat("\n>>> Laplace optimum: lam_T =", signif(lam[1], 7), "(x", signif(lam[1]/start_lam[1], 3), ") lam_R =", signif(lam[2], 7), "(x", signif(lam[2]/start_lam[2], 3), ") vs bam start\n")
Sd <- numeric(ncol(X)); Sd[iT] <- lam[1]; Sd[iR] <- lam[2]; rr <- solve_cond(Sd)
if (!inner_ok(rr)) stop("Selected lambda did not pass the inner-fit check: gradient=", rr$ginf, ", newton_gain=", rr$newton_gain)
be <- rr$beta; Hf <- rr$H
fwrite(rbind(g1[, .(grid = "coarse", a, b, M, grad_inf, newton_gain, penalized_nll, gain_tol_eff, gain_pass,
                    grad_inf_fixed, grad_inf_RE, ridge_min_abs, ridge_max_abs, it, ls_stop,
                    accepted_steps, initial_grad_inf, relative_grad, last_step_scale,
                    last_relative_decrease, status, failure_reason)],
            g2[, .(grid = "refine", a, b, M, grad_inf, newton_gain, penalized_nll, gain_tol_eff, gain_pass,
                    grad_inf_fixed, grad_inf_RE, ridge_min_abs, ridge_max_abs, it, ls_stop,
                    accepted_steps, initial_grad_inf, relative_grad, last_step_scale,
                    last_relative_decrease, status, failure_reason)]), GRID_FILE)
Hinv_diag <- Matrix::diag(solve(Hf))
edf_T <- length(iT) - sum(Hinv_diag[iT] * Sd[iT]); edf_R <- length(iR) - sum(Hinv_diag[iR] * Sd[iR])
cat("RE_TEMP: ref.df=", length(iT), " edf=", signif(edf_T, 4), " beta_rms=", signif(sqrt(mean(be[iT]^2)), 4),
    "\nRE_RH  : ref.df=", length(iR), " edf=", signif(edf_R, 4), " beta_rms=", signif(sqrt(mean(be[iR]^2)), 4), "\n")
saveRDS(list(lam = lam, beta = be, cn = cn, objective_version = "conditional_city_quarter_laplace_v5",
             outcome = OC$outcome, hour = OC$hour, night = OC$night,
             city_slope_groups = c("h00_02", "h03_05", "h06_08", "h09_11"), bam_start_lam = start_lam,
             inner_grad_inf = rr$ginf, inner_newton_gain = rr$newton_gain,
             edf_TEMP_dev = edf_T, edf_RH_dev = edf_R,
             grad_tol = GRAD_TOL, newton_gain_tol = NEWTON_GAIN_TOL),
        LAMBDA_FILE)
message("Saved ", basename(LAMBDA_FILE), " / ", basename(GRID_FILE))
