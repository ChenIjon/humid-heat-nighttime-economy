################################################################################
# Frozen NEI main-model engine (source-only).
#
# T1 structure is frozen from compare_time_control_models4.R:
# TEMP, RH, TEMP-RH, city-specific TEMP/RH deviations, precipitation, PM2.5,
# wind, month FE, city-specific date smooth (k=10), calendar type, city RE,
# and model-specific AR(1). Inference targets the headline exposure-weighted
# marginal effect (EWME) of compare_time_control_models3.R / models4.R
# (central-95% exposure grid, trapezoidal density weights, link-scale design
# difference over the exposure smooths). Circular calendar-date moving-block
# bootstrap refits estimate the standard error. For multiplicative outcomes,
# the reported 95% CI is centred on the original-model Link_EWME and transformed
# back to percent change. For logit-share outcomes, it is centred on the
# original-model percentage-point EWME on its reported scale.
# Families: nb(log) for counts, Gamma(log) for positive spending, gaussian()
# for continuous responses, and betar(logit) for 0-1 share responses.
#
# A wrapper prepares a panel with city, DATE, TEMP, RH, WS, PRECIP, PM25,
# calendar_type and its response variables, then calls main_model_run().
################################################################################

suppressPackageStartupMessages({ library(mgcv); library(dplyr) })

MAIN_MODEL_VERSION <- "frozen_T1_AR1_date_block_bootstrap_normal_CI_EWME_v3_2026-09-19"
MAIN_CITIES <- c("Beijing","Changsha","Chengdu","Fuzhou","Guangzhou","Guiyang",
  "Hangzhou","Ningbo","Shanghai","Chongqing","Shenyang","Wuhan","Hefei",
  "Huhehaote","Nanjing","Nanning","Qingdao","Shenzhen","Wulumuqi","Xining",
  "Yinchuan","Dalian","Xian","Kunming","Nanchang","Lanzhou","Zhengzhou",
  "Xiamen","Haikou","Haerbin","Shijiazhuang","Jinan","Lasa","Taiyuan",
  "Changchun","Tianjin")
MAIN_CALENDAR_LEVELS <- c("Monday","Tuesday","Wednesday","Thursday","Friday",
  "Saturday","Sunday","Statutory holiday","Make-up workday")

main_model_defaults <- function() list(
  # seed = 20260912 matches block_bootstrap_final_T1_EWME.R (SEED) so the
  # calendar-date block bootstrap here reproduces the final main-model results.
  rho_abs_bound = .95, block_length = 7L, n_boot = 1000L,
  seed = 20260912L, fix_sp = TRUE, force_refit = FALSE,
  export_bootstrap_curves = FALSE, curve_grid_n = 100L, n_cores = 1L)

main_model_require_outcomes <- function(outcomes) {
  required <- c("Outcome", "Response", "Family", "Response_scale")
  if (!is.data.frame(outcomes) || length(setdiff(required, names(outcomes)))) {
    stop("outcomes must contain: ", paste(required, collapse = ", "))
  }
  if (!all(outcomes$Family %in% c("nb", "gamma", "gaussian", "beta"))) stop("Invalid Family")
  if (!all(outcomes$Response_scale %in% c("multiplicative", "logit_share"))) stop("Invalid Response_scale")
}

main_model_prepare_panel <- function(data, responses) {
  required <- c("city","DATE","TEMP","RH","WS","PRECIP","PM25","calendar_type",responses)
  missing <- setdiff(required, names(data))
  if (length(missing)) stop("Panel lacks: ", paste(missing, collapse = ", "))
  d <- data
  d$DATE <- as.Date(d$DATE)
  if (anyNA(d$DATE) || anyDuplicated(d[, c("city","DATE")])) stop("DATE missing or city-date duplicates")
  d <- d[d$city %in% MAIN_CITIES, , drop = FALSE]
  d$city_factor <- factor(d$city, levels = MAIN_CITIES)
  d$date_numeric <- as.numeric(d$DATE - min(d$DATE))
  d$month <- factor(format(d$DATE, "%Y-%m"))
  d$calendar_type <- factor(d$calendar_type, levels = MAIN_CALENDAR_LEVELS)
  needed <- c(responses,"TEMP","RH","WS","PRECIP","PM25","calendar_type","city_factor","date_numeric","month")
  keep <- complete.cases(d[, needed])
  for (v in c(responses,"TEMP","RH","WS","PRECIP","PM25")) keep <- keep & is.finite(d[[v]])
  keep <- keep & d$RH >= 0 & d$RH <= 100 & d$PRECIP >= 0
  d <- d[keep, , drop = FALSE]
  d <- d[order(d$city_factor, d$DATE), , drop = FALSE]; rownames(d) <- NULL
  for (v in c("calendar_type","city_factor","month")) d[[v]] <- droplevels(d[[v]])
  d$AR_start <- c(TRUE, d$city_factor[-1L] != d$city_factor[-nrow(d)] |
    as.integer(d$DATE[-1L] - d$DATE[-nrow(d)]) != 1L)
  d
}

main_model_formula <- function(response) {
  city_temp <- sprintf('ti(TEMP,city_factor,bs=c("cr","re"),k=c(3,%d))', length(MAIN_CITIES))
  city_rh <- sprintf('ti(RH,city_factor,bs=c("cr","re"),k=c(3,%d))', length(MAIN_CITIES))
  terms <- c('ti(TEMP,bs="cr",k=5)','ti(RH,bs="cr",k=5)',
    'ti(TEMP,RH,bs=c("cr","cr"),k=c(5,5))', city_temp, city_rh,
    's(PRECIP,bs="cr",k=3)','s(PM25,bs="cr",k=3)','s(WS,bs="cr",k=5)',
    'factor(month)','s(date_numeric,by=city_factor,bs="cr",k=10)',
    'factor(calendar_type)','s(city_factor,bs="re")')
  as.formula(paste(response, "~", paste(terms, collapse = " + ")))
}

main_model_clamp_nb_deviance <- function(family) {
  # bam whitens AR(1) residuals with sqrt(dev.resids); clamp tiny negative NB
  # deviances to 0 so sqrt(-x) cannot abort a bootstrap refit.
  orig <- family$dev.resids
  family$dev.resids <- function(y, mu, wt, ...) {
    d <- orig(y, mu, wt, ...); d[!is.na(d) & d < 0] <- 0; d
  }
  family
}

main_model_family <- function(spec, theta = NA_real_) {
  if (spec$Family == "gaussian") return(gaussian())
  if (spec$Family == "beta") return(betar(link = "logit"))
  if (spec$Family == "gamma") return(Gamma(link = "log"))
  if (is.finite(theta) && theta > 0) {
    main_model_clamp_nb_deviance(nb(theta = theta, link = "log"))
  } else {
    nb(link = "log")
  }
}

main_model_theta <- function(fit, spec) {
  if (spec$Family != "nb") return(NA_real_)
  tryCatch(as.numeric(fit$family$getTheta(TRUE)), error = function(e) NA_real_)
}

main_model_rho <- function(fit, d, bound) {
  r <- as.numeric(residuals(fit, type = "working"))
  city_rhos <- unlist(lapply(split(seq_along(r), d$city), function(i) {
    left <- i[-length(i)]; right <- i[-1L]
    ok <- as.integer(d$DATE[right] - d$DATE[left]) == 1L
    if (sum(ok) < 4L) return(NA_real_)
    cor(r[left][ok], r[right][ok])
  }))
  city_rhos <- city_rhos[is.finite(city_rhos)]
  if (!length(city_rhos)) stop("No finite city working-residual lag-1 correlation")
  max(-bound, min(bound, median(city_rhos)))
}

# ---------------------------------------------------------------------------
# Headline EWME (exposure-weighted marginal effect), aligned with the EWME
# calculation in compare_time_control_models3.R / compare_time_control_models4.R
# and block_bootstrap_final_T1_EWME.R:
#   * exposure grid over the central 95% range (301 points);
#   * trapezoidal density weights;
#   * link-scale design difference from a p2.5 reference, restricted to the
#     selected exposure smooths.
# Effects: 100*(exp(eta)-1) percent change for log-link families; density-
# weighted change in predicted probability (percentage points) for logit-share.
# ---------------------------------------------------------------------------
main_model_norm_label <- function(x) gsub("\\s+", "", x)

main_model_term_columns <- function(fit, labels) {
  labs <- vapply(fit$smooth, function(x) main_model_norm_label(x$label), character(1))
  idx <- match(main_model_norm_label(labels), labs)
  if (anyNA(idx)) stop("Exposure smooth not found: ", paste(labels[is.na(idx)], collapse = ", "))
  unlist(lapply(fit$smooth[idx], function(x) seq.int(x$first.para, x$last.para)))
}

main_model_design_difference <- function(fit, newdata, refdata, labels) {
  # Explicit lpmatrix prediction avoids predict.bam's discrete tensor
  # reassembly issue for synthetic exposure grids that contain re margins.
  X <- mgcv::predict.gam(fit, newdata = newdata, type = "lpmatrix")
  Xr <- mgcv::predict.gam(fit, newdata = refdata, type = "lpmatrix")
  if (nrow(Xr) == 1L) Xr <- Xr[rep(1L, nrow(X)), , drop = FALSE]
  D <- X - Xr
  D[, -main_model_term_columns(fit, labels)] <- 0
  D
}

main_model_trapz_weights <- function(x) {
  z <- diff(x); w <- numeric(length(x)); w[1] <- z[1] / 2; w[length(x)] <- z[length(z)] / 2
  if (length(x) > 2) w[2:(length(x) - 1)] <- (z[-length(z)] + z[-1]) / 2
  w
}

main_model_ewme_row <- function(fit, D, Xr, w, labels, name, response_scale) {
  w <- w / sum(w)
  keep <- main_model_term_columns(fit, labels)
  cf <- coef(fit)
  if (response_scale == "logit_share") {
    eta_ref <- sum(Xr[1L, keep] * cf[keep])
    eta <- eta_ref + as.vector(D %*% cf)
    ewme <- 100 * sum(w * (plogis(eta) - plogis(eta_ref)))
    unit <- "percentage_points"
    link <- sum(w * (eta - eta_ref))
  } else {
    db <- as.vector(crossprod(w, D))
    link <- sum(db * cf)
    ewme <- 100 * (exp(link) - 1)
    unit <- "percent_change"
  }
  data.frame(Effect = name, EWME = ewme, Effect_unit = unit, Link_EWME = link,
    stringsAsFactors = FALSE)
}

# Exposure grid / reference / density weights / joint empirical rows, computed
# once on the full analysis sample and held fixed across bootstrap replicates
# (the "fixed original exposure specification" convention of the final
# main-model block bootstrap).
main_model_ewme_spec <- function(d, grid_n = 301L) {
  make_grid <- function(exposure) {
    lo <- as.numeric(quantile(d[[exposure]], .025))
    hi <- as.numeric(quantile(d[[exposure]], .975))
    x <- seq(lo, hi, length.out = grid_n)
    den <- density(d[[exposure]], from = lo, to = hi, n = grid_n)
    w <- main_model_trapz_weights(x) * approx(den$x, den$y, xout = x, rule = 2)$y
    list(lower = lo, upper = hi, x = x, weights = w / sum(w))
  }
  grids <- list(TEMP = make_grid("TEMP"), RH = make_grid("RH"))
  keep <- d$TEMP >= grids$TEMP$lower & d$TEMP <= grids$TEMP$upper &
    d$RH >= grids$RH$lower & d$RH <= grids$RH$upper
  list(grid_n = grid_n, grids = grids,
    joint_TEMP = d$TEMP[keep], joint_RH = d$RH[keep],
    joint_weights = rep(1 / sum(keep), sum(keep)))
}

# With fixed_spec supplied, the grids/weights/joint rows come from the original
# sample (spec), so only the fitted coefficients vary across replicates. Without
# it, the specification is recomputed from d (the historical behaviour).
main_model_ewme <- function(fit, d, spec, fixed_spec = NULL) {
  base <- d[1, , drop = FALSE]
  grid_n <- if (is.null(fixed_spec)) 301L else fixed_spec$grid_n
  one_dim <- function(exposure, labels, name) {
    if (is.null(fixed_spec)) {
      lo <- as.numeric(quantile(d[[exposure]], .025))
      hi <- as.numeric(quantile(d[[exposure]], .975))
      x <- seq(lo, hi, length.out = grid_n)
      den <- density(d[[exposure]], from = lo, to = hi, n = grid_n)
      w <- main_model_trapz_weights(x) * approx(den$x, den$y, xout = x, rule = 2)$y
    } else {
      g <- fixed_spec$grids[[exposure]]; lo <- g$lower; x <- g$x; w <- g$weights
    }
    nd <- base[rep(1L, length(x)), , drop = FALSE]
    rd <- base
    nd[[exposure]] <- x
    rd[[exposure]] <- lo
    D <- main_model_design_difference(fit, nd, rd, labels)
    Xr <- mgcv::predict.gam(fit, newdata = rd, type = "lpmatrix")
    main_model_ewme_row(fit, D, Xr, w, labels, name, spec$Response_scale)
  }
  out <- rbind(
    one_dim("TEMP", "ti(TEMP)", "TEMP_MAIN"),
    one_dim("RH", "ti(RH)", "RH_MAIN")
  )
  if (is.null(fixed_spec)) {
    tr <- as.numeric(quantile(d$TEMP, .025)); rr <- as.numeric(quantile(d$RH, .025))
    keep <- d$TEMP >= tr & d$TEMP <= quantile(d$TEMP, .975) &
      d$RH >= rr & d$RH <= quantile(d$RH, .975)
    nd <- d[keep, , drop = FALSE]
    joint_w <- rep(1, nrow(nd))
  } else {
    tr <- fixed_spec$grids$TEMP$lower; rr <- fixed_spec$grids$RH$lower
    nd <- base[rep(1L, length(fixed_spec$joint_TEMP)), , drop = FALSE]
    nd$TEMP <- fixed_spec$joint_TEMP; nd$RH <- fixed_spec$joint_RH
    joint_w <- fixed_spec$joint_weights
  }
  rd <- base
  rd$TEMP <- tr; rd$RH <- rr
  labels <- c("ti(TEMP)", "ti(RH)", "ti(TEMP,RH)")
  D <- main_model_design_difference(fit, nd, rd, labels)
  Xr <- mgcv::predict.gam(fit, newdata = rd, type = "lpmatrix")
  out <- rbind(out, main_model_ewme_row(fit, D, Xr, joint_w, labels,
    "JOINT_TEMP_RH_EMPIRICAL", spec$Response_scale))
  out$Outcome <- spec$Outcome
  out[, c("Outcome", "Effect", "EWME", "Effect_unit", "Link_EWME")]
}

# Main-effect response curves use the same definition as the historical T1
# comparison plots: retain only ti(TEMP) or ti(RH), excluding the interaction
# and city-specific deviation terms.  The uncertainty is added below from the
# calendar-date block bootstrap refits, never from model covariance Wald bands.
main_model_curve_grid <- function(d, n = 100L) {
  make_grid <- function(exposure) {
    x <- seq(min(d[[exposure]], na.rm = TRUE), max(d[[exposure]], na.rm = TRUE), length.out = n)
    data.frame(Exposure = exposure, Exposure_value = x,
      Reference_value = as.numeric(quantile(d[[exposure]], .025, na.rm = TRUE)),
      In_central_95pct = x >= quantile(d[[exposure]], .025, na.rm = TRUE) &
        x <= quantile(d[[exposure]], .975, na.rm = TRUE), stringsAsFactors = FALSE)
  }
  rbind(make_grid("TEMP"), make_grid("RH"))
}

main_model_response_curves <- function(fit, d, spec, curve_grid) {
  base <- d[1, , drop = FALSE]
  out <- lapply(split(curve_grid, curve_grid$Exposure), function(g) {
    exposure <- g$Exposure[1]
    nd <- base[rep(1L, nrow(g)), , drop = FALSE]
    rd <- base
    nd[[exposure]] <- g$Exposure_value
    rd[[exposure]] <- g$Reference_value[1]
    D <- main_model_design_difference(fit, nd, rd, paste0("ti(", exposure, ")"))
    link <- as.vector(D %*% stats::coef(fit))
    if (spec$Response_scale == "logit_share") {
      ref_lp <- as.numeric(mgcv::predict.gam(fit, newdata = rd, type = "link"))
      effect <- 100 * (plogis(ref_lp + link) - plogis(ref_lp))
      unit <- "percentage_points"
    } else {
      ref_lp <- NA_real_
      effect <- 100 * (exp(link) - 1)
      unit <- "percent_change"
    }
    data.frame(Outcome = spec$Outcome, g, Effect = effect, Link_effect = link,
      Pr_ref = plogis(ref_lp), Effect_unit = unit, stringsAsFactors = FALSE)
  })
  dplyr::bind_rows(out)
}

main_model_boot_sample <- function(d, date_levels, block_length, seed, city_count) {
  set.seed(seed)
  n_days <- length(date_levels); n_blocks <- ceiling(n_days / block_length)
  starts <- sample.int(n_days, n_blocks, replace = TRUE)
  selected <- unlist(lapply(starts, function(x) ((x - 1L + seq_len(block_length) - 1L) %% n_days) + 1L))[seq_len(n_days)]
  map <- data.frame(boot_day = seq_len(n_days), source_date = date_levels[selected],
    boot_block = rep(seq_len(n_blocks), each = block_length)[seq_len(n_days)])
  out <- merge(map, d, by.x = "source_date", by.y = "DATE", all.x = TRUE, sort = FALSE)
  if (any(table(out$boot_day) != city_count)) stop("Bootstrap sample lost city rows")
  out$DATE <- out$source_date; out$source_date <- NULL
  out <- out[order(out$city_factor, out$boot_day), , drop = FALSE]; rownames(out) <- NULL
  out$AR_start <- c(TRUE, out$city_factor[-1L] != out$city_factor[-nrow(out)] |
    out$boot_block[-1L] != out$boot_block[-nrow(out)])
  out
}

# Create a PSOCK cluster, retrying on the transient "port already in use"
# failures that happen when several independent runs start clusters at once.
main_model_make_cluster <- function(n, attempts = 20L) {
  for (k in seq_len(attempts)) {
    cl <- tryCatch(parallel::makeCluster(n, type = "PSOCK", outfile = ""),
      error = function(e) NULL)
    if (!is.null(cl)) return(cl)
    Sys.sleep(stats::runif(1, 0.5, 3))
  }
  stop("Could not create a PSOCK cluster after ", attempts, " attempts")
}

# One bootstrap replicate. Runs either in-process (n_cores = 1) or on a PSOCK
# worker; on a worker it reads the replicate context from the global
# .mm_context list that main_model_bootstrap_ewme exports (the parent session
# has already fitted bam models, so PSOCK is used instead of fork).
main_model_bootstrap_one <- function(i) {
  ctx <- .mm_context
  cache <- file.path(ctx$bundle_dir, sprintf("bootstrap_%04d.rds", i))
  if (!ctx$config$force_refit && file.exists(cache)) {
    old <- tryCatch(readRDS(cache), error = function(e) NULL)
    if (!is.null(old) && isTRUE(old$ok) &&
      (!isTRUE(ctx$config$export_bootstrap_curves) || !is.null(old$curve))) return(old)
  }
  ans <- tryCatch({
    b <- main_model_boot_sample(ctx$d, ctx$date_levels, ctx$config$block_length, ctx$seeds[i], ctx$city_count)
    # do.call forces the resampled data and AR.start vector before bam evaluates
    # its call, avoiding lazy evaluation of local objects (models3/4 pattern).
    bf <- do.call(mgcv::bam, list(formula = ctx$formula, data = b,
      family = main_model_family(ctx$spec, main_model_theta(ctx$final_fit, ctx$spec)),
      method = "fREML", discrete = TRUE, nthreads = 1, rho = ctx$final_fit$AR1.rho,
      AR.start = as.logical(b$AR_start), control = mgcv::gam.control(maxit = 200),
      sp = if (ctx$config$fix_sp) ctx$final_fit$sp else NULL))
    ew <- main_model_ewme(bf, b, ctx$spec, fixed_spec = ctx$ewme_spec); ew$Replicate <- i
    cv <- if (isTRUE(ctx$config$export_bootstrap_curves)) {
      x <- main_model_response_curves(bf, b, ctx$spec, ctx$curve_grid); x$Replicate <- i; x
    } else NULL
    list(ok = TRUE, Replicate = i, ewme = ew, curve = cv, error = NA_character_)
  }, error = function(e) list(ok = FALSE, Replicate = i, ewme = NULL, curve = NULL, error = conditionMessage(e)))
  saveRDS(ans, cache); ans
}

main_model_bootstrap_ewme <- function(final_fit, formula, d, spec, config, bundle_dir) {
  date_levels <- sort(unique(d$DATE)); city_count <- dplyr::n_distinct(d$city)
  if (length(date_levels) * city_count != nrow(d) || any(table(d$DATE) != city_count)) {
    stop("A balanced city-date panel is required for date-block bootstrap")
  }
  dir.create(bundle_dir, recursive = TRUE, showWarnings = FALSE)
  if (!is.finite(config$curve_grid_n) || config$curve_grid_n < 2L) stop("curve_grid_n must be at least 2")
  curve_grid <- if (isTRUE(config$export_bootstrap_curves)) main_model_curve_grid(d, config$curve_grid_n) else NULL
  # Fix the EWME specification on the full analysis sample and reuse it in every
  # replicate (matches the final main-model block bootstrap).
  ewme_spec <- main_model_ewme_spec(d)
  set.seed(config$seed); seeds <- sample.int(.Machine$integer.max, config$n_boot)
  n_cores <- config$n_cores
  if (is.null(n_cores) || length(n_cores) != 1L || !is.finite(n_cores) || n_cores < 1L) n_cores <- 1L
  n_cores <- min(as.integer(n_cores), config$n_boot)
  use_psock <- n_cores > 1L && .Platform$OS.type == "unix" &&
    requireNamespace("parallel", quietly = TRUE)
  # Expose one replicate's context to main_model_bootstrap_one. It is assigned in
  # the global environment so the same function works sequentially and on PSOCK
  # workers (which receive it via clusterExport).
  .mm_context <- list(d = d, formula = formula, spec = spec, config = config,
    final_fit = final_fit, seeds = seeds, date_levels = date_levels,
    city_count = city_count, bundle_dir = bundle_dir, curve_grid = curve_grid,
    ewme_spec = ewme_spec)
  if (exists(".mm_context", envir = globalenv(), inherits = FALSE)) {
    rm(".mm_context", envir = globalenv())
  }
  assign(".mm_context", .mm_context, envir = globalenv())
  on.exit(if (exists(".mm_context", envir = globalenv(), inherits = FALSE)) {
    rm(".mm_context", envir = globalenv())
  }, add = TRUE)
  job_failure <- function(i, x) list(ok = FALSE, Replicate = i, ewme = NULL, curve = NULL,
    error = if (inherits(x, "try-error")) as.character(x) else "bootstrap worker failed")
  if (use_psock) {
    cl <- main_model_make_cluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, { suppressPackageStartupMessages({ library(mgcv); library(dplyr) }); NULL })
    mm_export <- ls(envir = globalenv(), pattern = "^(main_model|MAIN_)")
    parallel::clusterExport(cl, c(mm_export, ".mm_context"), envir = globalenv())
    raw <- parallel::parLapplyLB(cl, seq_len(config$n_boot), main_model_bootstrap_one)
    jobs <- lapply(seq_len(config$n_boot), function(i) {
      x <- raw[[i]]
      if (is.list(x) && !is.null(x$Replicate)) x else job_failure(i, x)
    })
  } else {
    jobs <- lapply(seq_len(config$n_boot), main_model_bootstrap_one)
  }
  status <- bind_rows(lapply(jobs, function(x) data.frame(Outcome = spec$Outcome,
    Replicate = x$Replicate, Success = x$ok, Error = x$error)))
  success <- jobs[vapply(jobs, function(x) isTRUE(x$ok), logical(1))]
  if (!length(success)) stop("All bootstrap refits failed for ", spec$Outcome)
  reps <- bind_rows(lapply(success, function(x) x$ewme))
  point <- main_model_ewme(final_fit, d, spec, fixed_spec = ewme_spec)
  keys <- c("Outcome", "Effect", "Effect_unit")
  boot_summary <- reps %>% group_by(across(all_of(keys))) %>% summarise(
    N_boot_success = n(),
    Bootstrap_mean = mean(EWME), Bootstrap_SE = sd(EWME),
    Bootstrap_mean_link = mean(Link_EWME), Bootstrap_SE_link = sd(Link_EWME),
    .groups = "drop")
  ewme_summary <- left_join(point, boot_summary, by = keys) %>%
    mutate(
      CI_lower = if_else(
        Effect_unit == "percent_change",
        100 * (exp(Link_EWME - qnorm(.975) * Bootstrap_SE_link) - 1),
        EWME - qnorm(.975) * Bootstrap_SE
      ),
      CI_upper = if_else(
        Effect_unit == "percent_change",
        100 * (exp(Link_EWME + qnorm(.975) * Bootstrap_SE_link) - 1),
        EWME + qnorm(.975) * Bootstrap_SE
      ),
      Significant_95CI = CI_lower > 0 | CI_upper < 0,
      CI_method = if_else(
        Effect_unit == "percent_change",
        "original_link_estimate_plus_minus_1.96_calendar_date_block_bootstrap_SE",
        "original_EWME_estimate_plus_minus_1.96_calendar_date_block_bootstrap_SE"
      ),
      Bootstrap_mean_role = "diagnostic_only_not_for_reporting"
    )
  curve_summary <- NULL
  curve_replicates <- NULL
  if (isTRUE(config$export_bootstrap_curves)) {
    curve_replicates <- bind_rows(lapply(success, function(x) x$curve))
    curve_point <- main_model_response_curves(final_fit, d, spec, curve_grid)
    # Bootstrap summary keys exclude the exposure-invariant reference probability;
    # Pr_ref is taken from the full-sample point estimate below.
    curve_keys <- c("Outcome", "Exposure", "Exposure_value", "Reference_value",
                    "In_central_95pct", "Effect_unit")
    curve_boot <- curve_replicates %>% group_by(across(all_of(curve_keys))) %>% summarise(
      N_boot_success = n(), Bootstrap_SE_link = sd(Link_effect), .groups = "drop")
    # For percent_change outcomes, build the interval on the log link and
    # back-transform. For logit-share outcomes, build the interval on the logit
    # link around the point estimate's link effect and back-transform with
    # plogis at the reference probability, so the band is in percentage points
    # (never mix the link-scale SE with the percentage-point Effect scale).
    curve_summary <- left_join(curve_point, curve_boot, by = curve_keys) %>% mutate(
      # On the logit link the point effect is Link_effect; back-transform the
      # interval endpoints through plogis at the reference probability.
      CI_lower = if_else(Effect_unit == "percent_change",
        100 * (exp(Link_effect - qnorm(.975) * Bootstrap_SE_link) - 1),
        100 * (plogis(Link_effect - qnorm(.975) * Bootstrap_SE_link + qlogis(Pr_ref)) - Pr_ref)),
      CI_upper = if_else(Effect_unit == "percent_change",
        100 * (exp(Link_effect + qnorm(.975) * Bootstrap_SE_link) - 1),
        100 * (plogis(Link_effect + qnorm(.975) * Bootstrap_SE_link + qlogis(Pr_ref)) - Pr_ref)),
      Significant_95CI = CI_lower > 0 | CI_upper < 0,
      CI_method = "original_curve_plus_minus_1.96_calendar_date_block_bootstrap_SE_link"
    )
  }
  list(ewme = ewme_summary, status = status, replicates = reps,
    curves = curve_summary, curve_replicates = curve_replicates)
}

main_model_run <- function(data, outcomes, output_dir, config = main_model_defaults()) {
  main_model_require_outcomes(outcomes)
  if (!"Scale_factor" %in% names(outcomes)) outcomes$Scale_factor <- 1
  if (any(!is.finite(outcomes$Scale_factor) | outcomes$Scale_factor <= 0)) {
    stop("Scale_factor must be finite and positive")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "models"), recursive = TRUE, showWarnings = FALSE)
  scaled_data <- data
  for (i in seq_len(nrow(outcomes))) {
    scaled_data[[outcomes$Response[i]]] <- scaled_data[[outcomes$Response[i]]] / outcomes$Scale_factor[i]
  }
  d <- main_model_prepare_panel(scaled_data, outcomes$Response)
  write.csv(d, file.path(output_dir, "analysis_data_model_input.csv"), row.names = FALSE)
  write.csv(outcomes, file.path(output_dir, "outcome_specifications.csv"), row.names = FALSE)
  # Config for a single outcome's bootstrap. When several outcomes are fitted
  # concurrently, give each only a fair share of n_cores so total concurrency
  # never exceeds the cores budget (fit workers x bootstrap workers).
  n_fit_probe <- config$n_cores
  if (is.null(n_fit_probe) || length(n_fit_probe) != 1L || !is.finite(n_fit_probe) || n_fit_probe < 1L) n_fit_probe <- 1L
  n_fit_probe <- min(as.integer(n_fit_probe), nrow(outcomes))
  config_i <- config
  if (n_fit_probe > 1L) config_i$n_cores <- max(1L, as.integer(config$n_cores) %/% n_fit_probe)

  fit_one_outcome <- function(i) {
    spec <- outcomes[i, , drop = FALSE]; formula <- main_model_formula(spec$Response)
    pilot <- do.call(bam, list(formula = formula, data = d, family = main_model_family(spec),
      method = "fREML", discrete = TRUE, nthreads = 1, control = gam.control(maxit = 200)))
    rho <- main_model_rho(pilot, d, config$rho_abs_bound)
    final_fit <- do.call(bam, list(formula = formula, data = d,
      family = main_model_family(spec, main_model_theta(pilot, spec)),
      method = "fREML", discrete = TRUE, nthreads = 1, rho = rho,
      AR.start = as.logical(d$AR_start), control = gam.control(maxit = 200)))
    saveRDS(final_fit, file.path(output_dir, "models", paste0(spec$Outcome, "_T1_AR1_final.rds")))
    boot <- main_model_bootstrap_ewme(final_fit, formula, d, spec, config_i,
      file.path(output_dir, "bootstrap_bundles", spec$Outcome))
    write.csv(boot$status, file.path(output_dir, paste0(spec$Outcome, "_bootstrap_status.csv")), row.names = FALSE)
    write.csv(boot$replicates, file.path(output_dir, paste0(spec$Outcome, "_bootstrap_ewme_replicates.csv")), row.names = FALSE)
    list(ewme = boot$ewme, curves = boot$curves,
      diagnostics = data.frame(Outcome = spec$Outcome, Response = spec$Response,
        Family = spec$Family, Response_scale = spec$Response_scale, Rho = rho, N = nrow(d),
        Deviance_explained_pct = 100 * summary(final_fit)$dev.expl,
        Adjusted_R_squared = summary(final_fit)$r.sq, Converged = isTRUE(final_fit$converged)),
      status = data.frame(Outcome = spec$Outcome, Model = "T1_AR1", Success = TRUE,
        Formula = paste(deparse(formula, width.cutoff = 500), collapse = " ")))
  }
  # Outcomes are fitted independently; run them in parallel up to n_cores so the
  # pilot/AR(1) fit stage is not limited to a single thread. Each outcome still
  # uses one BAM thread, and within main_model_bootstrap_ewme the bootstrap keeps
  # using n_cores workers.
  n_fit <- config$n_cores
  if (is.null(n_fit) || length(n_fit) != 1L || !is.finite(n_fit) || n_fit < 1L) n_fit <- 1L
  n_fit <- min(as.integer(n_fit), nrow(outcomes))
  fits <- if (n_fit > 1L && .Platform$OS.type == "unix" && requireNamespace("parallel", quietly = TRUE)) {
    parallel::mclapply(seq_len(nrow(outcomes)), fit_one_outcome,
      mc.cores = n_fit, mc.preschedule = TRUE, mc.set.seed = TRUE, mc.cleanup = TRUE)
  } else {
    lapply(seq_len(nrow(outcomes)), fit_one_outcome)
  }
  fit_err <- vapply(fits, function(x) inherits(x, "try-error"), logical(1))
  if (any(fit_err)) stop("Outcome fit failed: ", as.character(fits[[which(fit_err)[1]]]))
  all_ewme <- lapply(fits, `[[`, "ewme"); all_curves <- lapply(fits, `[[`, "curves")
  all_diagnostics <- lapply(fits, `[[`, "diagnostics"); all_status <- lapply(fits, `[[`, "status")
  write.csv(bind_rows(all_ewme), file.path(output_dir, "ewme_bootstrap_CI.csv"), row.names = FALSE)
  if (isTRUE(config$export_bootstrap_curves)) {
    write.csv(bind_rows(all_curves), file.path(output_dir, "response_curves.csv"), row.names = FALSE)
  }
  write.csv(bind_rows(all_diagnostics), file.path(output_dir, "model_diagnostics.csv"), row.names = FALSE)
  write.csv(bind_rows(all_status), file.path(output_dir, "model_run_status.csv"), row.names = FALSE)
  invisible(list(data = d, outcomes = outcomes, ewme = bind_rows(all_ewme),
    curves = bind_rows(all_curves), diagnostics = bind_rows(all_diagnostics), config = config))
}

# ---------------------------------------------------------------------------
# Consolidated run report (one human-readable report.txt per output directory).
# fmt_table renders a data.frame as fixed-width text; report_block wraps it with
# a numbered section header, matching the Rain_Robustness report convention.
# ---------------------------------------------------------------------------
main_model_fmt_table <- function(df, digits = 4) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (!nrow(df)) return("(no rows)")
  cols <- lapply(df, function(x) {
    if (is.numeric(x)) {
      x <- if (all(is.na(x) | x == round(x))) formatC(x, format = "d") else formatC(x, format = "g", digits = digits)
    } else x <- as.character(x)
    x[is.na(x)] <- "NA"; x
  })
  widths <- pmax(nchar(names(df)), vapply(cols, function(x) max(nchar(x), 0L), integer(1)))
  pad <- function(x, w) formatC(x, width = w, flag = "-")
  hdr <- paste(mapply(pad, names(df), widths), collapse = "  ")
  body <- vapply(seq_len(nrow(df)), function(i)
    paste(mapply(pad, vapply(cols, `[[`, character(1), i), widths), collapse = "  "), character(1))
  c(hdr, strrep("-", nchar(hdr)), body)
}
main_model_report_block <- function(title, df, digits = 4) {
  c(strrep("-", 92), title, strrep("-", 92), main_model_fmt_table(df, digits), "")
}

# Writes report.txt into output_dir from the objects produced by main_model_run
# (result list) plus optional row/column descriptors. `scenario`, `sample_line`,
# `source_line`, `curve_ci` and `extra_blocks` let each wrapper add its context.
main_model_write_report <- function(result, output_dir, title,
                                    sample_line = NULL, source_line = NULL,
                                    exposure_line = NULL, scenario_line = NULL,
                                    extra_blocks = list(), run_version = MAIN_MODEL_VERSION,
                                    report_name = NULL) {
  if (is.null(report_name) || !nzchar(report_name)) report_name <- "report.txt"
  if (!grepl("[.]txt$", report_name)) report_name <- paste0(report_name, "_report.txt")
  cfg <- result$config; d <- result$data
  diag <- result$diagnostics; ewme <- result$ewme
  fit_report <- if (is.null(diag)) data.frame() else
    diag[, intersect(c("Outcome","Response","Family","Response_scale","Rho","N",
      "Deviance_explained_pct","Adjusted_R_squared","Converged"), names(diag)), drop = FALSE]
  ewme_report <- if (is.null(ewme) || !nrow(ewme)) data.frame() else
    ewme[, intersect(c("Outcome","Effect","Effect_unit","EWME","Link_EWME",
      "Bootstrap_SE_link","CI_lower","CI_upper","Significant_95CI","N_boot_success"), names(ewme)), drop = FALSE]
  n_cores <- cfg$n_cores; if (is.null(n_cores) || !is.finite(n_cores) || n_cores < 1) n_cores <- 1L
  header <- c(
    strrep("=", 92),
    title,
    strrep("=", 92),
    paste0("Run version   : ", run_version),
    paste0("Generated     : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    if (!is.null(source_line)) paste0("Source        : ", source_line),
    if (!is.null(exposure_line)) paste0("Exposure      : ", exposure_line),
    if (!is.null(scenario_line)) paste0("Scenario      : ", scenario_line),
    paste0("Sample        : ", if (!is.null(sample_line)) sample_line else
      paste0(nrow(d), " city-nights, ", length(unique(d$city)), " cities")),
    "Design        : frozen T1 AR(1): ti(TEMP)+ti(RH)+ti(TEMP,RH)+city-specific TEMP/RH+PRECIP+PM25+WS+month FE+city-specific date smooth(k=10)+calendar type+city RE",
    "AR(1)         : median city lag-1 working-residual rho from matching IID pilot; AR.start at city/date gaps",
    sprintf("Bootstrap     : B=%d, block=%d, fixed sp=%s, calendar-date circular moving block, normal CI on the link scale",
      cfg$n_boot, cfg$block_length, isTRUE(cfg$fix_sp)),
    paste0("CI convention : original-fit link-scale EWME +/- 1.96 * bootstrap SE; report scale per Effect_unit"),
    "", 
    main_model_report_block("1. MODEL FIT SUMMARY", fit_report, digits = 4),
    main_model_report_block("2. EWME (point estimate + bootstrap 95% CI)", ewme_report, digits = 5)
  )
  blocks <- header
  for (b in extra_blocks) blocks <- c(blocks, main_model_report_block(b$title, b$df, if (is.null(b$digits)) 4 else b$digits))
  blocks <- c(blocks, strrep("=", 92))
  report_path <- file.path(output_dir, report_name)
  writeLines(blocks, report_path)
  invisible(report_path)
}
