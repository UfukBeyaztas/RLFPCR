## Functions for the ADNI analysis: data loading, the LPFR, LFPCR and RLFPCR
## fits, prediction for subjects left out of a fit, and the tasks of the
## leave-one-subject-out analysis and the cross-validation. RLFPCR comes from
## the package code and classical LFPCA/LFPCR from sim_functions.R of the
## simulation study, so the analysis and the simulations use the same
## estimators.

load_rlfpcr_code <- function() {
  suppressPackageStartupMessages({
    library(mgcv)
    library(splines)
    library(readxl)
    library(refund)
  })
  for (f in list.files(PKG_DIR, pattern = "[.]R$", full.names = TRUE)) {
    source(f, local = globalenv())
  }
  source(file.path(SIM_DIR, "sim_functions.R"), local = globalenv())
  invisible(TRUE)
}

make_dirs <- function() {
  for (p in c(OUT_DIR, RAW_DIR, TAB_DIR, FIG_DIR)) {
    dir.create(p, showWarnings = FALSE, recursive = TRUE)
  }
}

## Data

## Reads the ADNI file, sorts the rows by subject and visit time (required by
## rlfpca_estimate), and applies min-max scaling to each PET feature.
load_adni <- function(file = DATA_FILE) {
  dat <- as.data.frame(read_excel(file))
  subject <- as.numeric(factor(dat$patient_id))
  ord <- order(subject, dat$time_year, seq_len(nrow(dat)))
  dat <- dat[ord, ]
  subject <- subject[ord]

  feat <- as.matrix(dat[, 2:1001])
  X <- apply(feat, 2, function(z) (z - min(z)) / (max(z) - min(z)))
  Z <- cbind(Age = dat$age_visit,
             Gender = as.numeric(factor(dat$ptgender)),
             Hand = as.numeric(factor(dat$pthand)),
             Education = dat$pteducat)
  list(X = X, y = as.numeric(dat$adas13), Z = Z,
       times = as.numeric(dat$time_year), subject = subject,
       patient_id = dat$patient_id,
       site = sub("_S_.*", "", dat$patient_id),
       gender = dat$ptgender, hand = dat$pthand,
       tps = TPS, n = nrow(X))
}

subset_adni <- function(d, rows) {
  list(X = d$X[rows, , drop = FALSE], y = d$y[rows],
       Z = d$Z[rows, , drop = FALSE], times = d$times[rows],
       subject = d$subject[rows], patient_id = d$patient_id[rows],
       site = d$site[rows], gender = d$gender[rows], hand = d$hand[rows],
       tps = d$tps, n = length(rows))
}

## Site indicators (reference: first site) for the site sensitivity analysis.
site_dummies <- function(site) {
  M <- model.matrix(~ factor(site))[, -1, drop = FALSE]
  colnames(M) <- paste0("site_", sort(unique(site))[-1])
  M
}

## Coefficient functions

## Gram-corrected coefficient functions: beta_0 from the phi_k^(0) with
## coefficients a, beta_1 from the phi_k^(1) with b, beta_U from the phi_l^U
## with c, beta_B from the joint subject-level dictionary with (a, b), and
## beta from the full dictionary.
component_betas <- function(phi0, phi1, phiU, a, b, cU, w, tol = GRAM_TOL) {
  phi0 <- as.matrix(phi0)
  phi1 <- as.matrix(phi1)
  phiU <- as.matrix(phiU)
  list(beta_0 = gram_reconstruct(phi0, a, w, tol),
       beta_1 = gram_reconstruct(phi1, b, w, tol),
       beta_U = gram_reconstruct(phiU, cU, w, tol),
       beta_B = gram_reconstruct(cbind(phi0, phi1), c(a, b), w, tol),
       beta = gram_reconstruct(cbind(phi0, phi1, phiU), c(a, b, cU), w, tol))
}

## Time-specific between-subject coefficient beta_B(s, t) = beta_0 + t beta_1,
## with t on the standardized time scale of the score model.
beta_B_at <- function(obj, t_raw) {
  tt <- (t_raw - obj$time_mean) / obj$time_sd
  obj$beta_0 + tt * obj$beta_1
}

## Integral and integrated squared difference on the grid.
trapz_int <- function(f, tps = TPS) sum(trap_weights(tps) * f)
isd <- function(f, g, tps = TPS) trapz_int((f - g)^2, tps)

## Fitting

scalar_table <- function(ptab, rows, labels) {
  data.frame(term = labels,
             estimate = as.numeric(ptab[rows, 1]),
             se = as.numeric(ptab[rows, 2]),
             p = as.numeric(ptab[rows, 4]),
             stringsAsFactors = FALSE)
}

## LPFR (Goldsmith et al., 2012). refund::lpfr uses the sum over the grid
## points, so BetaHat is divided by the grid spacing to put it on the
## integral scale of the score-based coefficients.
fit_lpfr_adni <- function(d, Z = d$Z) {
  t0 <- proc.time()[3]
  lp <- lpfr(Y = d$y, funcs = d$X, subj = d$subject, covariates = Z)
  p <- ncol(Z)
  ptab <- summary(lp$fit)$p.table
  beta_sum <- as.numeric(lp$BetaHat[[1]])
  psi <- as.matrix(lp$psi[[1]])
  list(method = "LPFR",
       fitted = as.numeric(lp$fitted.vals),
       scalar = scalar_table(ptab, 1 + seq_len(p), colnames(Z)),
       beta = beta_sum / mean(diff(d$tps)),
       pred_parts = list(coef_cov = as.numeric(lp$beta.covariates),
                         meanX = colMeans(d$X),
                         proj_beta = as.numeric(psi %*% crossprod(psi, beta_sum))),
       sec = as.numeric(proc.time()[3] - t0))
}

## Prediction for subjects not used in the fit (subject effect set to zero).
predict_lpfr_new <- function(obj, Xnew, Znew) {
  pp <- obj$pred_parts
  as.numeric(pp$coef_cov[1] + as.matrix(Znew) %*% pp$coef_cov[-1] +
               sweep(Xnew, 2, pp$meanX) %*% pp$proj_beta)
}

## Classical LFPCR (Greven et al., 2010; Gertheiss et al., 2013).
fit_lfpcr_adni <- function(d, cfg, Z = d$Z, N.X = NA, N.U = NA) {
  t0 <- proc.time()[3]
  lf <- lfpca_classical(d$X, d$subject, d$times, L = cfg$L, bf = cfg$bf,
                        k_time = cfg$k_time, N.X = N.X, N.U = N.U)
  Zs <- lfpcr_design(lf, d$subject)
  reg <- score_regression(d$y, Zs, d$times, d$subject, Z)
  g <- as.numeric(reg$gamma)
  NX <- lf$N.X
  NU <- lf$N.U
  a <- g[seq_len(NX)]
  b <- g[NX + seq_len(NX)]
  cU <- g[2 * NX + seq_len(NU)]
  ptab <- summary(reg$fit)$p.table
  qn <- paste0("q", seq_len(ncol(Z)))
  out <- c(list(method = "LFPCR",
                fitted = reg$fitted,
                fitted_pop = as.numeric(predict(reg$fit, exclude = "s(subject)")),
                scalar = scalar_table(ptab, qn, colnames(Z)),
                KX = NX, KU = NU, prop_var = lf$prop_var,
                a = a, b = b, cU = cU,
                xi = lf$xi, zeta = lf$zeta,
                time_mean = lf$time_mean, time_sd = lf$time_sd),
           component_betas(lf$phi.0, lf$phi.1, lf$phi.U, a, b, cU,
                           trap_weights(d$tps)))
  out$sec <- as.numeric(proc.time()[3] - t0)
  out$heavy <- list(lf = lf, reg = reg, scalar_names = qn,
                    subject_levels = levels(reg$fit$model$subject))
  out
}

predict_lfpcr_new <- function(obj, dnew, Znew = dnew$Z) {
  lf <- obj$heavy$lf
  n <- dnew$n
  D <- ncol(dnew$X)
  tstd <- (dnew$times - lf$time_mean) / lf$time_sd
  nd <- data.frame(d.vec = rep(seq_len(D), each = n), time.vec = rep(tstd, D))
  eta <- matrix(predict(lf$mean_fit, newdata = nd), n, D)
  sc <- blup_scores(dnew$X - eta, dnew$subject, tstd, lf$phi.0, lf$phi.1,
                    lf$phi.U, lf$lambda, lf$nu, lf$sigma2)
  idx <- match(dnew$subject, unique(dnew$subject[order(dnew$subject)]))
  xi <- sc$xi[idx, , drop = FALSE]
  Zs <- cbind(xi, xi * tstd, sc$zeta)
  colnames(Zs) <- c(paste0("xi", seq_len(lf$N.X)), paste0("xiT", seq_len(lf$N.X)),
                    paste0("ze", seq_len(lf$N.U)))
  newdf <- data.frame(times = dnew$times,
                      subject = factor(rep(obj$heavy$subject_levels[1], n),
                                       levels = obj$heavy$subject_levels))
  Q <- as.matrix(Znew)
  for (k in seq_len(ncol(Q))) newdf[[obj$heavy$scalar_names[k]]] <- Q[, k]
  for (k in seq_len(ncol(Zs))) newdf[[colnames(Zs)[k]]] <- Zs[, k]
  as.numeric(predict(obj$heavy$reg$fit, newdata = newdf, exclude = "s(subject)"))
}

## RLFPCR.
fit_rlfpcr_adni <- function(d, cfg, Z = d$Z, N.X = NA, N.U = NA) {
  t0 <- proc.time()[3]
  lf <- rlfpca_estimate(Y = d$X, subject = d$subject, Time = d$times,
                        L = cfg$L, N.X = N.X, N.U = N.U, smooth = TRUE,
                        bf = cfg$bf, df_time_mean = cfg$df_time_mean,
                        mean_iter = cfg$mean_iter, clean_cells = cfg$clean_cells,
                        clean_c = cfg$clean_c, score_iter = cfg$score_iter,
                        tps = d$tps)
  fr <- fit_rlfpcr(y = d$y, lfpca = lf, subject = d$subject, times = d$times,
                   scalar_cov = Z, max_iter = cfg$max_iter, psi = cfg$psi,
                   c = cfg$c, tps = d$tps, gram_tol = GRAM_TOL)
  ptab <- summary(fr$fit)$p.table
  tot <- sum(lf$lambda) + sum(lf$nu)
  out <- c(list(method = "RLFPCR",
                fitted = fr$fitted,
                fitted_pop = as.numeric(predict(fr$fit, exclude = "s(subject)")),
                scalar = scalar_table(ptab, colnames(Z), colnames(Z)),
                KX = lf$N.X, KU = lf$N.U,
                a = fr$a, b = fr$b, cU = fr$cU,
                xi = lf$xi, zeta = lf$zeta,
                time_mean = lf$time_mean, time_sd = lf$time_sd,
                curve_w = lf$curve_w, resp_w = fr$weights),
           component_betas(lf$phi.0, lf$phi.1, lf$phi.U, fr$a, fr$b, fr$cU,
                           trap_weights(d$tps)))
  out$sec <- as.numeric(proc.time()[3] - t0)
  out$heavy <- list(lf = lf, fr = fr, cfg = cfg, Xtrain = d$X, ttrain = d$times,
                    subject_levels = levels(fr$design$subject))
  out
}

## Block-coordinate robust scores (as GetBlockwise_scores) with the visit
## times already standardized by the training mean and sd.
robust_scores_tstd <- function(Ytil, subject, tstd, phi0, phi1, phiU, lam, nu,
                               sigma2, w_curve, max_iter = 8, tol = 1e-6,
                               ridge = 1e-8) {
  n <- nrow(Ytil)
  levs <- unique(subject[order(subject)])
  idx <- match(subject, levs)
  rows <- split(seq_len(n), idx)
  NX <- ncol(phi0)
  NU <- ncol(phiU)
  A <- lapply(seq_len(n), function(r) phi0 + tstd[r] * phi1)
  UtU <- crossprod(phiU)
  PenX <- sigma2 * diag(1 / pmax(lam, 1e-8), NX)
  PenU <- sigma2 * diag(1 / pmax(nu, 1e-8), NU)
  xi <- matrix(0, length(levs), NX)
  zeta <- matrix(0, n, NU)
  update_xi <- function() {
    for (i in seq_along(levs)) {
      M <- PenX
      rhs <- rep(0, NX)
      for (r in rows[[i]]) {
        res <- Ytil[r, ] - as.numeric(phiU %*% zeta[r, ])
        M <- M + w_curve[r] * crossprod(A[[r]])
        rhs <- rhs + w_curve[r] * crossprod(A[[r]], res)
      }
      xi[i, ] <<- solve(M + ridge * diag(NX), rhs)
    }
  }
  update_xi()
  for (it in seq_len(max_iter)) {
    xi_old <- xi
    zeta_old <- zeta
    for (r in seq_len(n)) {
      res <- Ytil[r, ] - as.numeric(A[[r]] %*% xi[idx[r], ])
      zeta[r, ] <- solve(w_curve[r] * UtU + PenU + ridge * diag(NU),
                         w_curve[r] * crossprod(phiU, res))
    }
    update_xi()
    if (max(abs(xi - xi_old), abs(zeta - zeta_old)) < tol) break
  }
  list(xi = xi, zeta = zeta)
}

## Prediction for new subjects: robust mean at the new visit times, curve
## weights and cellwise winsorization calibrated on the training residuals,
## robust scores, and the robust score regression without the subject effect.
predict_rlfpcr_new <- function(obj, dnew, Znew = dnew$Z) {
  h <- obj$heavy
  lf <- h$lf
  cfg <- h$cfg
  ttr <- (h$ttrain - lf$time_mean) / lf$time_sd
  tte <- (dnew$times - lf$time_mean) / lf$time_sd
  nsb <- ns(ttr, df = cfg$df_time_mean)
  theta <- qr.solve(cbind(1, nsb), lf$eta)
  eta_new <- cbind(1, predict(nsb, tte)) %*% theta

  Rtr <- h$Xtrain - lf$eta
  Rnew <- dnew$X - eta_new
  scale_tr <- sqrt(rowMeans(Rtr^2))
  scale_new <- sqrt(rowMeans(Rnew^2))
  u <- (scale_new - median(scale_tr)) / robust_mad(scale_tr)
  w_new <- pmax(psi_weight(u, c = 4.685, type = "bisquare"), 1e-4)

  if (isTRUE(cfg$clean_cells)) {
    med <- apply(Rtr, 2, median)
    sc <- apply(Rtr, 2, robust_mad)
    lo <- med - cfg$clean_c * sc
    hi <- med + cfg$clean_c * sc
    Rnew <- pmin(pmax(Rnew, matrix(lo, nrow(Rnew), ncol(Rnew), byrow = TRUE)),
                 matrix(hi, nrow(Rnew), ncol(Rnew), byrow = TRUE))
  }

  sc <- robust_scores_tstd(Rnew, dnew$subject, tte, as.matrix(lf$phi.0),
                           as.matrix(lf$phi.1), as.matrix(lf$phi.U),
                           lf$lambda, lf$nu, lf$sigma2, w_new,
                           max_iter = cfg$score_iter)
  idx <- match(dnew$subject, unique(dnew$subject[order(dnew$subject)]))
  xi <- sc$xi[idx, , drop = FALSE]
  newdf <- data.frame(times = dnew$times,
                      subject = factor(rep(h$subject_levels[1], dnew$n),
                                       levels = h$subject_levels))
  Q <- as.matrix(Znew)
  for (k in seq_len(ncol(Q))) newdf[[colnames(Q)[k]]] <- Q[, k]
  for (k in seq_len(lf$N.X)) {
    newdf[[paste0("xi", k)]] <- xi[, k]
    newdf[[paste0("xiT", k)]] <- xi[, k] * tte
  }
  for (k in seq_len(lf$N.U)) newdf[[paste0("ze", k)]] <- sc$zeta[, k]
  as.numeric(predict(h$fr$fit, newdata = newdf, exclude = "s(subject)"))
}

## Drops the large objects before saving.
light <- function(obj) {
  obj$heavy <- NULL
  obj
}

## Leave-one-subject-out and cross-validation tasks

safe_fit <- function(expr) {
  tryCatch(expr, error = function(e) structure(list(error = conditionMessage(e)),
                                               class = "fit_error"))
}

loso_record <- function(obj, dtest) {
  if (inherits(obj, "fit_error")) return(list(success = FALSE, error = obj$error))
  list(success = TRUE, KX = obj$KX, KU = obj$KU,
       beta_0 = obj$beta_0, beta_1 = obj$beta_1, beta_U = obj$beta_U,
       beta_B = obj$beta_B, beta = obj$beta,
       time_mean = obj$time_mean, time_sd = obj$time_sd,
       delta = obj$scalar$estimate, sec = obj$sec)
}

## One deleted subject: refit LFPCR and RLFPCR without it, keep the
## coefficient functions and predict the deleted subject's visits.
run_loso_task <- function(i_subj, cfg_name, d) {
  out_dir <- file.path(RAW_DIR, "loso", cfg_name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_file <- file.path(out_dir, sprintf("subj_%03d.rds", i_subj))
  if (file.exists(out_file)) return("exists")
  if (file.exists(STOP_FILE)) return("stopped")
  cfg <- CONFIGS[[cfg_name]]
  test <- d$subject == i_subj
  dtr <- subset_adni(d, which(!test))
  dte <- subset_adni(d, which(test))
  res <- list(config = cfg_name, subject_deleted = i_subj,
              patient_id = dte$patient_id[1], y_test = dte$y)
  f1 <- safe_fit(fit_lfpcr_adni(dtr, cfg$lfpcr))
  res$LFPCR <- loso_record(f1, dte)
  if (res$LFPCR$success) {
    res$LFPCR$pred_test <- tryCatch(predict_lfpcr_new(f1, dte),
                                    error = function(e) rep(NA, dte$n))
  }
  rm(f1)
  f2 <- safe_fit(fit_rlfpcr_adni(dtr, cfg$rlfpcr))
  res$RLFPCR <- loso_record(f2, dte)
  if (res$RLFPCR$success) {
    res$RLFPCR$pred_test <- tryCatch(predict_rlfpcr_new(f2, dte),
                                     error = function(e) rep(NA, dte$n))
  }
  rm(f2)
  saveRDS(res, out_file)
  "done"
}

## Subject-level folds: all visits of a subject are in the same fold.
cv_folds <- function(subject, K, seed) {
  set.seed(seed)
  levs <- sort(unique(subject))
  f <- sample(rep(seq_len(K), length.out = length(levs)))
  f[match(subject, levs)]
}

run_cv_task <- function(r, k, d) {
  out_dir <- file.path(RAW_DIR, "cv", MAIN_CONFIG)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_file <- file.path(out_dir, sprintf("rep%02d_fold%02d.rds", r, k))
  if (file.exists(out_file)) return("exists")
  if (file.exists(STOP_FILE)) return("stopped")
  cfg <- CONFIGS[[MAIN_CONFIG]]
  fold <- cv_folds(d$subject, CV_FOLDS, CV_SEED + r)
  test <- which(fold == k)
  dtr <- subset_adni(d, which(fold != k))
  dte <- subset_adni(d, test)
  pred <- data.frame(rep = r, fold = k, row = test, subject = dte$subject,
                     y = dte$y, LPFR = NA, LFPCR = NA, RLFPCR = NA)
  info <- list()
  f0 <- safe_fit(fit_lpfr_adni(dtr))
  if (!inherits(f0, "fit_error")) pred$LPFR <- predict_lpfr_new(f0, dte$X, dte$Z)
  info$LPFR <- if (inherits(f0, "fit_error")) f0$error else "ok"
  rm(f0)
  f1 <- safe_fit(fit_lfpcr_adni(dtr, cfg$lfpcr))
  if (!inherits(f1, "fit_error")) {
    pred$LFPCR <- predict_lfpcr_new(f1, dte)
    info$LFPCR_K <- c(f1$KX, f1$KU)
  } else info$LFPCR <- f1$error
  rm(f1)
  f2 <- safe_fit(fit_rlfpcr_adni(dtr, cfg$rlfpcr))
  if (!inherits(f2, "fit_error")) {
    pred$RLFPCR <- predict_rlfpcr_new(f2, dte)
    info$RLFPCR_K <- c(f2$KX, f2$KU)
  } else info$RLFPCR <- f2$error
  rm(f2)
  saveRDS(list(pred = pred, info = info), out_file)
  "done"
}

## Runs the tasks on a PSOCK cluster. Each worker loads the code and the
## data once; the task functions skip fits that are already saved.
run_parallel <- function(tasks, fun_name, n_workers = N_WORKERS) {
  library(parallel)
  if (length(tasks) == 0) return(invisible(character(0)))
  n_workers <- min(n_workers, length(tasks))
  cl <- makeCluster(n_workers)
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl, "ADNI_DIR", envir = globalenv())
  clusterEvalQ(cl, {
    setwd(ADNI_DIR)
    source("config.R")
    source("adni_functions.R")
    load_rlfpcr_code()
    d <- load_adni()
    NULL
  })
  log_file <- file.path(RAW_DIR, paste0(fun_name, "_errors.log"))
  out <- parLapplyLB(cl, tasks, function(tk, fun_name, log_file) {
    res <- tryCatch(do.call(fun_name, c(tk, list(d = d))),
                    error = function(e) {
                      cat(format(Sys.time()), paste(unlist(tk), collapse = " "),
                          conditionMessage(e), "\n", file = log_file, append = TRUE)
                      "error"
                    })
    res
  }, fun_name = fun_name, log_file = log_file)
  invisible(unlist(out))
}

## Number formatting for the tables

fmt_p <- function(p) ifelse(p < 0.001, "$<0.001$", sprintf("%.3f", p))
fmt <- function(x, k = 3) formatC(x, format = "f", digits = k)

write_tex_table <- function(lines, file) {
  writeLines(lines, file.path(TAB_DIR, file))
}
