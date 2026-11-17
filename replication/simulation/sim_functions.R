## Functions for the simulation study in Section 4 of the paper: the
## contamination settings, the competing methods and the error measures.
## RLFPCR and the data generator simulate_data() come from the package code,
## which run_simulations.R sources before this file.

scenario_table <- function() {
  data.frame(
    scenario = c("clean", "mixed05", "mixed10", "spike", "drift", "cluster",
                 "covariate", "batch20", "batch30", "instrument"),
    contamination = c("none", "mixed", "mixed", "spike", "drift", "cluster",
                      "covariate", "batch", "batch", "instrument"),
    prop = c(0, 0.05, 0.10, 0.10, 0.10, 0.10, 0.10, 0.20, 0.30, 0),
    label = c("Clean", "Mixed 5%", "Mixed 10%", "Local spikes 10%",
              "Subject drift 10%", "Clustered visits 10%",
              "Covariate outliers 10%", "Batch 20%", "Batch 30%",
              "Instrument drift"),
    stringsAsFactors = FALSE
  )
}

design_table <- function() {
  data.frame(design = 1:4, nI = c(100, 100, 250, 250), nV = c(4, 8, 4, 8))
}

## Seed of one replication. All settings of a replication share the same
## clean data; simulate_data() draws the contamination with seed + 100000.
sim_seed <- function(case, design, rep) {
  1000000 * case + 10000 * design + rep
}

## Integrated squared error on the grid (quadrature weights w).
ise_w <- function(f, g, w) sum(w * (f - g)^2)

## Projection of f onto the span of the columns of Phi.
proj_span <- function(Phi, f, w) {
  gram_reconstruct(Phi, crossprod(as.matrix(Phi), w * f), w, tol = 1e-10)
}

## Second stage of LFPCR, MAD-LFPCR, Bag-LFPCR and KFPCA-LS:
## y ~ s(time) + random subject intercept + scalar covariates + scores,
## fitted with bam() on the rows in `rows`. Fitted values are returned for
## all rows.
score_regression <- function(y, Z, times, subject, Q, rows = NULL) {
  n <- length(y)
  if (is.null(rows)) rows <- seq_len(n)
  df <- data.frame(y = y, times = times,
                   subject = factor(subject, levels = unique(subject)))
  qn <- character(0)
  if (!is.null(Q)) {
    Q <- as.matrix(Q)
    qn <- paste0("q", seq_len(ncol(Q)))
    for (k in seq_len(ncol(Q))) df[[qn[k]]] <- Q[, k]
  }
  zn <- colnames(Z)
  for (k in seq_along(zn)) df[[zn[k]]] <- Z[, k]
  fo <- as.formula(paste("y ~ s(times, bs = 'ps', m = 1) + s(subject, bs = 're')",
                         paste(c("", qn, zn), collapse = " + ")))
  fit <- bam(fo, data = df[rows, ], method = "fREML", drop.unused.levels = FALSE)
  cf <- coef(fit)
  list(fit = fit, fitted = as.numeric(predict(fit, newdata = df)),
       delta = as.numeric(cf[qn]), gamma = cf[zn])
}

## Classical LFPCA (Greven et al., 2010), the first stage of LFPCR
## (Gertheiss et al., 2013). The mean surface and the covariance functions
## are smoothed with bam(), and the scores are the BLUPs. Curves with
## keep = FALSE get weight zero. The numbers of components are chosen by the
## explained-variance threshold L unless N.X and N.U are given.

lfpca_classical <- function(Y, subject, Time, L = 0.95, bf = 10,
                            keep = rep(TRUE, nrow(Y)), N.X = NA, N.U = NA,
                            k_time = bf) {
  D <- ncol(Y)
  n <- nrow(Y)
  time_mean <- mean(Time)
  time_sd <- sd(Time)
  tstd <- (Time - time_mean) / time_sd

  Yk <- Y[keep, , drop = FALSE]
  sk <- subject[keep]
  tk <- tstd[keep]
  nk <- nrow(Yk)
  mdat <- data.frame(y = as.vector(Yk), d.vec = rep(seq_len(D), each = nk),
                     time.vec = rep(tk, D))
  mfit <- bam(y ~ te(d.vec, time.vec, k = c(bf, k_time)), data = mdat, method = "fREML",
              discrete = TRUE)
  ndat <- data.frame(d.vec = rep(seq_len(D), each = n), time.vec = rep(tstd, D))
  eta <- matrix(predict(mfit, newdata = ndat), n, D)
  Ytil <- Y - eta
  Ytk <- Ytil[keep, , drop = FALSE]

  J.vec <- as.numeric(table(factor(sk, levels = unique(sk))))
  I <- length(J.vec)
  m <- sum(J.vec^2)
  start <- c(0, cumsum(J.vec))[seq_len(I)]
  ind1 <- unlist(lapply(seq_len(I), function(i) rep(start[i] + seq_len(J.vec[i]), each = J.vec[i])))
  ind2 <- unlist(lapply(seq_len(I), function(i) rep(start[i] + seq_len(J.vec[i]), J.vec[i])))
  Xd <- cbind(1, tk[ind2], tk[ind1], tk[ind1] * tk[ind2], (ind1 == ind2) * 1)
  A <- solve(crossprod(Xd), t(Xd))
  beta <- do.call(cbind, lapply(seq_len(D - 1), function(s) {
    A %*% (Ytk[ind1, s] * Ytk[ind2, (s + 1):D, drop = FALSE])
  }))
  sel <- ind1 >= ind2
  Xss <- cbind(Xd[, 1], Xd[, 2] + Xd[, 3], Xd[, 4], Xd[, 5])[sel, ]
  css <- sapply(seq_len(D), function(s) Ytk[ind1[sel], s] * Ytk[ind2[sel], s])
  betass <- solve(crossprod(Xss), crossprod(Xss, css))

  low <- outer(seq_len(D), seq_len(D), ">")
  G.0 <- G.1 <- G.U <- G.01 <- H.01 <- matrix(0, D, D)
  G.0[low] <- beta[1, ]; G.0 <- G.0 + t(G.0)
  G.1[low] <- beta[4, ]; G.1 <- G.1 + t(G.1)
  G.U[low] <- beta[5, ]; G.U <- G.U + t(G.U)
  H.01[low] <- beta[2, ]
  G.01[low] <- beta[3, ]; G.01 <- G.01 + t(H.01)
  diag(G.0) <- betass[1, ]
  diag(G.1) <- betass[3, ]
  diag(G.01) <- betass[2, ]
  diago <- betass[4, ]
  diag(G.U) <- NA

  K.0 <- smooth_sym_bam(G.0, bf = bf)
  K.1 <- smooth_sym_bam(G.1, bf = bf)
  K.01 <- smooth_asym_bam(G.01, bf = bf)
  K.U <- smooth_KU_bam(G.U, bf = bf)
  K.X <- rbind(cbind(K.0, K.01), cbind(t(K.01), K.1))
  sigma2 <- max(mean(diago - diag(K.U)), 0)
  eX <- eigen(K.X, symmetric = TRUE)
  eU <- eigen(K.U, symmetric = TRUE)
  lam <- eX$values
  nu <- eU$values
  tot <- sum(lam[lam > 0]) + sum(nu[nu > 0]) + sigma2
  if (is.na(N.X) || is.na(N.U)) {
    NX <- NU <- 0
    prop <- 0
    while (prop < L) {
      if (lam[NX + 1] >= nu[NU + 1]) NX <- NX + 1 else NU <- NU + 1
      prop <- (sum(lam[seq_len(NX)]) + sum(nu[seq_len(NU)])) / tot
    }
  } else {
    NX <- N.X
    NU <- N.U
  }
  phiX <- eX$vectors[, seq_len(NX), drop = FALSE]
  phiU <- eU$vectors[, seq_len(NU), drop = FALSE]
  phi0 <- phiX[seq_len(D), , drop = FALSE]
  phi1 <- phiX[D + seq_len(D), , drop = FALSE]

  sc <- blup_scores(Ytil, subject, tstd, phi0, phi1, phiU,
                    lam[seq_len(NX)], nu[seq_len(NU)], sigma2,
                    w_curve = as.numeric(keep))
  list(phi.0 = phi0, phi.1 = phi1, phi.U = phiU, xi = sc$xi, zeta = sc$zeta,
       N.X = NX, N.U = NU, tstd = tstd,
       lambda = lam[seq_len(NX)], nu = nu[seq_len(NU)], sigma2 = sigma2,
       prop_var = (sum(lam[seq_len(NX)]) + sum(nu[seq_len(NU)])) / tot,
       mean_fit = mfit, time_mean = time_mean, time_sd = time_sd)
}

## BLUP scores with curve weights w_curve, computed exactly for each subject
## after profiling out the visit scores.
blup_scores <- function(Ytil, subject, tstd, phi0, phi1, phiU, lam, nu, sigma2,
                        w_curve = rep(1, nrow(Ytil)), ridge = 1e-8) {
  levs <- unique(subject[order(subject)])
  idx <- match(subject, levs)
  NX <- ncol(phi0)
  NU <- ncol(phiU)
  PenX <- sigma2 * diag(1 / pmax(lam, 1e-8), NX) + ridge * diag(NX)
  PenU <- sigma2 * diag(1 / pmax(nu, 1e-8), NU) + ridge * diag(NU)
  UtU <- crossprod(phiU)
  xi <- matrix(0, length(levs), NX)
  zeta <- matrix(0, nrow(Ytil), NU)
  for (i in seq_along(levs)) {
    rows <- which(idx == i)
    M <- PenX
    rhs <- rep(0, NX)
    keepM <- vector("list", length(rows))
    for (j in seq_along(rows)) {
      r <- rows[j]
      wr <- w_curve[r]
      Ar <- phi0 + tstd[r] * phi1
      Minv <- solve(wr * UtU + PenU)
      AU <- crossprod(Ar, phiU)
      yr <- Ytil[r, ]
      M <- M + wr * crossprod(Ar) - wr^2 * AU %*% Minv %*% t(AU)
      rhs <- rhs + wr * crossprod(Ar, yr) - wr^2 * AU %*% (Minv %*% crossprod(phiU, yr))
      keepM[[j]] <- list(Ar = Ar, Minv = Minv, yr = yr, wr = wr)
    }
    xi[i, ] <- solve(M, rhs)
    for (j in seq_along(rows)) {
      k <- keepM[[j]]
      zeta[rows[j], ] <- k$wr * k$Minv %*% crossprod(phiU, k$yr - k$Ar %*% xi[i, ])
    }
  }
  list(xi = xi, zeta = zeta)
}

## Score covariates of the second stage: xi, xi times standardized time,
## and zeta.
lfpcr_design <- function(lf, subject) {
  idx <- match(subject, unique(subject[order(subject)]))
  xi <- lf$xi[idx, , drop = FALSE]
  Z <- cbind(xi, xi * lf$tstd, lf$zeta)
  colnames(Z) <- c(paste0("xi", seq_len(lf$N.X)), paste0("xiT", seq_len(lf$N.X)),
                   paste0("ze", seq_len(lf$N.U)))
  Z
}

## Fitted values, scalar effects and Gram-corrected coefficient functions.
lfpcr_output <- function(lf, reg, w) {
  NX <- lf$N.X
  NU <- lf$N.U
  g <- as.numeric(reg$gamma)
  a <- g[seq_len(NX)]
  b <- g[NX + seq_len(NX)]
  cU <- g[2 * NX + seq_len(NU)]
  PhiB <- cbind(lf$phi.0, lf$phi.1)
  Phi <- cbind(PhiB, lf$phi.U)
  list(fitted = reg$fitted, delta = reg$delta,
       beta = gram_reconstruct(Phi, c(a, b, cU), w),
       beta_B = gram_reconstruct(PhiB, c(a, b), w),
       beta_U = gram_reconstruct(lf$phi.U, cU, w),
       span = list(full = Phi, B = PhiB, U = lf$phi.U),
       KX = NX, KU = NU)
}

## The seven methods. Each takes a simulated data set d and quadrature
## weights w, and returns fitted values, scalar effects and coefficient
## functions.

m_lpfr <- function(d, w) {
  fit <- lpfr(Y = d$y, subj = d$subject, covariates = d$scalar_cov, funcs = d$X)
  bh <- if (!is.null(fit$betaHat)) fit$betaHat[[1]] else fit$BetaHat[[1]]
  list(fitted = as.numeric(fit$fitted.vals),
       delta = as.numeric(fit$beta.covariates)[-1],
       beta = as.numeric(bh), beta_B = NULL, beta_U = NULL,
       span = NULL, KX = NA, KU = NA)
}

m_lfpcr <- function(d, w, bf = 10) {
  lf <- lfpca_classical(d$X, d$subject, d$times, L = 0.95, bf = bf)
  reg <- score_regression(d$y, lfpcr_design(lf, d$subject), d$times,
                          d$subject, d$scalar_cov)
  lfpcr_output(lf, reg, w)
}

## Curve screening of MAD-LFPCR: robust z-score of the L2 distance to the
## pointwise median.
flag_curves_mad <- function(X, w, cut = 3) {
  med <- apply(X, 2, median)
  dist <- sqrt(as.numeric((sweep(X, 2, med)^2) %*% w))
  z <- (dist - median(dist)) / mad(dist)
  z > cut
}

## Curve screening of Bag-LFPCR: curves outside the functional bagplot
## (Hyndman and Shang, 2010).
flag_curves_bagplot <- function(X, tps) {
  Yt <- t(X)
  colnames(Yt) <- seq_len(nrow(X))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  fb <- suppressWarnings(rainbow::fboxplot(rainbow::fds(x = tps, y = Yt),
                                           plot.type = "functional", type = "bag",
                                           projmethod = "PCAproj", factor = 1.96))
  out <- rep(FALSE, nrow(X))
  if (length(fb) > 0) out[as.integer(fb)] <- TRUE
  out
}

## Response screening: residuals more than `cut` MADs from their median.
flag_resid_mad <- function(r, cut = 3) {
  abs(r - median(r)) / mad(r) > cut
}

## MAD-LFPCR and Bag-LFPCR: drop the flagged curves, fit LFPCR, drop the
## visits with outlying residuals and refit the second stage.
m_screen_lfpcr <- function(d, w, rule = c("mad", "bagplot"), bf = 10) {
  rule <- match.arg(rule)
  drop_curve <- if (rule == "mad") flag_curves_mad(d$X, w) else
    flag_curves_bagplot(d$X, d$tps)
  keep <- !drop_curve
  lf <- lfpca_classical(d$X, d$subject, d$times, L = 0.95, bf = bf, keep = keep)
  Z <- lfpcr_design(lf, d$subject)
  reg <- score_regression(d$y, Z, d$times, d$subject, d$scalar_cov,
                          rows = which(keep))
  drop_y <- rep(FALSE, length(d$y))
  drop_y[keep] <- flag_resid_mad(d$y[keep] - reg$fitted[keep])
  reg <- score_regression(d$y, Z, d$times, d$subject, d$scalar_cov,
                          rows = which(keep & !drop_y))
  out <- lfpcr_output(lf, reg, w)
  out$n_curve_drop <- sum(drop_curve)
  out$n_resp_drop <- sum(drop_y)
  out
}

## Kendall FPCA (Zhong et al., 2022) for dense curves. The eigenfunctions
## are those of the U-statistic estimate of the Kendall tau function
## K(s, s') = E[(X - X*)(s) (X - X*)(s') / ||X - X*||^2], the scores are
## projections of the curves centred at the pointwise median, and the number
## of components is chosen from the MAD variances of the scores.
kendall_fpca <- function(X, w, L = 0.95, kmax = 20) {
  sw <- sqrt(w)
  Xw <- sweep(X, 2, sw, "*")
  sq <- rowSums(Xw^2)
  D2 <- outer(sq, sq, "+") - 2 * tcrossprod(Xw)
  diag(D2) <- Inf
  Wm <- 1 / pmax(D2, 1e-12)
  diag(Wm) <- 0
  K <- crossprod(Xw, rowSums(Wm) * Xw) - crossprod(Xw, Wm %*% Xw)
  e <- eigen((K + t(K)) / 2, symmetric = TRUE)
  kmax <- min(kmax, ncol(X))
  phi <- sweep(e$vectors[, seq_len(kmax), drop = FALSE], 1, sw, "/")
  med <- apply(X, 2, median)
  sc <- sweep(X, 2, med) %*% (w * phi)
  rv <- apply(sc, 2, mad)^2
  K_sel <- which(cumsum(rv) / sum(rv) >= L)[1]
  list(phi = phi[, seq_len(K_sel), drop = FALSE],
       scores = sc[, seq_len(K_sel), drop = FALSE], K = K_sel)
}

m_kfpca_ls <- function(d, w) {
  kf <- kendall_fpca(d$X, w)
  Z <- kf$scores
  colnames(Z) <- paste0("sc", seq_len(kf$K))
  reg <- score_regression(d$y, Z, d$times, d$subject, d$scalar_cov)
  list(fitted = reg$fitted, delta = reg$delta,
       beta = gram_reconstruct(kf$phi, as.numeric(reg$gamma), w),
       beta_B = NULL, beta_U = NULL, span = list(full = kf$phi),
       KX = kf$K, KU = NA)
}

## rLMM: robustlmm::rlmer (Koller, 2016) with the curves reduced to G bin
## integrals, y ~ T + T^2 + q + z_1 + ... + z_G + (1 | subject). The
## coefficient of z_g estimates beta(s) on bin g.
m_rlmm_grid <- function(d, w, G = 10) {
  D <- ncol(d$X)
  bin <- cut(seq_len(D), G, labels = FALSE)
  Zb <- sapply(seq_len(G), function(g) as.numeric(d$X[, bin == g, drop = FALSE] %*%
                                                    w[bin == g]))
  colnames(Zb) <- paste0("z", seq_len(G))
  Q <- as.matrix(d$scalar_cov)
  colnames(Q) <- paste0("q", seq_len(ncol(Q)))
  df <- data.frame(y = d$y, T1 = d$times, T2 = d$times^2,
                   subject = factor(d$subject), Q, Zb)
  fo <- as.formula(paste("y ~ T1 + T2 +", paste(c(colnames(Q), colnames(Zb)),
                                                collapse = " + "), "+ (1 | subject)"))
  fit <- robustlmm::rlmer(fo, data = df)
  fe <- lme4::fixef(fit)
  ind <- sapply(seq_len(G), function(g) as.numeric(bin == g))
  list(fitted = as.numeric(fitted(fit)), delta = as.numeric(fe[colnames(Q)]),
       beta = as.numeric(ind %*% fe[colnames(Zb)]), beta_B = NULL, beta_U = NULL,
       span = list(full = ind), KX = NA, KU = NA)
}

## RLFPCR with the settings of the paper. The sensitivity study changes c,
## clean_cells and clean_c.
m_rlfpcr <- function(d, w, bf = 8, c = 4.685, clean_cells = TRUE, clean_c = 4.5) {
  lf <- rlfpca_estimate(Y = d$X, subject = d$subject, Time = d$times, L = 0.95,
                        smooth = TRUE, bf = bf, df_time_mean = 5, mean_iter = 3,
                        clean_cells = clean_cells, clean_c = clean_c, score_iter = 8,
                        c = c, tps = d$tps)
  fit <- fit_rlfpcr(y = d$y, lfpca = lf, subject = d$subject, times = d$times,
                    scalar_cov = d$scalar_cov, max_iter = 4, psi = "bisquare",
                    c = c, tps = d$tps)
  PhiB <- cbind(lf$phi.0, lf$phi.1)
  list(fitted = fit$fitted, delta = fit$delta, beta = fit$beta_hat,
       beta_B = fit$beta_B, beta_U = fit$beta_U,
       span = list(full = cbind(PhiB, lf$phi.U), B = PhiB, U = lf$phi.U),
       KX = lf$N.X, KU = lf$N.U,
       n_curve_drop = sum(lf$curve_w < 0.5), n_resp_drop = sum(fit$weights < 0.5))
}

method_list <- function() {
  list(
    "LPFR" = m_lpfr,
    "LFPCR" = m_lfpcr,
    "MAD-LFPCR" = function(d, w) m_screen_lfpcr(d, w, "mad"),
    "Bag-LFPCR" = function(d, w) m_screen_lfpcr(d, w, "bagplot"),
    "KFPCA-LS" = m_kfpca_ls,
    "rLMM" = m_rlmm_grid,
    "RLFPCR" = m_rlfpcr
  )
}

## Errors of one fit: MSE against the true conditional mean, MSE_delta for
## the scalar effects, and ISEs against the coefficient functions that
## generate the response. The *_proj versions use the projections of these
## functions onto the span retained by the method.
evaluate_fit <- function(out, d, w) {
  tr <- d$truth
  res <- c(mse = mean((out$fitted - d$mu_true)^2),
           mse_delta = mean((out$delta - d$beta_scalar_true)^2),
           ise = NA, ise_proj = NA, ise_B = NA, ise_B_proj = NA,
           ise_U = NA, ise_U_proj = NA)
  if (d$case == "1") {
    res["ise"] <- ise_w(out$beta, tr$beta_true, w)
    if (!is.null(out$span$full)) {
      res["ise_proj"] <- ise_w(out$beta, proj_span(out$span$full, tr$beta_true, w), w)
    }
  } else if (!is.null(out$beta_B)) {
    res["ise_B"] <- ise_w(out$beta_B, tr$beta_B_true, w)
    res["ise_U"] <- ise_w(out$beta_U, tr$beta_U_true, w)
    res["ise_B_proj"] <- ise_w(out$beta_B, proj_span(out$span$B, tr$beta_B_true, w), w)
    res["ise_U_proj"] <- ise_w(out$beta_U, proj_span(out$span$U, tr$beta_U_true, w), w)
  }
  res
}

## One replication: simulate the data and fit the methods. A failed fit
## gives a row of NA with the error message.
run_one <- function(case, design, scenario, rep, methods = names(method_list()),
                    keep_curves = FALSE) {
  dt <- design_table()[design, ]
  sc <- scenario_table()
  sc <- sc[sc$scenario == scenario, ]
  d <- simulate_data(case = as.character(case), nI = dt$nI, nV = dt$nV,
                     contamination = sc$contamination,
                     prop = max(sc$prop, 0.01),
                     seed = sim_seed(case, design, rep))
  w <- trap_weights(d$tps)
  ml <- method_list()[methods]
  rows <- list()
  curves <- list()
  for (mn in names(ml)) {
    t0 <- proc.time()[3]
    out <- tryCatch(ml[[mn]](d, w), error = function(e) e)
    el <- proc.time()[3] - t0
    base <- data.frame(case = case, nI = dt$nI, nV = dt$nV, scenario = scenario,
                       rep = rep, method = mn, stringsAsFactors = FALSE)
    if (inherits(out, "error")) {
      met <- as.list(rep(NA_real_, 8))
      names(met) <- c("mse", "mse_delta", "ise", "ise_proj", "ise_B",
                      "ise_B_proj", "ise_U", "ise_U_proj")
      rows[[mn]] <- cbind(base, as.data.frame(met), KX = NA, KU = NA,
                          n_curve_drop = NA, n_resp_drop = NA, time = el,
                          error = conditionMessage(out), stringsAsFactors = FALSE)
    } else {
      met <- evaluate_fit(out, d, w)
      rows[[mn]] <- cbind(base, as.data.frame(as.list(met)),
                          KX = out$KX, KU = out$KU,
                          n_curve_drop = if (is.null(out$n_curve_drop)) NA else out$n_curve_drop,
                          n_resp_drop = if (is.null(out$n_resp_drop)) NA else out$n_resp_drop,
                          time = el, error = "", stringsAsFactors = FALSE)
      if (keep_curves) {
        curves[[mn]] <- list(beta = out$beta, beta_B = out$beta_B,
                             beta_U = out$beta_U)
      }
    }
  }
  list(metrics = do.call(rbind, rows), curves = curves,
       truth = if (keep_curves) d$truth else NULL,
       tps = if (keep_curves) d$tps else NULL)
}
