fit_rlfpcr <- function(y, lfpca, subject, times,
                       scalar_cov = NULL,
                       max_iter = 4,
                       psi = "bisquare",
                       c = 4.685,
                       tps = NULL,
                       gram_tol = 0.05) {
  subj_levels <- unique(subject[order(subject)])
  subject_index <- match(subject, subj_levels)

  xi_long <- as.matrix(lfpca$xi[subject_index, , drop = FALSE])
  zeta_long <- as.matrix(lfpca$zeta)

  phi0 <- as.matrix(lfpca$phi.0)
  phi1 <- as.matrix(lfpca$phi.1)
  phiU <- as.matrix(lfpca$phi.U)
  D <- nrow(phi0)

  t_mean <- if (is.null(lfpca$time_mean)) mean(times) else lfpca$time_mean
  t_sd <- if (is.null(lfpca$time_sd)) sd(times) else lfpca$time_sd
  tstd <- (as.numeric(times) - t_mean) / t_sd

  df <- data.frame(
    y = as.numeric(y),
    times = as.numeric(times),
    subject = factor(subject, levels = subj_levels)
  )

  scalar_names <- character(0)
  if (!is.null(scalar_cov)) {
    Q <- as.matrix(scalar_cov)
    if (nrow(Q) != nrow(df)) stop("scalar_cov must have the same number of rows as y.")
    if (is.null(colnames(Q))) colnames(Q) <- paste0("q", seq_len(ncol(Q)))
    for (k in seq_len(ncol(Q))) df[[colnames(Q)[k]]] <- Q[, k]
    scalar_names <- colnames(Q)
  }

  NX <- ncol(xi_long)
  NU <- ncol(zeta_long)
  a_names <- paste0("xi", seq_len(NX))
  b_names <- paste0("xiT", seq_len(NX))
  cU_names <- paste0("ze", seq_len(NU))

  for (k in seq_len(NX)) {
    df[[a_names[k]]] <- xi_long[, k]
    df[[b_names[k]]] <- xi_long[, k] * tstd
  }
  for (k in seq_len(NU)) df[[cU_names[k]]] <- zeta_long[, k]

  rhs <- c("s(times, bs = 'ps', m = 1)", "s(subject, bs = 're')",
           scalar_names, a_names, b_names, cU_names)
  formula <- as.formula(paste("y ~", paste(rhs, collapse = " + ")))

  obs_w <- rep(1, nrow(df))
  for (it in seq_len(max_iter)) {
    fit <- bam(formula, data = df, weights = obs_w, method = "fREML",
               drop.unused.levels = FALSE)
    res <- residuals(fit, type = "response")
    u <- res / robust_mad(res)
    obs_w_new <- pmax(psi_weight(u, c = c, type = psi), 1e-4)
    done <- max(abs(obs_w_new - obs_w)) < 1e-5
    obs_w <- obs_w_new
    if (done) break
  }

  fit <- bam(formula, data = df, weights = obs_w, method = "fREML",
             drop.unused.levels = FALSE)
  coefs <- coef(fit)

  pick <- function(nm) {
    out <- rep(0, length(nm))
    names(out) <- nm
    hit <- intersect(names(coefs), nm)
    if (length(hit) > 0) out[hit] <- coefs[hit]
    out
  }
  delta <- pick(scalar_names)
  a <- pick(a_names)
  b <- pick(b_names)
  cU <- pick(cU_names)

  if (is.null(tps)) tps <- if (is.null(lfpca$tps)) seq_len(D) else lfpca$tps
  w <- trap_weights(tps)

  beta_hat <- gram_reconstruct(cbind(phi0, phi1, phiU), c(a, b, cU), w, gram_tol)
  beta_B <- gram_reconstruct(cbind(phi0, phi1), c(a, b), w, gram_tol)
  beta_0 <- gram_reconstruct(phi0, a, w, gram_tol)
  beta_1 <- gram_reconstruct(phi1, b, w, gram_tol)
  beta_U <- gram_reconstruct(phiU, cU, w, gram_tol)

  list(
    fit = fit,
    fitted = as.numeric(fitted(fit)),
    beta_hat = beta_hat,
    beta_B = beta_B,
    beta_0 = beta_0,
    beta_1 = beta_1,
    beta_U = beta_U,
    delta = as.numeric(delta),
    a = as.numeric(a),
    b = as.numeric(b),
    cU = as.numeric(cU),
    weights = obs_w,
    tps = tps,
    design = df
  )
}
