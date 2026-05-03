fit_rlfpcr <- function(y, lfpca, subject, times,
                       scalar_cov = NULL,
                       max_iter = 4,
                       psi = "bisquare",
                       c = 4.685) {
  subject_index <- match(subject, unique(subject[order(subject)]))
  
  xi_long <- as.matrix(lfpca$xi[subject_index, , drop = FALSE])
  zeta_long <- as.matrix(lfpca$zeta)
  
  phi0 <- as.matrix(lfpca$phi.0)
  phi1 <- as.matrix(lfpca$phi.1)
  phiU <- as.matrix(lfpca$phi.U)
  
  df <- data.frame(
    y = as.numeric(y),
    times = as.numeric(times),
    subject = as.factor(subject)
  )
  
  scalar_names <- character(0)
  if (!is.null(scalar_cov)) {
    Q <- as.matrix(scalar_cov)
    if (nrow(Q) != nrow(df)) stop("scalar_cov must have the same number of rows as y.")
    if (is.null(colnames(Q))) colnames(Q) <- paste0("q", 1:ncol(Q))
    
    for (k in 1:ncol(Q)) {
      df[[colnames(Q)[k]]] <- Q[, k]
    }
    scalar_names <- colnames(Q)
  }
  
  NX <- ncol(xi_long)
  NU <- ncol(zeta_long)
  
  if (NX > 0) {
    for (k in 1:NX) {
      df[[paste0("xi", k)]] <- xi_long[, k]
      df[[paste0("xiT", k)]] <- xi_long[, k] * df$times
    }
  }
  
  if (NU > 0) {
    for (k in 1:NU) df[[paste0("ze", k)]] <- zeta_long[, k]
  }
  
  rhs <- c("s(times, bs='ps', m=1)", "s(subject, bs='re')")
  if (length(scalar_names) > 0) rhs <- c(rhs, scalar_names)
  if (NX > 0) rhs <- c(rhs, paste0("xi", 1:NX), paste0("xiT", 1:NX))
  if (NU > 0) rhs <- c(rhs, paste0("ze", 1:NU))
  
  formula <- as.formula(paste("y ~", paste(rhs, collapse = " + ")))
  
  obs_w <- rep(1, nrow(df))
  for (it in 1:max_iter) {
    fit <- gam(formula, data = df, weights = obs_w, method = "REML")
    res <- residuals(fit, type = "response")
    s <- robust_mad(res)
    u <- res / s
    obs_w_new <- psi_weight(u, c = c, type = psi)
    obs_w_new <- pmax(obs_w_new, 1e-4)
    
    if (max(abs(obs_w_new - obs_w)) < 1e-5) {
      obs_w <- obs_w_new
      break
    }
    obs_w <- obs_w_new
  }
  
  fit <- gam(formula, data = df, weights = obs_w, method = "REML")
  coefs <- coef(fit)
  
  delta <- rep(0, length(scalar_names))
  names(delta) <- scalar_names
  hit_delta <- intersect(names(coefs), scalar_names)
  if (length(hit_delta) > 0) delta[match(hit_delta, scalar_names)] <- coefs[hit_delta]
  
  a_names  <- paste0("xi", 1:NX)
  b_names  <- paste0("xiT", 1:NX)
  cU_names <- paste0("ze", 1:NU)
  
  a <- rep(0, NX)
  b <- rep(0, NX)
  cU <- rep(0, NU)
  
  names(a) <- a_names
  names(b) <- b_names
  names(cU) <- cU_names
  
  hit_a <- intersect(names(coefs), a_names)
  hit_b <- intersect(names(coefs), b_names)
  hit_c <- intersect(names(coefs), cU_names)
  
  if (length(hit_a) > 0) a[match(hit_a, a_names)] <- coefs[hit_a]
  if (length(hit_b) > 0) b[match(hit_b, b_names)] <- coefs[hit_b]
  if (length(hit_c) > 0) cU[match(hit_c, cU_names)] <- coefs[hit_c]
  
  beta_hat <- rep(0, nrow(phi0))
  if (NX > 0) {
    beta_hat <- beta_hat +
      as.numeric(phi0 %*% matrix(a, ncol = 1)) +
      as.numeric(phi1 %*% matrix(b, ncol = 1))
  }
  if (NU > 0) {
    beta_hat <- beta_hat + as.numeric(phiU %*% matrix(cU, ncol = 1))
  }
  
  list(
    fit = fit,
    beta_hat = beta_hat,
    delta = as.numeric(delta),
    a = as.numeric(a),
    b = as.numeric(b),
    cU = as.numeric(cU),
    weights = obs_w,
    design = df
  )
}
