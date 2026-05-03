simulate_data <- function(case = c("1", "2"),
                          nI = 100,
                          nV = 4,
                          tps = seq(0.5, 99.5, by = 1),
                          sigma_x = 0.01,
                          sigma_y = 0.80,
                          sigma_b = 0.60,
                          lambda = c(2.0, 1.2, 0.7, 0.4),
                          nu = c(1.2, 0.7, 0.4),
                          scalar_cov = NULL,
                          beta_scalar = c(0.6, 0.8, -0.5, 0.2),
                          outlier = FALSE,
                          prop_curve = 0.05,
                          prop_y = 0.05,
                          prop_subject = 0.02,
                          amp_subject = 4,
                          amp_visit = 4,
                          amp_y = 5,
                          seed = NULL) {
  
  case <- match.arg(case)
  if (!is.null(seed)) set.seed(seed)
  
  phi0_basis <- function(tps, k, tmin = min(tps), tmax = max(tps)) {
    tpsc <- 2 * pi * (tps - tmin) / (tmax - tmin)
    (if (k %% 2 == 1) sin(ceiling(k / 2) * tpsc) else cos((k / 2) * tpsc)) / sqrt(tmax - tmin)
  }
  
  phi1_basis <- function(tps, k, tmin = min(tps), tmax = max(tps)) {
    tpsc <- 2 * pi * (tps - tmin) / (tmax - tmin)
    if (k == 1) {
      rep(1 / sqrt(2), length(tpsc)) / sqrt(tmax - tmin)
    } else {
      (if (k %% 2 == 1) cos((ceiling(k / 2) + 1) * tpsc) else sin((k / 2 + 2) * tpsc)) / sqrt(tmax - tmin)
    }
  }
  
  phiU_basis <- function(tps, k, tmin = min(tps), tmax = max(tps)) {
    tpsc <- 2 * (tps - tmin) / (tmax - tmin) - 1
    as.function(legendre.polynomials(k, normalized = TRUE)[[k]])(tpsc) / sqrt((tmax - tmin) / 2)
  }
  
  generate_random_times <- function(j) {
    tt <- c(0, cumsum(runif(j - 1, min = 0.3, max = 1)))
    tt - mean(tt)
  }
  
  normalize_curve <- function(v, w) {
    v <- as.numeric(v)
    nrm <- sqrt(sum(w * v^2))
    if (nrm < 1e-12) return(v)
    v / nrm
  }
  
  D <- length(tps)
  w <- trap_weights(tps)
  
  Jvec <- rep(nV, nI)
  n <- sum(Jvec)
  subject <- rep(1:nI, each = nV)
  
  times <- unlist(lapply(Jvec, FUN = generate_random_times))
  times <- times / sqrt(var(times))
  subject_index <- match(subject, 1:nI)
  
  scalar_info <- make_scalar_covariates(
    n = n,
    subject = subject,
    times = times,
    scalar_cov = scalar_cov,
    beta_scalar = beta_scalar
  )
  Q <- scalar_info$scalar_cov
  beta_scalar_true <- scalar_info$beta_scalar_true
  scalar_names <- scalar_info$scalar_names
  
  scalar_part <- if (!is.null(Q)) as.numeric(Q %*% beta_scalar_true) else rep(0, n)
  
  eta_function <- function(tps, times) {
    (times / 5 - tps / max(tps))^2 / 4
  }
  
  eta_mat_full <- t(outer(tps, times, eta_function))
  
  # CASE 1
  
  if (case == "1") {
    nX <- length(lambda)
    nU <- length(nu)
    
    phi0m <- matrix(NA, nX, D)
    phi1m <- matrix(NA, nX, D)
    phiUm <- matrix(NA, nU, D)
    
    for (k in 1:nX) {
      phi0m[k, ] <- phi0_basis(tps, k, tmin = 0, tmax = 100)
      phi1m[k, ] <- phi1_basis(tps, k, tmin = 0, tmax = 100)
    }
    for (k in 1:nU) {
      phiUm[k, ] <- phiU_basis(tps, k, tmin = 0, tmax = 100)
    }
    
    xi <- matrix(rnorm(nX * nI), nI, nX) %*% diag(sqrt(lambda), nX)
    zeta <- matrix(rnorm(nU * n), n, nU) %*% diag(sqrt(nu), nU)
    
    X0 <- xi[rep(1:nI, Jvec), , drop = FALSE] %*% phi0m
    X1 <- xi[rep(1:nI, Jvec), , drop = FALSE] %*% phi1m
    U_comp <- zeta %*% phiUm
    
    time_mat <- matrix(rep(times, D), nrow = n, ncol = D)
    B_comp <- 0.90 * (X0 + X1 * time_mat)
    U_comp <- 1.00 * U_comp
    eta_mat <- 0.05 * eta_mat_full
    
    W_latent <- eta_mat + B_comp + U_comp
    X <- W_latent + matrix(rnorm(n * D, 0, sigma_x), n, D)
    
    beta_true <- 1.10 * phi0_basis(tps, 1, 0, 100) -
      0.80 * phi0_basis(tps, 2, 0, 100) +
      0.60 * phiU_basis(tps, 1, 0, 100)
    beta_true <- normalize_curve(beta_true, w)
    
    ref <- rep(rnorm(nI, 0, sigma_b), each = nV)
    mu_true <- as.numeric(
      ref + scalar_part + 1.40 * W_latent %*% (w * beta_true)
    )
    
    y <- mu_true + rnorm(n, 0, sigma_y)
    
    truth <- list(
      beta_true = beta_true,
      beta_B_true = NULL,
      beta_U_true = NULL
    )
    
    shape_subject <- phi0_basis(tps, 1, 0, 100) -
      0.8 * phi1_basis(tps, 2, 0, 100) +
      0.7 * phiU_basis(tps, 1, 0, 100)
    
    shape_visit <- -0.9 * phiU_basis(tps, 2, 0, 100) +
      0.7 * phi0_basis(tps, 3, 0, 100)
    
    shape_subject <- normalize_curve(shape_subject, w)
    shape_visit <- normalize_curve(shape_visit, w)
  }
  
  # CASE 2
  
  if (case == "2") {
    nX <- length(lambda)
    nU <- length(nu)
    
    phi0m <- matrix(NA, nX, D)
    phi1m <- matrix(NA, nX, D)
    phiUm <- matrix(NA, nU, D)
    
    for (k in 1:nX) {
      phi0m[k, ] <- phi0_basis(tps, k, tmin = 0, tmax = 100)
      phi1m[k, ] <- phi1_basis(tps, k, tmin = 0, tmax = 100)
    }
    for (k in 1:nU) {
      phiUm[k, ] <- phiU_basis(tps, k, tmin = 0, tmax = 100)
    }
    
    xi <- matrix(rnorm(nX * nI), nI, nX) %*% diag(sqrt(lambda), nX)
    zeta <- matrix(rnorm(nU * n), n, nU) %*% diag(sqrt(nu), nU)
    
    X0 <- xi[rep(1:nI, Jvec), , drop = FALSE] %*% phi0m
    X1 <- xi[rep(1:nI, Jvec), , drop = FALSE] %*% phi1m
    U_comp <- zeta %*% phiUm
    
    time_mat <- matrix(rep(times, D), nrow = n, ncol = D)
    
    B_comp <- 1.10 * (X0 + X1 * time_mat)
    U_comp <- 0.70 * U_comp
    eta_mat <- 0.03 * eta_mat_full
    
    W_latent <- eta_mat + B_comp + U_comp
    X <- W_latent + matrix(rnorm(n * D, 0, sigma_x), n, D)
    
    beta_B_true <- 1.00 * phi0_basis(tps, 1, 0, 100) -
      0.70 * phi1_basis(tps, 1, 0, 100) +
      0.45 * phi0_basis(tps, 2, 0, 100)
    
    beta_U_true <- 1.00 * phiU_basis(tps, 1, 0, 100) -
      0.65 * phiU_basis(tps, 2, 0, 100) +
      0.45 * phiU_basis(tps, 3, 0, 100)
    
    beta_B_true <- normalize_curve(beta_B_true, w)
    beta_U_true <- normalize_curve(beta_U_true, w)
    
    ref <- rep(rnorm(nI, 0, sigma_b), each = nV)
    
    mu_true <- as.numeric(
      ref + scalar_part +
        1.35 * B_comp %*% (w * beta_B_true) +
        1.10 * U_comp %*% (w * beta_U_true)
    )
    
    y <- mu_true + rnorm(n, 0, sigma_y)
    
    truth <- list(
      beta_true = NULL,
      beta_B_true = beta_B_true,
      beta_U_true = beta_U_true
    )
    
    shape_subject <- phi0_basis(tps, 1, 0, 100) -
      0.9 * phi1_basis(tps, 2, 0, 100) +
      0.5 * phi0_basis(tps, 3, 0, 100)
    
    shape_visit <- -0.9 * phiU_basis(tps, 1, 0, 100) +
      0.7 * phiU_basis(tps, 2, 0, 100)
    
    shape_subject <- normalize_curve(shape_subject, w)
    shape_visit <- normalize_curve(shape_visit, w)
  }
  
  # STRUCTURED OUTLIERS
  
  X_clean <- X
  y_clean <- y
  
  idx_subject_out <- integer(0)
  idx_visit_out <- integer(0)
  idx_y_out <- integer(0)
  
  if (outlier) {
    x_scale <- sd(as.vector(X_clean))
    y_scale <- sd(y_clean)
    
    n_subject_out <- max(1, floor(prop_subject * nI))
    subject_out_ids <- sample(1:nI, n_subject_out)
    idx_subject_out <- which(subject %in% subject_out_ids)
    
    subj_signs <- sample(c(-1, 1), n_subject_out, replace = TRUE)
    
    for (hh in seq_along(subject_out_ids)) {
      rows_h <- which(subject == subject_out_ids[hh])
      X[rows_h, ] <- X[rows_h, ] +
        subj_signs[hh] * amp_subject * x_scale *
        matrix(shape_subject, nrow = length(rows_h), ncol = D, byrow = TRUE)
    }
    
    n_visit_out <- max(1, floor(prop_curve * n))
    remaining_rows <- setdiff(seq_len(n), idx_subject_out)
    idx_visit_out <- sample(remaining_rows, min(n_visit_out, length(remaining_rows)))
    
    visit_signs <- sample(c(-1, 1), length(idx_visit_out), replace = TRUE)
    
    for (hh in seq_along(idx_visit_out)) {
      rr <- idx_visit_out[hh]
      X[rr, ] <- X[rr, ] + visit_signs[hh] * amp_visit * x_scale * shape_visit
    }
    
    n_y_out <- max(1, floor(prop_y * n))
    idx_y_out <- sample(seq_len(n), n_y_out)
    y[idx_y_out] <- y[idx_y_out] +
      amp_y * y_scale * sign(rnorm(length(idx_y_out))) *
      (1 + abs(rt(length(idx_y_out), df = 3)))
  }
  
  list(
    case = case,
    nI = nI,
    n = n,
    Jvec = Jvec,
    subject = subject,
    times = times,
    tps = tps,
    X = X,
    y = y,
    X_clean = X_clean,
    y_clean = y_clean,
    scalar_cov = Q,
    beta_scalar_true = beta_scalar_true,
    scalar_names = scalar_names,
    mu_true = mu_true,
    W_latent = W_latent,
    B_comp = B_comp,
    U_comp = U_comp,
    eta_mat = eta_mat,
    truth = truth,
    outlier = outlier,
    idx_subject_out = idx_subject_out,
    idx_visit_out = idx_visit_out,
    idx_y_out = idx_y_out
  )
}
