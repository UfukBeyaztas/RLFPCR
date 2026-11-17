simulate_data <- function(case = c("1", "2"),
                          nI = 100,
                          nV = 4,
                          tps = seq(0.5, 99.5, by = 1),
                          sigma_x = 0.01,
                          sigma_y = 0.80,
                          sigma_b = 0.60,
                          lambda = c(2.0, 1.2, 0.7, 0.4),
                          nu = c(1.2, 0.7, 0.4),
                          scalar_cov = 4,
                          beta_scalar = c(0.6, 0.8, -0.5, 0.2),
                          contamination = c("none", "mixed", "spike", "drift",
                                            "cluster", "covariate", "batch",
                                            "instrument"),
                          prop = 0.05,
                          amp_subject = 4,
                          amp_visit = 4,
                          amp_y = 5,
                          amp_spike = 3,
                          spike_width = 2,
                          amp_drift = 4,
                          amp_cluster = 4,
                          amp_cov = 5,
                          amp_batch = 4,
                          drift_rate = 2,
                          seed = NULL) {

  case <- match.arg(case)
  contamination <- match.arg(contamination)
  if (!is.null(seed)) set.seed(seed)

  tmin <- 0
  tmax <- 100
  phi0_basis <- function(k) {
    u <- 2 * pi * (tps - tmin) / (tmax - tmin)
    (if (k %% 2 == 1) sin(ceiling(k / 2) * u) else cos((k / 2) * u)) /
      sqrt(tmax - tmin)
  }
  phi1_basis <- function(k) {
    u <- 2 * pi * (tps - tmin) / (tmax - tmin)
    if (k == 1) return(rep(1 / sqrt(2), length(u)) / sqrt(tmax - tmin))
    (if (k %% 2 == 1) cos((ceiling(k / 2) + 1) * u) else sin((k / 2 + 2) * u)) /
      sqrt(tmax - tmin)
  }
  phiU_basis <- function(k) {
    x <- 2 * (tps - tmin) / (tmax - tmin) - 1
    p <- switch(k,
                rep(1, length(x)),
                x,
                (3 * x^2 - 1) / 2,
                (5 * x^3 - 3 * x) / 2)
    p * sqrt((2 * k - 1) / 2) / sqrt((tmax - tmin) / 2)
  }
  normalize_curve <- function(v) {
    nrm <- sqrt(sum(w * v^2))
    if (nrm < 1e-12) v else v / nrm
  }

  D <- length(tps)
  w <- trap_weights(tps)
  nX <- length(lambda)
  nU <- length(nu)
  if (nU > 4) stop("At most four visit-level basis functions are available.")

  Jvec <- rep(nV, nI)
  n <- sum(Jvec)
  subject <- rep(seq_len(nI), each = nV)
  times <- unlist(lapply(Jvec, function(j) c(0, cumsum(runif(j - 1, 0.3, 1)))))
  tstd <- (times - mean(times)) / sd(times)

  scalar_info <- make_scalar_covariates(n = n, subject = subject, times = times,
                                        scalar_cov = scalar_cov,
                                        beta_scalar = beta_scalar)
  Q <- scalar_info$scalar_cov
  beta_scalar_true <- scalar_info$beta_scalar_true
  scalar_part <- if (is.null(Q)) rep(0, n) else as.numeric(Q %*% beta_scalar_true)

  phi0m <- t(sapply(seq_len(nX), phi0_basis))
  phi1m <- t(sapply(seq_len(nX), phi1_basis))
  phiUm <- t(sapply(seq_len(nU), phiU_basis))
  if (nX == 1) {
    phi0m <- matrix(phi0m, 1)
    phi1m <- matrix(phi1m, 1)
  }
  if (nU == 1) phiUm <- matrix(phiUm, 1)

  xi <- matrix(rnorm(nX * nI), nI, nX) %*% diag(sqrt(lambda), nX)
  zeta <- matrix(rnorm(nU * n), n, nU) %*% diag(sqrt(nu), nU)
  xi_long <- xi[subject, , drop = FALSE]
  X0 <- xi_long %*% phi0m
  X1 <- xi_long %*% phi1m
  U0 <- zeta %*% phiUm

  c_mu <- if (case == "1") 0.05 else 0.03
  eta_mat <- c_mu * t(outer(tps, times, function(s, tt) (tt / 5 - s / 100)^2))

  if (case == "1") {
    B_comp <- 0.90 * (X0 + X1 * tstd)
    U_comp <- 1.00 * U0
  } else {
    B_comp <- 1.10 * (X0 + X1 * tstd)
    U_comp <- 0.70 * U0
  }
  W_latent <- eta_mat + B_comp + U_comp
  X <- W_latent + matrix(rnorm(n * D, 0, sigma_x), n, D)

  ref <- rep(rnorm(nI, 0, sigma_b), each = nV)

  if (case == "1") {
    beta_shape <- normalize_curve(1.10 * phi0_basis(1) - 0.80 * phi0_basis(2) +
                                    0.60 * phiU_basis(1))
    beta_true <- 1.40 * beta_shape
    mu_true <- as.numeric(ref + scalar_part + W_latent %*% (w * beta_true))
    truth <- list(beta_true = beta_true, beta_B_true = NULL, beta_U_true = NULL,
                  beta_shape = beta_shape, beta_B_shape = NULL,
                  beta_U_shape = NULL, scale = 1.40)
    shape_subject <- phi0_basis(1) - 0.8 * phi1_basis(2) + 0.7 * phiU_basis(1)
    shape_visit <- -0.9 * phiU_basis(2) + 0.7 * phi0_basis(3)
  } else {
    beta_B_shape <- normalize_curve(1.00 * phi0_basis(1) - 0.70 * phi1_basis(1) +
                                      0.45 * phi0_basis(2))
    beta_U_shape <- normalize_curve(1.00 * phiU_basis(1) - 0.65 * phiU_basis(2) +
                                      0.45 * phiU_basis(3))
    beta_B_true <- 1.35 * beta_B_shape
    beta_U_true <- 1.10 * beta_U_shape
    mu_true <- as.numeric(ref + scalar_part + B_comp %*% (w * beta_B_true) +
                            U_comp %*% (w * beta_U_true))
    truth <- list(beta_true = NULL, beta_B_true = beta_B_true,
                  beta_U_true = beta_U_true, beta_shape = NULL,
                  beta_B_shape = beta_B_shape, beta_U_shape = beta_U_shape,
                  scale = c(B = 1.35, U = 1.10))
    shape_subject <- phi0_basis(1) - 0.9 * phi1_basis(2) + 0.5 * phi0_basis(3)
    shape_visit <- -0.9 * phiU_basis(1) + 0.7 * phiU_basis(2)
  }
  shape_subject <- normalize_curve(shape_subject)
  shape_visit <- normalize_curve(shape_visit)
  shape_sys <- normalize_curve(phiU_basis(1))

  y <- mu_true + rnorm(n, 0, sigma_y)

  X_clean <- X
  y_clean <- y
  Q_clean <- Q
  idx_subject_out <- integer(0)
  idx_visit_out <- integer(0)
  idx_y_out <- integer(0)
  idx_cov_out <- integer(0)

  if (contamination != "none") {
    if (!is.null(seed)) set.seed(seed + 100000)
    x_scale <- sd(as.vector(X_clean))
    y_scale <- sd(y_clean)
    n_visit <- max(1, floor(prop * n))
    n_subj <- max(1, floor(prop * nI))
    add_row <- function(rows, curve) {
      X[rows, ] <<- X[rows, ] + matrix(curve, length(rows), D, byrow = TRUE)
    }

    if (contamination == "mixed") {
      subj_out <- sample(seq_len(nI), n_subj)
      idx_subject_out <- which(subject %in% subj_out)
      sgn <- sample(c(-1, 1), n_subj, replace = TRUE)
      for (h in seq_along(subj_out)) {
        add_row(which(subject == subj_out[h]),
                sgn[h] * amp_subject * x_scale * shape_subject)
      }
      rest <- setdiff(seq_len(n), idx_subject_out)
      idx_visit_out <- rest[sample.int(length(rest), min(n_visit, length(rest)))]
      sgn <- sample(c(-1, 1), length(idx_visit_out), replace = TRUE)
      for (h in seq_along(idx_visit_out)) {
        add_row(idx_visit_out[h], sgn[h] * amp_visit * x_scale * shape_visit)
      }
      idx_y_out <- sample(seq_len(n), n_visit)
      y[idx_y_out] <- y[idx_y_out] + amp_y * y_scale *
        sign(rnorm(n_visit)) * (1 + abs(rt(n_visit, df = 3)))
    }

    if (contamination == "spike") {
      idx_visit_out <- sample(seq_len(n), n_visit)
      lo <- min(tps) + 3 * spike_width
      hi <- max(tps) - 3 * spike_width
      for (r in idx_visit_out) {
        cen <- runif(1, lo, hi)
        bump <- exp(-0.5 * ((tps - cen) / spike_width)^2)
        add_row(r, sample(c(-1, 1), 1) * amp_spike * x_scale * bump)
      }
    }

    if (contamination == "drift") {
      subj_out <- sample(seq_len(nI), n_subj)
      idx_subject_out <- which(subject %in% subj_out)
      for (h in subj_out) {
        rows <- which(subject == h)
        frac <- times[rows] / max(times[rows])
        sg <- sample(c(-1, 1), 1)
        for (j in seq_along(rows)) {
          add_row(rows[j], sg * frac[j] * amp_drift * x_scale * shape_sys)
        }
      }
    }

    if (contamination == "cluster") {
      idx_visit_out <- sample(seq_len(n), n_visit)
      add_row(idx_visit_out, amp_cluster * x_scale * shape_visit)
    }

    if (contamination == "covariate") {
      if (is.null(Q)) stop("Covariate contamination needs scalar covariates.")
      col <- ncol(Q)
      idx_cov_out <- sample(seq_len(n), n_visit)
      Q[idx_cov_out, col] <- Q[idx_cov_out, col] +
        sample(c(-1, 1), n_visit, replace = TRUE) * amp_cov * sd(Q[, col])
    }

    if (contamination == "batch") {
      idx_visit_out <- sample(seq_len(n), n_visit)
      add_row(idx_visit_out, amp_batch * x_scale * shape_sys)
    }

    if (contamination == "instrument") {
      X <- X + drift_rate * x_scale * outer(times, shape_sys)
    }
  }

  list(
    case = case,
    contamination = contamination,
    nI = nI,
    n = n,
    Jvec = Jvec,
    subject = subject,
    times = times,
    tstd = tstd,
    tps = tps,
    X = X,
    y = y,
    X_clean = X_clean,
    y_clean = y_clean,
    scalar_cov = Q,
    scalar_cov_clean = Q_clean,
    beta_scalar_true = beta_scalar_true,
    scalar_names = scalar_info$scalar_names,
    mu_true = mu_true,
    W_latent = W_latent,
    B_comp = B_comp,
    U_comp = U_comp,
    eta_mat = eta_mat,
    truth = truth,
    idx_subject_out = idx_subject_out,
    idx_visit_out = idx_visit_out,
    idx_y_out = idx_y_out,
    idx_cov_out = idx_cov_out
  )
}
