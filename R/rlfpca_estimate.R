rlfpca_estimate <- function(Y, subject, Time, L = 0.90,
                            N.X = NA, N.U = NA,
                            smooth = FALSE, bf = 8,
                            df_time_mean = 5,
                            mean_iter = 3,
                            clean_cells = TRUE,
                            clean_c = 4.5,
                            score_iter = 8) {
  D <- ncol(Y)
  n <- nrow(Y)
  time <- (Time - mean(Time)) / sqrt(var(Time))
  J.vec <- sapply(unique(subject[order(subject)]),
                  function(subj) sum(subject == subj))
  I <- length(J.vec)
  m <- sum(J.vec^2)
  
  t0 <- proc.time()[3]
  mean_fit <- robust_eta(
    Y = Y, Time = Time,
    df_time = df_time_mean,
    max_iter = mean_iter
  )
  eta.matrix <- mean_fit$eta
  curve_w <- mean_fit$curve_w
  t1 <- proc.time()[3]
  
  Y.tilde <- Y - eta.matrix
  if (clean_cells) Y.tilde <- clean_residuals(Y.tilde, c = clean_c)
  
  G.0 <- G.01 <- H.01 <- G.1 <- G.U <- matrix(0, D, D)
  diago <- rep(NA, D)
  
  i1 <- function(i) {
    before <- sum(J.vec[1:(i - 1)]) * (i > 1)
    rep(before + (1:J.vec[i]), each = J.vec[i])
  }
  i2 <- function(i) {
    before <- sum(J.vec[1:(i - 1)]) * (i > 1)
    rep(before + (1:J.vec[i]), J.vec[i])
  }
  ind1 <- as.vector(unlist(sapply(1:I, i1)))
  ind2 <- as.vector(unlist(sapply(1:I, i2)))
  
  X <- cbind(rep(1, m), time[ind2], time[ind1],
             time[ind1] * time[ind2], (ind1 == ind2) * 1)
  
  pair_w <- curve_w[ind1] * curve_w[ind2]
  sqrt_pair_w <- sqrt(pair_w)
  
  t2 <- proc.time()[3]
  Cmat <- matrix(unlist(lapply(1:(D - 1), function(s) {
    Y.tilde[ind1, s] * Y.tilde[ind2, (s + 1):D]
  })), length(ind1), D * (D - 1) / 2)
  
  Xw <- X * sqrt_pair_w
  Cw <- Cmat * sqrt_pair_w
  beta <- solve(crossprod(Xw) + 1e-8 * diag(ncol(Xw)), t(Xw) %*% Cw)
  
  keep <- ind1 >= ind2
  Xss <- cbind(X[, 1], X[, 2] + X[, 3], X[, 4], X[, 5])[keep, ]
  css <- sapply(1:D, function(s) {
    Y.tilde[ind1[keep], s] * Y.tilde[ind2[keep], s]
  })
  
  sqrt_pair_w_ss <- sqrt(curve_w[ind1[keep]] * curve_w[ind2[keep]])
  Xssw <- Xss * sqrt_pair_w_ss
  cssw <- css * sqrt_pair_w_ss
  betass <- solve(crossprod(Xssw) + 1e-8 * diag(ncol(Xssw)), t(Xssw) %*% cssw)
  t3 <- proc.time()[3]
  
  G.0[outer(1:D, 1:D, function(s, t) (s > t))] <- beta[1, ]
  G.0 <- G.0 + t(G.0)
  G.1[outer(1:D, 1:D, function(s, t) (s > t))] <- beta[4, ]
  G.1 <- G.1 + t(G.1)
  G.U[outer(1:D, 1:D, function(s, t) (s > t))] <- beta[5, ]
  G.U <- G.U + t(G.U)
  H.01[outer(1:D, 1:D, function(s, t) (s > t))] <- beta[2, ]
  G.01[outer(1:D, 1:D, function(s, t) (s > t))] <- beta[3, ]
  G.01 <- G.01 + t(H.01)
  
  diag(G.0) <- betass[1, ]
  diag(G.1) <- betass[3, ]
  diag(G.01) <- betass[2, ]
  diago <- betass[4, ]
  diag(G.U) <- rep(NA, D)
  
  t4 <- proc.time()[3]
  if (smooth) {
    K.0 <- smooth_sym_bam(G.0, bf = bf)
    K.1 <- smooth_sym_bam(G.1, bf = bf)
    K.01 <- smooth_asym_bam(G.01, bf = bf)
    K.U <- smooth_KU_bam(G.U, bf = bf)
  } else {
    K.U <- G.U
    K.0 <- G.0
    K.01 <- G.01
    K.1 <- G.1
    diag(K.U) <- diag_KU_bam(G.U, bf = bf)
  }
  t5 <- proc.time()[3]
  
  K.X <- rbind(cbind(K.0, K.01), cbind(t(K.01), K.1))
  sigma2.hat <- max(mean(diago - diag(K.U)), 0)
  
  eigX <- eigen(K.X, symmetric = TRUE)
  eigU <- eigen(K.U, symmetric = TRUE)
  
  lambda.hat <- eigX$values
  nu.hat <- eigU$values
  total.variance <- sum(lambda.hat * (lambda.hat > 0)) +
    sum(nu.hat * (nu.hat > 0)) + sigma2.hat
  
  if (is.na(N.X) | is.na(N.U)) {
    prop <- N.X <- N.U <- 0
    while (prop < L) {
      if (lambda.hat[N.X + 1] >= nu.hat[N.U + 1]) N.X <- N.X + 1 else N.U <- N.U + 1
      prop <- (sum(lambda.hat[1:N.X]) + sum(nu.hat[1:N.U])) / total.variance
    }
  }
  
  lambda.hat <- lambda.hat[1:N.X]
  nu.hat <- nu.hat[1:N.U]
  
  phi.X <- eigX$vectors[, 1:N.X, drop = FALSE]
  phi.U <- eigU$vectors[, 1:N.U, drop = FALSE]
  
  phi.0 <- phi.X[1:D, , drop = FALSE]
  phi.1 <- phi.X[(D + 1):(2 * D), , drop = FALSE]
  
  t6 <- proc.time()[3]
  score_fit <- GetBlockwise_scores(
    Y.tilde = Y.tilde,
    subject = subject,
    Time = Time,
    phi.0 = phi.0,
    phi.1 = phi.1,
    phi.U = phi.U,
    lambda.hat = lambda.hat,
    nu.hat = nu.hat,
    sigma2.hat = sigma2.hat,
    w_curve = curve_w,
    max_iter = score_iter,
    tol = 1e-6
  )
  t7 <- proc.time()[3]
  
  list(
    phi.0 = phi.0,
    phi.1 = phi.1,
    phi.U = phi.U,
    xi = score_fit$xi,
    zeta = score_fit$zeta,
    N.X = N.X,
    N.U = N.U,
    eta = eta.matrix,
    Y.tilde = Y.tilde,
    sigma2 = sigma2.hat,
    lambda = lambda.hat,
    nu = nu.hat,
    curve_w = curve_w
  )
}
