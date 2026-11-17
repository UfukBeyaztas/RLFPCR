GetBlockwise_scores <- function(Y.tilde, subject, Time,
                                phi.0, phi.1, phi.U,
                                lambda.hat, nu.hat, sigma2.hat,
                                w_curve = NULL,
                                max_iter = 8,
                                tol = 1e-6,
                                ridge = 1e-8) {
  Y.tilde <- as.matrix(Y.tilde)
  n <- nrow(Y.tilde)
  time <- (Time - mean(Time)) / sqrt(var(Time))
  subj_levels <- unique(subject[order(subject)])
  I <- length(subj_levels)
  subj_index <- match(subject, subj_levels)
  
  NX <- ncol(phi.0)
  NU <- ncol(phi.U)
  
  if (is.null(w_curve)) w_curve <- rep(1, n)
  w_curve <- as.numeric(w_curve)
  
  subj_rows <- split(seq_len(n), subj_index)
  
  A_list <- vector("list", n)
  AtA_list <- vector("list", n)
  AtY_list <- vector("list", n)
  
  for (r in seq_len(n)) {
    A_r <- phi.0 + time[r] * phi.1
    A_list[[r]] <- A_r
    AtA_list[[r]] <- crossprod(A_r)
    AtY_list[[r]] <- crossprod(A_r, Y.tilde[r, ])
  }
  
  UtU <- crossprod(phi.U)
  
  PenX <- sigma2.hat * diag(1 / pmax(lambda.hat, 1e-8), NX)
  PenU <- sigma2.hat * diag(1 / pmax(nu.hat, 1e-8), NU)
  
  xi <- matrix(0, I, NX)
  zeta <- matrix(0, n, NU)
  
  for (i in seq_len(I)) {
    rows_i <- subj_rows[[i]]
    M_i <- PenX
    rhs_i <- rep(0, NX)
    
    for (r in rows_i) {
      wr <- w_curve[r]
      M_i <- M_i + wr * AtA_list[[r]]
      rhs_i <- rhs_i + wr * AtY_list[[r]]
    }
    
    xi[i, ] <- solve(M_i + ridge * diag(NX), rhs_i)
  }
  
  for (it in seq_len(max_iter)) {
    xi_old <- xi
    zeta_old <- zeta
    
    for (r in seq_len(n)) {
      i <- subj_index[r]
      wr <- w_curve[r]
      resid_r <- Y.tilde[r, ] - as.numeric(A_list[[r]] %*% xi[i, ])
      M_r <- wr * UtU + PenU + ridge * diag(NU)
      rhs_r <- wr * crossprod(phi.U, resid_r)
      zeta[r, ] <- solve(M_r, rhs_r)
    }
    
    for (i in seq_len(I)) {
      rows_i <- subj_rows[[i]]
      M_i <- PenX
      rhs_i <- rep(0, NX)
      
      for (r in rows_i) {
        wr <- w_curve[r]
        resid_r <- Y.tilde[r, ] - as.numeric(phi.U %*% zeta[r, ])
        M_i <- M_i + wr * AtA_list[[r]]
        rhs_i <- rhs_i + wr * crossprod(A_list[[r]], resid_r)
      }
      
      xi[i, ] <- solve(M_i + ridge * diag(NX), rhs_i)
    }
    
    dx <- max(abs(xi - xi_old))
    dz <- max(abs(zeta - zeta_old))
    
    if (max(dx, dz) < tol) break
  }
  
  list(
    xi = xi,
    zeta = zeta
  )
}
