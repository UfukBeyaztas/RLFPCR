robust_eta <- function(Y, Time,
                       df_time = 5,
                       max_iter = 3,
                       psi = "bisquare",
                       c = 4.685,
                       ridge = 1e-8) {
  Y <- as.matrix(Y)
  n <- nrow(Y)
  time <- (Time - mean(Time)) / sqrt(var(Time))
  
  Zt <- cbind(1, ns(time, df = df_time))
  q <- ncol(Zt)
  
  w_curve <- rep(1, n)
  
  for (it in seq_len(max_iter)) {
    sqrtw <- sqrt(w_curve)
    Zw <- Zt * sqrtw
    Yw <- Y * sqrtw
    
    theta <- solve(crossprod(Zw) + ridge * diag(q), crossprod(Zw, Yw))
    eta_hat <- Zt %*% theta
    
    R <- Y - eta_hat
    curve_scale <- sqrt(rowMeans(R^2))
    u <- (curve_scale - median(curve_scale, na.rm = TRUE)) / robust_mad(curve_scale)
    w_new <- psi_weight(u, c = c, type = psi)
    w_new <- pmax(w_new, 1e-4)
    
    if (max(abs(w_new - w_curve), na.rm = TRUE) < 1e-6) {
      w_curve <- w_new
      break
    }
    w_curve <- w_new
  }
  
  sqrtw <- sqrt(w_curve)
  Zw <- Zt * sqrtw
  Yw <- Y * sqrtw
  theta <- solve(crossprod(Zw) + ridge * diag(q), crossprod(Zw, Yw))
  eta_hat <- Zt %*% theta
  
  list(
    eta = eta_hat,
    curve_w = w_curve,
    Zt = Zt,
    theta = theta
  )
}