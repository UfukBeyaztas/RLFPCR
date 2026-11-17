gram_reconstruct <- function(Phi, coef, w, tol = 0.05) {
  Phi <- as.matrix(Phi)
  if (ncol(Phi) == 0) return(rep(0, nrow(Phi)))
  G <- crossprod(Phi, w * Phi)
  e <- eigen(G, symmetric = TRUE)
  keep <- e$values > tol * max(e$values)
  V <- e$vectors[, keep, drop = FALSE]
  d <- V %*% (crossprod(V, as.numeric(coef)) / e$values[keep])
  as.numeric(Phi %*% d)
}
