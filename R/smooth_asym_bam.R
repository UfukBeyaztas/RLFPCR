smooth_asym_bam <- function(G, bf = 8) {
  D <- nrow(G)
  row.vec <- rep(seq_len(D), each = D)
  col.vec <- rep(seq_len(D), D)
  dat <- data.frame(z = as.vector(G), row.vec = row.vec, col.vec = col.vec)
  
  fit <- bam(
    z ~ te(row.vec, col.vec, k = bf),
    data = dat,
    method = "fREML",
    discrete = TRUE
  )
  matrix(as.numeric(predict(fit, newdata = dat)), D, D)
}
