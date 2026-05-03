diag_KU_bam <- function(GU, bf = 8) {
  D <- nrow(GU)
  idx <- which(row(GU) != col(GU))
  dat <- data.frame(
    z = GU[idx],
    row.vec = row(GU)[idx],
    col.vec = col(GU)[idx]
  )
  
  fit <- bam(
    z ~ te(row.vec, col.vec, k = bf),
    data = dat,
    method = "fREML",
    discrete = TRUE
  )
  newd <- data.frame(row.vec = seq_len(D), col.vec = seq_len(D))
  as.numeric(predict(fit, newdata = newd))
}