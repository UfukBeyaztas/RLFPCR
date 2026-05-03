clean_residuals <- function(R, c = 4.5) {
  R <- as.matrix(R)
  out <- R
  for (j in seq_len(ncol(R))) {
    medj <- median(R[, j], na.rm = TRUE)
    scj <- robust_mad(R[, j])
    lo <- medj - c * scj
    hi <- medj + c * scj
    out[, j] <- pmin(pmax(R[, j], lo), hi)
  }
  out
}
