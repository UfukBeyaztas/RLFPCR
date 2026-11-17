trap_weights <- function(tps) {
  tps <- as.numeric(tps)
  D <- length(tps)
  if (D < 2) stop("Need at least two grid points.")
  dt <- diff(tps)
  w <- numeric(D)
  w[1] <- dt[1] / 2
  w[D] <- dt[D - 1] / 2
  if (D > 2) {
    for (j in 2:(D - 1)) {
      w[j] <- (dt[j - 1] + dt[j]) / 2
    }
  }
  w
}
