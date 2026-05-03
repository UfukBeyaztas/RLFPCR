robust_mad <- function(x, constant = 1.4826) {
  med <- median(x, na.rm = TRUE)
  constant * median(abs(x - med), na.rm = TRUE) + 1e-8
}
