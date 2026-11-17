psi_weight <- function(u, c = 4.685, type = c("bisquare", "huber")) {
  type <- match.arg(type)
  au <- abs(u)
  
  if (type == "huber") {
    w <- ifelse(au <= c, 1, c / pmax(au, 1e-8))
  } else {
    w <- (1 - (u / c)^2)^2
    w[au >= c] <- 0
  }
  
  w[!is.finite(w)] <- 0
  w
}