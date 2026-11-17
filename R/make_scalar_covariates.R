make_scalar_covariates <- function(n, subject, times,
                                   scalar_cov = NULL,
                                   beta_scalar = c(0.6, 0.8, -0.5, 0.2)) {
  if (is.null(scalar_cov)) {
    return(list(
      scalar_cov = NULL,
      beta_scalar_true = numeric(0),
      scalar_names = character(0)
    ))
  }
  
  if (length(scalar_cov) == 1 && is.numeric(scalar_cov)) {
    p <- as.integer(scalar_cov)
    if (p < 1 || p > 4) stop("scalar_cov must be NULL or an integer between 1 and 4.")
    
    nI <- length(unique(subject))
    subject_index <- match(subject, unique(subject))
    
    q1_subj <- scale(rnorm(nI))[, 1] 
    q2_subj <- rbinom(nI, 1, 0.5)  
    
    visit_order <- ave(subject, subject, FUN = seq_along)
    q3_visit <- as.numeric(visit_order > 1)
    q4_visit <- scale(rnorm(n))[, 1]
    
    Q <- matrix(0, n, p)
    if (p >= 1) Q[, 1] <- q1_subj[subject_index]
    if (p >= 2) Q[, 2] <- q2_subj[subject_index]
    if (p >= 3) Q[, 3] <- q3_visit
    if (p >= 4) Q[, 4] <- q4_visit
    
    colnames(Q) <- paste0("q", 1:p)
    
    beta_scalar_true <- beta_scalar[1:p]
    names(beta_scalar_true) <- colnames(Q)
    
    return(list(
      scalar_cov = Q,
      beta_scalar_true = beta_scalar_true,
      scalar_names = colnames(Q)
    ))
  }
  
  Q <- as.matrix(scalar_cov)
  if (nrow(Q) != n) stop("Provided scalar_cov must have n rows.")
  
  p <- ncol(Q)
  if (p > length(beta_scalar)) {
    stop("Length of beta_scalar must be at least the number of scalar covariates.")
  }
  
  if (is.null(colnames(Q))) {
    colnames(Q) <- paste0("q", 1:p)
  }
  
  beta_scalar_true <- beta_scalar[1:p]
  names(beta_scalar_true) <- colnames(Q)
  
  list(
    scalar_cov = Q,
    beta_scalar_true = beta_scalar_true,
    scalar_names = colnames(Q)
  )
}
