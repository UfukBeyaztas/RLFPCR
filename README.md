# RLFPCR <img src="https://img.shields.io/badge/R-%3E=3.5.0-1f425f.svg" alt="R (>= 3.5.0)" align="right" height="20"/>

**RLFPCR** provides tools for fitting **robust longitudinal functional principal component regression** models for **longitudinal scalar-on-function data**. The package is designed for settings where a scalar response and a densely observed functional predictor are repeatedly measured over time for each subject, and where the data may contain atypical curves, contaminated visits, or outlying scalar responses.

The package implements a robust two-stage estimation framework. In the first stage, the longitudinal functional predictor is decomposed into a smooth mean surface, subject-specific functional intercept component, subject-specific functional slope component, and visit-specific functional deviation component using a robust longitudinal FPCA procedure. In the second stage, the resulting robust longitudinal scores are used in a robust additive mixed score regression model for the repeated scalar response.

The package includes simulation-based data generation, robust mean-surface estimation, robust covariance decomposition, robust eigenanalysis, blockwise score extraction, robust score-space regression, and coefficient-function reconstruction.

---

## 🚀 Key Features

- **Robust longitudinal scalar-on-function regression:** fits regression models where the response is scalar and the predictor is a repeatedly observed function.

- **Longitudinal functional decomposition:** decomposes the functional predictor into
  - a smooth population mean surface,
  - subject-specific functional intercept component,
  - subject-specific functional slope component,
  - visit-specific functional deviation component.

- **Robust mean-surface estimation:** estimates the longitudinal mean surface using curvewise robust weights to reduce the influence of atypical functional observations.

- **Robust covariance decomposition:** estimates subject-level and visit-level covariance components using weighted cross-product regression.

- **Robust longitudinal FPCA:** obtains subject-level and visit-level eigenfunctions and eigenvalues from robust covariance operators.

- **Blockwise score extraction:** estimates subject-level and visit-level longitudinal FPCA scores using efficient alternating blockwise updates.

- **Robust score-space regression:** fits the scalar response model using robust iteratively reweighted additive mixed regression.

- **Coefficient-function reconstruction:** reconstructs the estimated functional coefficient from the robust longitudinal FPCA basis.

- **Outlier-resistant estimation:** designed to reduce sensitivity to curvewise outliers, contaminated visits, and heavy-tailed or outlying scalar responses.

- **Built-in simulation tools:** generates clean and contaminated longitudinal scalar-on-function datasets under different functional-effect mechanisms.

- **User manual included:** the package is accompanied by the manual file `RLFPCR_1.0.0.pdf`.

---

## 📦 Installation

You can install the development version of **RLFPCR** from GitHub:

```r
install.packages("remotes")
remotes::install_github("UfukBeyaztas/RLFPCR")

Then load the package:

library(RLFPCR)

📘 Main Functions

The package contains three main user-facing functions:

Function	Description
simulate_data()	Simulates longitudinal scalar-on-function data under clean or contaminated settings.
rlfpca_estimate()	Performs robust longitudinal functional principal component analysis.
fit_rlfpcr()	Fits the robust longitudinal functional principal component regression model.

# The following shows the simulation example for Case-1 when no outliers are present in the data

sim_case1_clean <- simulate_data(
  case = "1",
  nI = 100,
  nV = 4,
  tps = seq(0.5, 99.5, by = 1),
  sigma_x = 0.01,
  sigma_y = sqrt(2),
  sigma_b = sqrt(2),
  scalar_cov = 4,
  beta_scalar = c(0.6, 0.8, -0.5, 0.2),
  outlier = FALSE
)

X <- sim_case1_clean$X
y <- sim_case1_clean$y
subject <- sim_case1_clean$subject
times <- sim_case1_clean$times
tps <- sim_case1_clean$tps
mu_true <- sim_case1_clean$mu_true
beta_true <- sim_case1_clean$truth$beta_true
scalar_cov <- sim_case1_clean$scalar_cov
beta_scalar_true <- sim_case1_clean$beta_scalar_true
w_int <- RLFPCR:::trap_weights(tps)



rlfpca_case1_clean <- rlfpca_estimate(
  Y = X,
  subject = subject,
  Time = times,
  L = 0.95,
  smooth = TRUE,
  bf = 8,
  df_time_mean = 5,
  mean_iter = 3,
  clean_cells = TRUE,
  clean_c = 4.5,
  score_iter = 8
)

rlfpcr_case1_clean <- fit_rlfpcr(
  y = y,
  lfpca = rlfpca_case1_clean,
  subject = subject,
  times = times,
  scalar_cov = scalar_cov,
  max_iter = 4,
  psi = "bisquare",
  c = 4.685
)

pred_rlfpcr <- as.numeric(predict(rlfpcr_case1_clean$fit))
beta_rlfpcr <- rlfpcr_case1_clean$beta_hat
delta_robust <- rlfpcr_case1_clean$delta
scalar_mse_rlfpcr <- mean((beta_scalar_true - delta_robust)^2)

mse_rlfpcr <- mean((mu_true - pred_rlfpcr)^2)
ise_rlfpcr <- sum((beta_rlfpcr - beta_true)^2 * w_int)

# The following shows the simulation example for Case-1 (contaminated)

sim_case1_out <- simulate_data(
  case = "1",
  nI = 100,
  nV = 4,
  tps = seq(0.5, 99.5, by = 1),
  sigma_x = 0.01,
  sigma_y = sqrt(2),
  sigma_b = sqrt(2),
  scalar_cov = 4,
  beta_scalar = c(0.6, 0.8, -0.5, 0.2),
  outlier = TRUE,
  prop_curve = 0.05,
  prop_y = 0.05,
  prop_subject = 0.05,
  amp_subject = 4,
  amp_visit = 4,
  amp_y = 5
)

X <- sim_case1_out$X
y <- sim_case1_out$y
subject <- sim_case1_out$subject
times <- sim_case1_out$times
tps <- sim_case1_out$tps
mu_true <- sim_case1_out$mu_true
beta_true <- sim_case1_out$truth$beta_true
scalar_cov <- sim_case1_out$scalar_cov
beta_scalar_true <- sim_case1_out$beta_scalar_true
w_int <- RLFPCR:::trap_weights(tps)


rlfpca_case1_out <- rlfpca_estimate(
  Y = X,
  subject = subject,
  Time = times,
  L = 0.95,
  smooth = TRUE,
  bf = 8,
  df_time_mean = 5,
  mean_iter = 3,
  clean_cells = TRUE,
  clean_c = 4.5,
  score_iter = 8
)

rlfpcr_case1_out <- fit_rlfpcr(
  y = y,
  lfpca = rlfpca_case1_out,
  subject = subject,
  times = times,
  scalar_cov = scalar_cov,
  max_iter = 4,
  psi = "bisquare",
  c = 4.685
)

pred_rlfpcr <- as.numeric(predict(rlfpcr_case1_out$fit))
beta_rlfpcr <- rlfpcr_case1_out$beta_hat
delta_robust <- rlfpcr_case1_out$delta
scalar_mse_rlfpcr <- mean((beta_scalar_true - delta_robust)^2)

mse_rlfpcr <- mean((mu_true - pred_rlfpcr)^2)
ise_rlfpcr <- sum((beta_rlfpcr - beta_true)^2 * w_int)


# The following shows the simulation example for Case-2 when no outliers are present in the data

sim_case2_clean <- simulate_data(
  case = "2",
  nI = 100,
  nV = 4,
  tps = seq(0.5, 99.5, by = 1),
  sigma_x = 0.01,
  sigma_y = 0.80,
  sigma_b = 0.60,
  lambda = c(2.0, 1.2, 0.7, 0.4),
  nu = c(1.2, 0.7, 0.4),
  scalar_cov = 4,
  beta_scalar = c(0.6, 0.8, -0.5, 0.2),
  outlier = FALSE
)

X <- sim_case2_clean$X
y <- sim_case2_clean$y
subject <- sim_case2_clean$subject
times <- sim_case2_clean$times
tps <- sim_case2_clean$tps
mu_true <- sim_case2_clean$mu_true
beta_B_true <- sim_case2_clean$truth$beta_B_true
beta_U_true <- sim_case2_clean$truth$beta_U_true
scalar_cov <- sim_case2_clean$scalar_cov
beta_scalar_true <- sim_case2_clean$beta_scalar_true
w_int <- RLFPCR:::trap_weights(tps)

rlfpca_case2_clean <- rlfpca_estimate(
  Y = X,
  subject = subject,
  Time = times,
  L = 0.95,
  smooth = TRUE,
  bf = 8,
  df_time_mean = 5,
  mean_iter = 3,
  clean_cells = TRUE,
  clean_c = 4.5,
  score_iter = 8
)

rlfpcr_case2_clean <- fit_rlfpcr(
  y = y,
  lfpca = rlfpca_case2_clean,
  subject = subject,
  times = times,
  scalar_cov = scalar_cov,
  max_iter = 4,
  psi = "bisquare",
  c = 4.685
)

pred_rlfpcr <- as.numeric(predict(rlfpcr_case2_clean$fit))
mse_rlfpcr <- mean((mu_true - pred_rlfpcr)^2)
delta_robust <- rlfpcr_case2_clean$delta
scalar_mse_rlfpcr[sim] <- mean((beta_scalar_true - delta_robust)^2)

beta_B_rlfpcr <- rep(0, nrow(rlfpca_case2_clean$phi.0))
beta_U_rlfpcr <- rep(0, nrow(rlfpca_case2_clean$phi.0))


beta_B_rlfpcr <- beta_B_rlfpcr +
  as.numeric(rlfpca_case2_clean$phi.0 %*% rlfpcr_case2_clean$a) +
  as.numeric(rlfpca_case2_clean$phi.1 %*% rlfpcr_case2_clean$b)


beta_U_rlfpcr <- as.numeric(rlfpca_case2_clean$phi.U %*% rlfpcr_case2_clean$cU)


ise_B_rlfpcr <- sum((beta_B_rlfpcr - beta_B_true)^2 * w_int)
ise_U_rlfpcr <- sum((beta_U_rlfpcr - beta_U_true)^2 * w_int)


# The following shows the simulation example for Case-2 (contaminated)

sim_case2_out <- simulate_data(
  case = "2",
  nI = 100,
  nV = 4,
  tps = seq(0.5, 99.5, by = 1),
  sigma_x = 0.01,
  sigma_y = 0.80,
  sigma_b = 0.60,
  lambda = c(2.0, 1.2, 0.7, 0.4),
  nu = c(1.2, 0.7, 0.4),
  scalar_cov = 4,
  beta_scalar = c(0.6, 0.8, -0.5, 0.2),
  outlier = TRUE,
  prop_subject = 0.05,
  prop_curve = 0.05,
  prop_y = 0.05,
  amp_subject = 4,
  amp_visit = 4,
  amp_y = 5
)

X <- sim_case2_out$X
y <- sim_case2_out$y
subject <- sim_case2_out$subject
times <- sim_case2_out$times
tps <- sim_case2_out$tps
mu_true <- sim_case2_out$mu_true
beta_B_true <- sim_case2_out$truth$beta_B_true
beta_U_true <- sim_case2_out$truth$beta_U_true
scalar_cov <- sim_case2_out$scalar_cov
beta_scalar_true <- sim_case2_out$beta_scalar_true
w_int <- RLFPCR:::trap_weights(tps)


rlfpca_case2_out <- rlfpca_estimate(
  Y = X,
  subject = subject,
  Time = times,
  L = 0.95,
  smooth = TRUE,
  bf = 8,
  df_time_mean = 5,
  mean_iter = 3,
  clean_cells = TRUE,
  clean_c = 4.5,
  score_iter = 8
)

rlfpcr_case2_out <- fit_rlfpcr(
  y = y,
  lfpca = rlfpca_case2_out,
  subject = subject,
  times = times,
  scalar_cov = scalar_cov,
  max_iter = 4,
  psi = "bisquare",
  c = 4.685
)

pred_rlfpcr <- as.numeric(predict(rlfpcr_case2_out$fit))
mse_rlfpcr <- mean((mu_true - pred_rlfpcr)^2)
delta_robust <- rlfpcr_case2_out$delta
scalar_mse_rlfpcr <- mean((beta_scalar_true - delta_robust)^2)

beta_B_rlfpcr <- rep(0, nrow(rlfpca_case2_out$phi.0))
beta_U_rlfpcr <- rep(0, nrow(rlfpca_case2_out$phi.0))

beta_B_rlfpcr <- beta_B_rlfpcr +
  as.numeric(rlfpca_case2_out$phi.0 %*% rlfpcr_case2_out$a) +
  as.numeric(rlfpca_case2_out$phi.1 %*% rlfpcr_case2_out$b)


beta_U_rlfpcr <- as.numeric(rlfpca_case2_out$phi.U %*% rlfpcr_case2_out$cU)


ise_B_rlfpcr <- sum((beta_B_rlfpcr - beta_B_true)^2 * w_int)
ise_U_rlfpcr <- sum((beta_U_rlfpcr - beta_U_true)^2 * w_int)

