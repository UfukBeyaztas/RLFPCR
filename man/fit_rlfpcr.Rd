\name{fit_rlfpcr}
\alias{fit_rlfpcr}
\title{
Fit Robust Longitudinal Functional Principal Component Regression
}
\description{
Fits the second-stage robust longitudinal functional principal component
regression model. The function uses robust longitudinal FPCA scores obtained
from \code{\link{rlfpca_estimate}}, fits a robust additive mixed score
regression for a repeated scalar response, and reconstructs the functional
coefficients through the Gram system of the estimated eigenfunctions.
}
\usage{
fit_rlfpcr(y, lfpca, subject, times,
           scalar_cov = NULL,
           max_iter = 4,
           psi = "bisquare",
           c = 4.685,
           tps = NULL,
           gram_tol = 0.05)
}
\arguments{
  \item{y}{
Numeric vector of repeated scalar responses, one value per subject-visit
observation.
}
  \item{lfpca}{
An object returned by \code{\link{rlfpca_estimate}}.
}
  \item{subject}{
Vector of subject identifiers of the same length as \code{y}.
}
  \item{times}{
Numeric vector of visit times on the original scale, the same values that were
passed to \code{\link{rlfpca_estimate}}. The smooth time effect uses these
values; the score interactions use the standardized times
\eqn{t_{ij}=(T_{ij}-\bar T)/s_T} computed with the mean and standard deviation
stored in \code{lfpca}.
}
  \item{scalar_cov}{
Optional matrix or data frame of scalar explanatory variables with one row per
observation. Missing column names are set to \code{q1}, \code{q2}, ...
}
  \item{max_iter}{
Maximum number of iteratively reweighted least-squares iterations. The default
is \code{4}.
}
  \item{psi}{
Weight function for the response residuals, \code{"bisquare"} or
\code{"huber"}.
}
  \item{c}{
Tuning constant of the weight function. The default \code{4.685} is the usual
bisquare constant.
}
  \item{tps}{
Grid of the functional predictor, used for the quadrature weights of the Gram
matrices. If \code{NULL}, the grid stored in \code{lfpca} is used.
}
  \item{gram_tol}{
Relative tolerance of the pseudo-inverse of the Gram matrix: eigenvalues of
\eqn{\widehat G} smaller than \code{gram_tol} times the largest eigenvalue are
treated as zero. The default is \code{0.05}.
}
}
\details{
Substituting the truncated decomposition of \code{\link{rlfpca_estimate}} into
the scalar-on-function model gives the score regression

\deqn{
Y_{ij} = \alpha(T_{ij}) + q_{ij}^{\top}\delta + u_i
+ \sum_{k=1}^{K_X} a_k \xi_{ik} + \sum_{k=1}^{K_X} b_k t_{ij}\xi_{ik}
+ \sum_{\ell=1}^{K_U} c_{\ell}\zeta_{ij\ell} + e_{ij},
}

with \eqn{a_k=\langle\beta_0,\phi_k^{(0)}\rangle},
\eqn{b_k=\langle\beta_1,\phi_k^{(1)}\rangle} and
\eqn{c_\ell=\langle\beta_U,\phi_\ell^U\rangle}. The smooth time effect is
\code{s(times, bs = "ps", m = 1)} and the subject random intercept is
\code{s(subject, bs = "re")}. The model is fitted with \code{bam} (fast REML)
and observation weights \eqn{w_{ij}=\psi(r_{ij}/\hat\tau)/(r_{ij}/\hat\tau)},
where \eqn{\hat\tau} is the MAD of the residuals.

The estimated coefficients are inner products of the coefficient function with
the eigenfunctions. Because the dictionary
\eqn{g=(\phi^{(0)}_1,\ldots,\phi^{(0)}_{K_X},\phi^{(1)}_1,\ldots,
\phi^{(1)}_{K_X},\phi^U_1,\ldots,\phi^U_{K_U})} is not orthonormal, a
coefficient function is recovered from its inner products by solving the Gram
system,

\deqn{
\widehat\beta(s) = \sum_m \widehat d_m \widehat g_m(s), \qquad
\widehat d = \widehat G^{+}\widehat c, \qquad
\widehat G_{mn} = \int \widehat g_m(s)\widehat g_n(s)ds,
}

with trapezoidal quadrature on \code{tps}. The result is the orthogonal
projection of the coefficient onto the span of the dictionary and does not
depend on how the eigenvectors are scaled. \eqn{\widehat G^{+}} discards
eigenvalues below \code{gram_tol} times the largest one. This matters for the
single coefficient \code{beta_hat}: the subject-level and visit-level
eigenfunctions often span nearly the same functions, so \eqn{\widehat G} is
close to singular and the inner products along the nearly null directions carry
only noise. The same reconstruction is applied to the sub-dictionaries to give
the component coefficients.
}
\value{
A list with components
  \item{fit}{The final weighted \code{bam} fit.}
  \item{fitted}{Fitted values of the response.}
  \item{beta_hat}{Single coefficient function reconstructed from the full
dictionary (all of \eqn{a}, \eqn{b} and \eqn{c}).}
  \item{beta_B}{Subject-level coefficient reconstructed from
\eqn{(\phi^{(0)},\phi^{(1)})} and \eqn{(a,b)}; this is the coefficient of the
subject-level component \eqn{B_{i0}+t_{ij}B_{i1}} when both parts share one
coefficient.}
  \item{beta_0, beta_1}{Coefficients of the functional intercept and slope,
reconstructed from \eqn{\phi^{(0)}} with \eqn{a} and from \eqn{\phi^{(1)}} with
\eqn{b}.}
  \item{beta_U}{Visit-level coefficient reconstructed from \eqn{\phi^U} and
\eqn{c}.}
  \item{delta}{Estimated scalar coefficients.}
  \item{a, b, cU}{Estimated score coefficients.}
  \item{weights}{Final robust observation weights.}
  \item{tps}{Grid used for the reconstruction.}
  \item{design}{Data frame used in the score regression.}
}
\seealso{
\code{\link{rlfpca_estimate}},
\code{\link{simulate_data}}
}
\examples{
dat <- simulate_data(case = "2", nI = 100, nV = 4,
                     contamination = "mixed", prop = 0.05, seed = 1)

lfpca <- rlfpca_estimate(Y = dat$X, subject = dat$subject, Time = dat$times,
                         L = 0.95, smooth = TRUE, bf = 8, df_time_mean = 5,
                         mean_iter = 3, clean_cells = TRUE, clean_c = 4.5,
                         score_iter = 8, tps = dat$tps)

fit <- fit_rlfpcr(y = dat$y, lfpca = lfpca, subject = dat$subject,
                  times = dat$times, scalar_cov = dat$scalar_cov,
                  max_iter = 4, psi = "bisquare", c = 4.685)

w <- RLFPCR:::trap_weights(dat$tps)
c(MSE = mean((fit$fitted - dat$mu_true)^2),
  ISE_B = sum(w * (fit$beta_B - dat$truth$beta_B_true)^2),
  ISE_U = sum(w * (fit$beta_U - dat$truth$beta_U_true)^2))

plot(dat$tps, dat$truth$beta_B_true, type = "l", lwd = 2,
     xlab = "Functional domain", ylab = "Coefficient")
lines(dat$tps, fit$beta_B, col = "blue")
}
