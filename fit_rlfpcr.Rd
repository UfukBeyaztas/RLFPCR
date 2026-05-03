\name{fit_rlfpcr}
\alias{fit_rlfpcr}
\title{
Fit Robust Longitudinal Functional Principal Component Regression
}
\description{
Fits the second-stage robust longitudinal functional principal component
regression model. The function uses robust longitudinal FPCA scores obtained
from \code{\link{rlfpca_estimate}} and fits a robust additive mixed score
regression for a repeated scalar response.
}
\usage{
fit_rlfpcr(y, lfpca, subject, times,
           scalar_cov = NULL,
           max_iter = 4,
           psi = "bisquare",
           c = 4.685)
}
\arguments{
  \item{y}{
A numeric vector containing the repeated scalar response values. Its length must
be equal to the total number of subject-visit observations.
}

  \item{lfpca}{
An object returned by \code{\link{rlfpca_estimate}}. It must contain the robust
longitudinal FPCA scores and eigenfunctions, including at least
\code{xi}, \code{zeta}, \code{phi.0}, \code{phi.1}, and \code{phi.U}.
}

  \item{subject}{
A vector of subject identifiers. The length of \code{subject} must be the same
as the length of \code{y}. Repeated observations belonging to the same subject
must share the same identifier.
}

  \item{times}{
A numeric vector of visit times. The length of \code{times} must be the same
as the length of \code{y}. The time variable is used both in the smooth
population time effect and in the interaction between subject-level scores and
visit time.
}

  \item{scalar_cov}{
An optional matrix or data frame of scalar explanatory variables. The number of
rows must be equal to the length of \code{y}. If column names are not supplied,
default names \code{q1}, \code{q2}, ... are assigned internally. The default is
\code{NULL}.
}

  \item{max_iter}{
Maximum number of iteratively reweighted least-squares iterations used in the
robust score-space regression. The default is \code{4}.
}

  \item{psi}{
Type of robust score/weight function used for the scalar-response residuals.
Currently supported choices are \code{"bisquare"} and \code{"huber"}, as
implemented by the internal weight function. The default is \code{"bisquare"}.
}

  \item{c}{
Tuning constant for the robust weight function. The default value is
\code{4.685}, corresponding to the usual Tukey bisquare tuning constant.
}
}
\details{
The function implements the final regression stage of robust longitudinal
functional principal component regression (RLFPCR). Let \eqn{Y_{ij}} denote the
scalar response for subject \eqn{i} at visit \eqn{j}, let \eqn{T_{ij}} denote
the visit time, and let \eqn{X_{ij}(s)} denote the longitudinal functional
predictor. The underlying scalar-on-function model is

\deqn{
Y_{ij}
=
\alpha(T_{ij})
+
q_{ij}^{\top}\delta
+
r_{ij}^{\top}u_i
+
\int_{\mathcal S}
\{X_{ij}(s)-\mu(s,T_{ij})\}\beta(s) ds
+
e_{ij}.
}

The robust longitudinal FPCA step decomposes the centered functional predictor
as

\deqn{
X_{ij}(s)-\mu(s,T_{ij})
\approx
\sum_{k=1}^{K_X}
\xi_{ik}\{\phi_k^{(0)}(s)+t_{ij}\phi_k^{(1)}(s)\}
+
\sum_{\ell=1}^{K_U}
\zeta_{ij\ell}\phi_{\ell}^{U}(s),
}

where \eqn{\xi_{ik}} are subject-level scores, \eqn{\zeta_{ij\ell}} are
visit-level scores, \eqn{\phi_k^{(0)}} and \eqn{\phi_k^{(1)}} describe the
subject-level baseline and temporal components, and \eqn{\phi_{\ell}^{U}}
describes the visit-level deviation component.

Substitution of this truncated representation into the scalar-on-function model
gives the finite-dimensional score regression

\deqn{
Y_{ij}
=
\alpha(T_{ij})
+
q_{ij}^{\top}\delta
+
r_{ij}^{\top}u_i
+
\sum_{k=1}^{K_X} a_k \xi_{ik}
+
\sum_{k=1}^{K_X} b_k t_{ij}\xi_{ik}
+
\sum_{\ell=1}^{K_U} c_{\ell}\zeta_{ij\ell}
+
e_{ij}.
}

The function fits this model using \code{gam}. The smooth population time
effect is represented by \code{s(times, bs = "ps", m = 1)}, and subject-specific
dependence in the scalar response is represented by the random-effect smooth
\code{s(subject, bs = "re")}. Robustness is introduced by iteratively updating
observation weights according to

\deqn{
w_{ij}
=
\frac{\psi(r_{ij}/\hat \tau)}{r_{ij}/\hat \tau},
}

where \eqn{r_{ij}} is the scalar-response residual and \eqn{\hat\tau} is a
MAD-type robust scale estimate.

After estimating the score coefficients, the functional coefficient is
reconstructed as

\deqn{
\widehat\beta(s)
=
\sum_{k=1}^{K_X}\widehat a_k \widehat\phi_k^{(0)}(s)
+
\sum_{k=1}^{K_X}\widehat b_k \widehat\phi_k^{(1)}(s)
+
\sum_{\ell=1}^{K_U}\widehat c_{\ell}\widehat\phi_{\ell}^{U}(s).
}

Thus, the returned coefficient function decomposes the functional effect into
a subject-level baseline component, a subject-level longitudinal/time component,
and a visit-level deviation component.
}
\value{
A list with the following components:

  \item{fit}{
The final weighted \code{gam} object fitted after the robust IRLS
iterations.
}

  \item{beta_hat}{
A numeric vector containing the reconstructed functional coefficient
\eqn{\widehat\beta(s)} evaluated on the same functional grid as the eigenfunctions
stored in \code{lfpca}.
}

  \item{delta}{
A numeric vector of estimated coefficients for the scalar explanatory variables.
If \code{scalar_cov = NULL}, this is a numeric vector of length zero.
}

  \item{a}{
Estimated coefficients corresponding to the subject-level baseline scores
\eqn{\xi_{ik}}.
}

  \item{b}{
Estimated coefficients corresponding to the subject-level longitudinal score
interactions \eqn{t_{ij}\xi_{ik}}.
}

  \item{cU}{
Estimated coefficients corresponding to the visit-level deviation scores
\eqn{\zeta_{ij\ell}}.
}

  \item{weights}{
Final robust observation weights from the IRLS procedure.
}

  \item{design}{
The data frame used internally to fit the score-space regression model. It
contains the scalar response, visit times, subject indicator, optional scalar
explanatory variables, and the constructed longitudinal FPCA score covariates.
}
}
\note{
This function should usually be applied after \code{\link{rlfpca_estimate}}.
The theoretical development assumes fixed truncation dimensions, whereas the
default implementation of \code{\link{rlfpca_estimate}} may select the number of
components by an explained-variance criterion.
}
\examples{

dat <- simulate_data(
  case = "1",
  nI = 100,
  nV = 4,
  tps = seq(0.5, 99.5, by = 1),
  sigma_x = 0.01,
  sigma_y = sqrt(2),
  sigma_b = sqrt(2),
  scalar_cov = 4,
  beta_scalar = c(0.6, 0.8, -0.5, 0.2),
  outlier = TRUE
)

X <- dat$X
y <- dat$y
subject <- dat$subject
times <- dat$times
tps <- dat$tps
scalar_cov <- dat$scalar_cov

lfpca <- rlfpca_estimate(
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

fit <- fit_rlfpcr(
  y = y,
  lfpca = lfpca,
  subject = subject,
  times = times,
  scalar_cov = scalar_cov,
  max_iter = 4,
  psi = "bisquare",
  c = 4.685
)

str(fit)
plot(dat$tps, fit$beta_hat, type = "l",
     xlab = "Functional domain",
     ylab = expression(hat(beta)(s)))
}
