\name{simulate_data}
\alias{simulate_data}
\title{
Simulate Longitudinal Functional Regression Data
}
\description{
Generates the simulation designs of the RLFPCR study: repeatedly observed
functional predictors, repeated scalar responses, scalar explanatory variables,
subject-level random intercepts and, optionally, one of seven contamination
mechanisms in the functional predictor, the scalar response or the scalar
covariates.
}
\usage{
simulate_data(case = c("1", "2"),
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
              contamination = c("none", "mixed", "spike", "drift",
                                "cluster", "covariate", "batch",
                                "instrument"),
              prop = 0.05,
              amp_subject = 4,
              amp_visit = 4,
              amp_y = 5,
              amp_spike = 3,
              spike_width = 2,
              amp_drift = 4,
              amp_cluster = 4,
              amp_cov = 5,
              amp_batch = 4,
              drift_rate = 2,
              seed = NULL)
}
\arguments{
  \item{case}{
Simulation case, \code{"1"} or \code{"2"}. Case \code{"1"} uses a single
coefficient function acting on the whole latent predictor; Case \code{"2"} uses
separate coefficient functions for the subject-level and visit-level components.
}
  \item{nI}{
Number of subjects. The default is \code{100}.
}
  \item{nV}{
Number of visits per subject (balanced design). The default is \code{4}.
}
  \item{tps}{
Grid on which the curves are observed. The basis functions are defined on
\eqn{[0, 100]}. The default is \code{seq(0.5, 99.5, by = 1)}.
}
  \item{sigma_x}{
Standard deviation of the measurement error of the functional predictor.
}
  \item{sigma_y}{
Standard deviation of the response error. The default is \code{0.80}.
}
  \item{sigma_b}{
Standard deviation of the subject random intercept. The default is \code{0.60}.
}
  \item{lambda}{
Variances of the subject-level scores; its length is the number of subject-level
components (at most four in the default bases).
}
  \item{nu}{
Variances of the visit-level scores; its length (at most four) is the number of
visit-level components.
}
  \item{scalar_cov}{
\code{NULL}, an integer between 1 and 4 (number of scalar covariates generated
internally), or a matrix of scalar covariates. The default is \code{4}.
}
  \item{beta_scalar}{
True scalar coefficients. The default is \code{c(0.6, 0.8, -0.5, 0.2)}.
}
  \item{contamination}{
Contamination mechanism, see Details. The default \code{"none"} gives clean data.
}
  \item{prop}{
Contamination proportion. For \code{"mixed"} it is the proportion of
contaminated subjects, of contaminated visits and of contaminated responses; for
\code{"spike"}, \code{"cluster"}, \code{"covariate"} and \code{"batch"} the
proportion of affected visits; for \code{"drift"} the proportion of affected
subjects. It is not used by \code{"instrument"}.
}
  \item{amp_subject, amp_visit, amp_y}{
Magnitudes of the subject-level, visit-level and response perturbations of the
\code{"mixed"} mechanism, in units of the standard deviation of the clean
predictor values and of the clean responses.
}
  \item{amp_spike, spike_width}{
Peak height (in predictor standard deviations) and Gaussian width (in grid
units) of the local spikes.
}
  \item{amp_drift, amp_cluster, amp_batch}{
Magnitudes of the subject drift, the clustered visit outliers and the batch
shift, in predictor standard deviations.
}
  \item{amp_cov}{
Size of the shift added to the last scalar covariate, in standard deviations of
that covariate.
}
  \item{drift_rate}{
Shift per unit of follow-up time in the \code{"instrument"} mechanism, in
predictor standard deviations.
}
  \item{seed}{
Optional seed. The clean data are generated from \code{seed} and the
contamination from \code{seed + 100000}, so that all contamination mechanisms
applied with the same seed share the same clean data set.
}
}
\details{
Visit times are \eqn{T_{i1}=0} and \eqn{T_{ij}=\sum_{l<j}U_{il}} with
\eqn{U_{il}} uniform on \eqn{(0.3, 1)}; \eqn{t_{ij}} is the standardized time.
The latent predictor is

\deqn{
W_{ij}(s) = \mu(s,T_{ij}) + c_B\sum_{k}\xi_{ik}\{\phi_k^{(0)}(s)+t_{ij}\phi_k^{(1)}(s)\}
+ c_U\sum_{\ell}\zeta_{ij\ell}\phi_\ell^U(s),
}

with \eqn{\mu(s,T)=c_\mu(T/5-s/100)^2}, \eqn{\phi_k^{(0)}} and \eqn{\phi_k^{(1)}}
Fourier functions, \eqn{\phi_\ell^U} normalized Legendre polynomials of degree
\eqn{\ell-1}, \eqn{(c_\mu,c_B,c_U)=(0.05,0.90,1.00)} in Case 1 and
\eqn{(0.03,1.10,0.70)} in Case 2. The observed curves are
\eqn{X_{ij}(s_d)=W_{ij}(s_d)+\varepsilon_{ij}(s_d)}.

In Case 1 the response is
\eqn{Y_{ij}=b_i+q_{ij}^{\top}\delta+\int W_{ij}(s)\beta^\ast(s)ds+e_{ij}} with
\eqn{\beta^\ast=1.40\beta} and \eqn{\beta} a normalized combination of
\eqn{\phi_1^{(0)}}, \eqn{\phi_2^{(0)}} and \eqn{\phi_1^U}. In Case 2,
\eqn{Y_{ij}=b_i+q_{ij}^{\top}\delta+\int B_{ij}(s)\beta_B^\ast(s)ds+
\int U_{ij}(s)\beta_U^\ast(s)ds+e_{ij}} with \eqn{\beta_B^\ast=1.35\beta_B} and
\eqn{\beta_U^\ast=1.10\beta_U}. The returned \code{truth} contains these
effective coefficients \eqn{\beta^\ast}, \eqn{\beta_B^\ast} and
\eqn{\beta_U^\ast}, which are the targets of estimation, together with the
normalized shapes.

Contamination is added after the clean data are generated and never changes
\code{mu_true}:
\code{"mixed"} adds a random-sign subject-level template to all curves of a
proportion of subjects, a random-sign visit-level template to a proportion of
the remaining curves, and heavy-tailed shifts to a proportion of the responses;
\code{"spike"} adds a narrow Gaussian bump with random location and sign;
\code{"drift"} adds to the curves of selected subjects a constant offset that
grows linearly from zero at baseline to its full size at the last visit;
\code{"cluster"} adds the same visit-level template with the same sign to the
selected visits;
\code{"covariate"} shifts the last scalar covariate of the selected visits;
\code{"batch"} adds the same constant offset (a global intensity shift) to the
selected visits, as a second scanner would;
\code{"instrument"} adds to every curve a constant offset proportional to the
time since baseline. The offset of these three mechanisms is the normalized
constant function, which is not orthogonal to the coefficient functions, so the
shifts change the conditional mean implied by the observed curves.
}
\value{
A list with components \code{case}, \code{contamination}, \code{nI}, \code{n},
\code{Jvec}, \code{subject}, \code{times} (time since baseline), \code{tstd}
(standardized time), \code{tps}, \code{X}, \code{y}, \code{X_clean},
\code{y_clean}, \code{scalar_cov} (observed covariates), \code{scalar_cov_clean},
\code{beta_scalar_true}, \code{scalar_names}, \code{mu_true} (true conditional
mean of the clean response), \code{W_latent}, \code{B_comp}, \code{U_comp},
\code{eta_mat}, \code{truth}, and the index vectors \code{idx_subject_out},
\code{idx_visit_out}, \code{idx_y_out} and \code{idx_cov_out} of the
contaminated observations.
}
\seealso{
\code{\link{rlfpca_estimate}},
\code{\link{fit_rlfpcr}}
}
\examples{
dat <- simulate_data(case = "2", nI = 100, nV = 4,
                     contamination = "mixed", prop = 0.05, seed = 1)
str(dat$truth)
length(dat$idx_visit_out)

plot(dat$tps, dat$truth$beta_B_true, type = "l",
     xlab = "Functional domain", ylab = "Effective coefficient")
lines(dat$tps, dat$truth$beta_U_true, lty = 2)

batch <- simulate_data(case = "1", contamination = "batch", prop = 0.2, seed = 1)
all.equal(batch$mu_true, simulate_data(case = "1", seed = 1)$mu_true)
}
