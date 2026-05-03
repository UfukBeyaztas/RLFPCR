\name{simulate_data}
\alias{simulate_data}
\title{
Simulate Longitudinal Functional Regression Data
}
\description{
Generates simulated data for robust longitudinal scalar-on-function regression.
The function creates repeatedly observed functional predictors, repeated scalar
responses, optional scalar explanatory variables, subject-level random effects,
and optional structured outliers in the functional predictor and scalar response.
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
              scalar_cov = NULL,
              beta_scalar = c(0.6, 0.8, -0.5, 0.2),
              outlier = FALSE,
              prop_curve = 0.05,
              prop_y = 0.05,
              prop_subject = 0.02,
              amp_subject = 4,
              amp_visit = 4,
              amp_y = 5,
              seed = NULL)
}
\arguments{
  \item{case}{
A character string specifying the simulation case. Available choices are
\code{"1"} and \code{"2"}. Case \code{"1"} uses a single global functional
coefficient, whereas Case \code{"2"} uses separate subject-level and visit-level
functional coefficient components.
}

  \item{nI}{
Number of subjects. The default is \code{100}.
}

  \item{nV}{
Number of visits per subject. The current implementation uses a balanced
longitudinal design with \code{nV} repeated observations for each subject. The
default is \code{4}.
}

  \item{tps}{
Numeric vector giving the grid points on which the functional predictor is
observed. The default is \code{seq(0.5, 99.5, by = 1)}.
}

  \item{sigma_x}{
Standard deviation of the measurement error added to the latent functional
predictor. The default is \code{0.01}.
}

  \item{sigma_y}{
Standard deviation of the scalar response error. The default is \code{0.80}.
}

  \item{sigma_b}{
Standard deviation of the subject-specific random intercept in the scalar
response model. The default is \code{0.60}.
}

  \item{lambda}{
Numeric vector of eigenvalues controlling the variance of the subject-level
longitudinal functional scores. The length of \code{lambda} determines the
number of subject-level components generated. The default is
\code{c(2.0, 1.2, 0.7, 0.4)}.
}

  \item{nu}{
Numeric vector of eigenvalues controlling the variance of the visit-level
functional scores. The length of \code{nu} determines the number of visit-level
components generated. The default is \code{c(1.2, 0.7, 0.4)}.
}

  \item{scalar_cov}{
Optional scalar explanatory variables. If \code{NULL}, no scalar explanatory
variables are included. If a single integer between \code{1} and \code{4} is
given, that number of scalar explanatory variables is generated internally. If a
matrix or data frame is supplied, it is used as the scalar covariate matrix.
}

  \item{beta_scalar}{
Numeric vector of true scalar coefficients used when scalar explanatory
variables are included. The default is \code{c(0.6, 0.8, -0.5, 0.2)}.
}

  \item{outlier}{
Logical; if \code{TRUE}, structured contamination is added to the data. The
default is \code{FALSE}.
}

  \item{prop_curve}{
Proportion of subject-visit functional observations contaminated by visit-level
functional outliers when \code{outlier = TRUE}. The default is \code{0.05}.
}

  \item{prop_y}{
Proportion of scalar responses contaminated by response outliers when
\code{outlier = TRUE}. The default is \code{0.05}.
}

  \item{prop_subject}{
Proportion of subjects contaminated by subject-level functional outliers when
\code{outlier = TRUE}. The default is \code{0.02}.
}

  \item{amp_subject}{
Magnitude of the subject-level functional outlier perturbation. Larger values
produce more severe subject-level functional contamination. The default is
\code{4}.
}

  \item{amp_visit}{
Magnitude of the visit-level functional outlier perturbation. Larger values
produce more severe visit-level functional contamination. The default is
\code{4}.
}

  \item{amp_y}{
Magnitude of the scalar response outlier perturbation. Larger values produce
more severe response contamination. The default is \code{5}.
}

  \item{seed}{
Optional random seed for reproducibility. The default is \code{NULL}.
}
}
\details{
The function generates data from a longitudinal scalar-on-function regression
framework in which both the scalar response and the functional predictor are
observed repeatedly over visits. For subject \eqn{i} and visit \eqn{j}, the
general regression structure is

\deqn{
Y_{ij}
=
\alpha(T_{ij})
+
q_{ij}^{\top}\delta
+
b_i
+
\int_{\mathcal S}
\{X_{ij}(s)-\mu(s,T_{ij})\}\beta(s)\,ds
+
e_{ij},
}

where \eqn{Y_{ij}} is the scalar response, \eqn{T_{ij}} is the visit time,
\eqn{q_{ij}} is an optional vector of scalar explanatory variables,
\eqn{\delta} is the corresponding scalar coefficient vector, \eqn{b_i} is a
subject-specific random intercept, \eqn{X_{ij}(s)} is the functional predictor,
\eqn{\mu(s,T_{ij})} is the longitudinal mean surface, and \eqn{\beta(s)} is the
functional coefficient.

The latent functional predictor is generated according to the longitudinal
functional decomposition

\deqn{
W_{ij}(s)
=
\mu(s,T_{ij})
+
\sum_{k=1}^{K_X}
\xi_{ik}
\{\phi_k^{(0)}(s)+t_{ij}\phi_k^{(1)}(s)\}
+
\sum_{\ell=1}^{K_U}
\zeta_{ij\ell}\phi_\ell^U(s),
}

where \eqn{\xi_{ik}} are subject-level scores, \eqn{\zeta_{ij\ell}} are
visit-level scores, \eqn{\phi_k^{(0)}} and \eqn{\phi_k^{(1)}} represent
subject-level intercept and slope basis functions, and \eqn{\phi_\ell^U}
represents visit-level basis functions. The observed functional predictor is

\deqn{
X_{ij}(s_d)
=
W_{ij}(s_d)
+
\varepsilon_{ij}^{X}(s_d),
}

where \eqn{\varepsilon_{ij}^{X}(s_d)} is independent Gaussian measurement error
with standard deviation \code{sigma_x}.

In Case \code{"1"}, the scalar response depends on the full latent functional
signal through a single global coefficient function,

\deqn{
Y_{ij}
=
b_i
+
q_{ij}^{\top}\delta
+
1.40
\int_{\mathcal S}
W_{ij}(s)\beta(s)\,ds
+
e_{ij}.
}

In Case \code{"2"}, the response depends separately on the subject-level and
visit-level functional components,

\deqn{
Y_{ij}
=
b_i
+
q_{ij}^{\top}\delta
+
1.35
\int_{\mathcal S}
B_{ij}(s)\beta_B(s)\,ds
+
1.10
\int_{\mathcal S}
U_{ij}(s)\beta_U(s)\,ds
+
e_{ij},
}

where \eqn{B_{ij}(s)} denotes the subject-level longitudinal component and
\eqn{U_{ij}(s)} denotes the visit-level deviation component.

If \code{outlier = TRUE}, three types of structured contamination may be
introduced. First, a subset of subjects receives a subject-level functional
perturbation, affecting all visits for those subjects. Second, a subset of
individual subject-visit curves receives a visit-level functional perturbation.
Third, a subset of scalar responses receives heavy-tailed additive
contamination. The clean versions of the functional predictor and scalar
response are also returned, which allows direct comparison between clean and
contaminated data.

The output of this function is designed to be used directly with
\code{\link{rlfpca_estimate}} and \code{\link{fit_rlfpcr}}.
}
\value{
A list with the following components:

  \item{case}{
The simulation case used, either \code{"1"} or \code{"2"}.
}

  \item{nI}{
Number of subjects.
}

  \item{n}{
Total number of subject-visit observations, equal to \code{nI * nV}.
}

  \item{Jvec}{
A vector containing the number of visits for each subject.
}

  \item{subject}{
Vector of subject identifiers of length \code{n}.
}

  \item{times}{
Vector of standardized visit times of length \code{n}.
}

  \item{tps}{
Functional grid points on which the curves are observed.
}

  \item{X}{
Observed functional predictor matrix. Rows correspond to subject-visit
observations and columns correspond to functional grid points. If
\code{outlier = TRUE}, this matrix contains the contaminated functional
observations.
}

  \item{y}{
Observed scalar response vector. If \code{outlier = TRUE}, this vector contains
the contaminated scalar responses.
}

  \item{X_clean}{
Clean functional predictor matrix before structured outliers are added.
}

  \item{y_clean}{
Clean scalar response vector before response outliers are added.
}

  \item{scalar_cov}{
Matrix of scalar explanatory variables, or \code{NULL} if no scalar explanatory
variables are used.
}

  \item{beta_scalar_true}{
True scalar coefficient vector used to generate the scalar covariate effect.
}

  \item{scalar_names}{
Names of the scalar explanatory variables.
}

  \item{mu_true}{
True conditional mean of the scalar response before adding response noise.
}

  \item{W_latent}{
Latent functional signal before adding functional measurement error.
}

  \item{B_comp}{
Subject-level longitudinal functional component.
}

  \item{U_comp}{
Visit-level functional deviation component.
}

  \item{eta_mat}{
Mean surface component evaluated at the observed visit times and functional grid
points.
}

  \item{truth}{
A list containing the true coefficient functions. In Case \code{"1"},
\code{truth$beta_true} contains the global coefficient function, while
\code{truth$beta_B_true} and \code{truth$beta_U_true} are \code{NULL}. In Case
\code{"2"}, \code{truth$beta_B_true} and \code{truth$beta_U_true} contain the
subject-level and visit-level coefficient functions, while
\code{truth$beta_true} is \code{NULL}.
}

  \item{outlier}{
Logical indicator showing whether structured contamination was added.
}

  \item{idx_subject_out}{
Indices of subject-visit observations belonging to contaminated subjects. This
is an empty integer vector when \code{outlier = FALSE}.
}

  \item{idx_visit_out}{
Indices of subject-visit observations contaminated by visit-level functional
outliers. This is an empty integer vector when \code{outlier = FALSE}.
}

  \item{idx_y_out}{
Indices of scalar responses contaminated by response outliers. This is an empty
integer vector when \code{outlier = FALSE}.
}
}
\note{
The current implementation generates balanced longitudinal data with the same
number of visits for each subject. The functional predictor is observed densely
on a common grid \code{tps}. The function is intended primarily for simulation
studies, reproducible examples, and validation of the RLFPCR estimation
pipeline.
}

\seealso{
\code{\link{rlfpca_estimate}},
\code{\link{fit_rlfpcr}}
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

str(dat)

matplot(dat$tps, t(dat$X[1:5, ]), type = "l",
        xlab = "Functional domain",
        ylab = "Functional predictor")

plot(dat$tps, dat$truth$beta_true, type = "l",
     xlab = "Functional domain",
     ylab = "True coefficient function")

dat_out <- simulate_data(
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

length(dat_out$idx_subject_out)
length(dat_out$idx_visit_out)
length(dat_out$idx_y_out)
}
