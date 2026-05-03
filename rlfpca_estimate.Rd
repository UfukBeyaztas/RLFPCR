\name{rlfpca_estimate}
\alias{rlfpca_estimate}
\title{
Robust Longitudinal Functional Principal Component Analysis
}
\description{
Estimates a robust longitudinal functional principal component decomposition
for repeatedly observed functional predictors. The function robustly estimates
the mean surface, subject-level and visit-level covariance components,
eigenfunctions, eigenvalues, and longitudinal FPCA scores.
}
\usage{
rlfpca_estimate(Y, subject, Time, L = 0.90,
                N.X = NA, N.U = NA,
                smooth = FALSE, bf = 8,
                df_time_mean = 5,
                mean_iter = 3,
                clean_cells = TRUE,
                clean_c = 4.5,
                score_iter = 8)
}
\arguments{
  \item{Y}{
A numeric matrix containing the longitudinal functional predictor. Rows
correspond to subject-visit observations and columns correspond to the
functional grid points.
}

  \item{subject}{
A vector of subject identifiers. The length of \code{subject} must be equal to
\code{nrow(Y)}. Repeated observations from the same subject must have the same
identifier.
}

  \item{Time}{
A numeric vector of visit times. The length of \code{Time} must be equal to
\code{nrow(Y)}. The time variable is internally standardized before estimating
the longitudinal functional decomposition.
}

  \item{L}{
Target proportion of total variation explained by the retained subject-level
and visit-level components when \code{N.X} or \code{N.U} is not supplied. The
default is \code{0.90}.
}

  \item{N.X}{
Number of subject-level longitudinal functional principal components to retain.
If \code{NA}, the number is selected jointly with \code{N.U} using the pooled
explained-variance rule determined by \code{L}. The default is \code{NA}.
}

  \item{N.U}{
Number of visit-level functional principal components to retain. If \code{NA},
the number is selected jointly with \code{N.X} using the pooled
explained-variance rule determined by \code{L}. The default is \code{NA}.
}

  \item{smooth}{
Logical; if \code{TRUE}, the estimated covariance surfaces are smoothed using
bivariate penalized spline smoothing. If \code{FALSE}, the raw robust covariance
estimates are used, with the diagonal of the visit-level covariance estimated
separately. The default is \code{FALSE}.
}

  \item{bf}{
Basis dimension used in the bivariate smoothing step for covariance surfaces.
This argument is passed to the internal smoothing routines. The default is
\code{8}.
}

  \item{df_time_mean}{
Degrees of freedom used for the spline basis in the robust mean surface
estimation step. The default is \code{5}.
}

  \item{mean_iter}{
Maximum number of robust reweighting iterations used to estimate the mean
surface. The default is \code{3}.
}

  \item{clean_cells}{
Logical; if \code{TRUE}, unusually large cellwise residual values in the
centered functional data matrix are winsorized before covariance estimation.
The default is \code{TRUE}.
}

  \item{clean_c}{
Tuning constant used in the cellwise residual cleaning step. The default is
\code{4.5}.
}

  \item{score_iter}{
Maximum number of alternating blockwise score-update iterations used to estimate
the subject-level and visit-level scores. The default is \code{8}.
}
}
\details{
The function implements the robust longitudinal functional principal component
analysis step used in robust longitudinal functional principal component
regression. Let \eqn{X_{ij}(s)} denote the functional predictor for subject
\eqn{i} at visit \eqn{j}, and let \eqn{T_{ij}} denote the corresponding visit
time. The longitudinal functional predictor is represented as

\deqn{
X_{ij}(s)-\mu(s,T_{ij})
=
B_{i0}(s)
+
t_{ij}B_{i1}(s)
+
U_{ij}(s)
+
\varepsilon_{ij}^{X}(s),
}

where \eqn{\mu(s,T_{ij})} is the mean surface, \eqn{B_{i0}(s)} is the
subject-specific functional intercept process, \eqn{B_{i1}(s)} is the
subject-specific functional slope process, \eqn{U_{ij}(s)} is the visit-specific
functional deviation process, and \eqn{\varepsilon_{ij}^{X}(s)} is a residual
measurement-error process.

The function first estimates the mean surface robustly by iteratively weighted
least squares. At each iteration, curvewise residual norms are computed and
robust curve weights are updated. These weights reduce the influence of
atypical functional observations when estimating the mean surface and the
covariance components.

After centering the curves, the function estimates the subject-level and
visit-level covariance components from the cross-product regression identity

\deqn{
E\{ \widetilde X_{ij}(s)\widetilde X_{ik}(s') \mid T_{ij},T_{ik} \}
=
G_0(s,s')
+
t_{ik}G_{10}(s,s')
+
t_{ij}G_{01}(s,s')
+
t_{ij}t_{ik}G_1(s,s')
+
I(j=k)G_U(s,s').
}

The kernels \eqn{G_0}, \eqn{G_1}, and \eqn{G_{01}} describe the subject-level
functional intercept/slope covariance structure, whereas \eqn{G_U} describes
the visit-level covariance structure. Robustness is introduced through
pairwise weights formed from the curve weights.

The estimated subject-level covariance operator is constructed as the block
operator

\deqn{
\widehat K_X
=
\begin{pmatrix}
\widehat K_0 & \widehat K_{01} \\
\widehat K_{10} & \widehat K_1
\end{pmatrix},
}

and the visit-level covariance operator is denoted by \eqn{\widehat K_U}.
Eigenanalysis of these two covariance operators yields subject-level
eigenfunctions

\deqn{
\left\{
\left(
\widehat\phi_k^{(0)},\widehat\phi_k^{(1)}
\right)
:
k=1,\ldots,K_X
\right\}
}

and visit-level eigenfunctions

\deqn{
\left\{
\widehat\phi_\ell^U:
\ell=1,\ldots,K_U
\right\}.
}

If \code{N.X} and \code{N.U} are not supplied, the function selects them
sequentially by pooling the positive subject-level and visit-level eigenvalues
and retaining components until the target explained variation \code{L} is
reached.

Finally, the robust longitudinal scores are estimated from the truncated
decomposition

\deqn{
\widetilde X_{ij}(s)
\approx
\sum_{k=1}^{K_X}
\xi_{ik}
\{
\widehat\phi_k^{(0)}(s)
+
t_{ij}\widehat\phi_k^{(1)}(s)
\}
+
\sum_{\ell=1}^{K_U}
\zeta_{ij\ell}\widehat\phi_\ell^U(s).
}

The score estimation is performed using alternating blockwise updates for the
subject-level scores \eqn{\xi_{ik}} and visit-level scores \eqn{\zeta_{ij\ell}}.
This avoids forming a single large global score system and exploits the
longitudinal block structure of the data.

The returned object is designed to be passed directly to
\code{\link{fit_rlfpcr}}, which uses the estimated scores and eigenfunctions in
the robust scalar-on-function regression stage.
}
\value{
A list with the following components:

  \item{phi.0}{
A matrix containing the estimated subject-level functional intercept
eigenfunctions \eqn{\widehat\phi_k^{(0)}}. Rows correspond to functional grid
points and columns correspond to retained subject-level components.
}

  \item{phi.1}{
A matrix containing the estimated subject-level functional slope eigenfunctions
\eqn{\widehat\phi_k^{(1)}}. Rows correspond to functional grid points and
columns correspond to retained subject-level components.
}

  \item{phi.U}{
A matrix containing the estimated visit-level eigenfunctions
\eqn{\widehat\phi_\ell^U}. Rows correspond to functional grid points and
columns correspond to retained visit-level components.
}

  \item{xi}{
A matrix of estimated subject-level scores. Rows correspond to subjects and
columns correspond to retained subject-level components.
}

  \item{zeta}{
A matrix of estimated visit-level scores. Rows correspond to subject-visit
observations and columns correspond to retained visit-level components.
}

  \item{N.X}{
The number of retained subject-level longitudinal functional principal
components.
}

  \item{N.U}{
The number of retained visit-level functional principal components.
}

  \item{eta}{
A matrix containing the estimated robust mean surface evaluated at the observed
visit times and functional grid points. It has the same dimension as \code{Y}.
}

  \item{Y.tilde}{
The robustly centered functional data matrix, obtained after subtracting
\code{eta} from \code{Y}. If \code{clean_cells = TRUE}, this component also
reflects the cellwise residual cleaning step.
}

  \item{sigma2}{
Estimated residual measurement-error variance of the functional predictor.
}

  \item{lambda}{
Estimated eigenvalues of the retained subject-level covariance operator.
}

  \item{nu}{
Estimated eigenvalues of the retained visit-level covariance operator.
}

  \item{curve_w}{
Final robust curve weights obtained from the mean-surface estimation step.
}
}
\note{
The theoretical results in the accompanying methodology are developed for fixed
truncation dimensions \eqn{K_X} and \eqn{K_U}. In the implementation,
\code{N.X} and \code{N.U} may also be selected adaptively using the pooled
explained-variance rule controlled by \code{L}.

The function assumes that the functional predictors are densely observed on a
common grid. Sparse or irregularly sampled functional data should be
preprocessed before using this function.
}
\seealso{
\code{\link{fit_rlfpcr}},
\code{\link{simulate_data}}
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

str(lfpca)

matplot(dat$tps, lfpca$phi.0, type = "l",
        xlab = "Functional domain",
        ylab = "Estimated subject-level eigenfunctions")
}
