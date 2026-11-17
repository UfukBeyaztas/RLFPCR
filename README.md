# RLFPCR

RLFPCR is an R package for robust longitudinal functional principal component
regression. It is meant for studies in which a scalar outcome and a densely
observed curve, such as a profile of imaging features, are recorded at
repeated visits of each subject, and some of the curves or outcomes may be
atypical. The method builds on the longitudinal functional principal component
regression of Greven et al. (2010) and Gertheiss et al. (2013) and makes each
of its steps resistant to atypical curves and outlying responses.

This is version 1.1.0, the version used in the paper cited below. The manual
is in `RLFPCR_1.1.0.pdf`.

## Installation

```r
install.packages("remotes")
remotes::install_github("UfukBeyaztas/RLFPCR")
```

RLFPCR needs R 3.5.0 or later and imports mgcv, splines and stats.

## Model and functions

For subject i at visit j, with visit time T_ij and standardized time t_ij, the
curve is decomposed as

    X_ij(s) = mu(s, T_ij) + B_i0(s) + t_ij B_i1(s) + U_ij(s) + error,

where B_i0 and B_i1 are the subject-level functional intercept and slope and
U_ij is the visit-level deviation. With xi_ik and zeta_ijl the subject-level
and visit-level scores, the regression of the outcome on the curve becomes the
score regression

    Y_ij = alpha(T_ij) + q_ij' delta + u_i + sum_k a_k xi_ik
           + sum_k b_k t_ij xi_ik + sum_l c_l zeta_ijl + e_ij,

with a smooth time effect alpha, scalar covariates q_ij and a random subject
intercept u_i. The coefficient functions are recovered from a_k, b_k and c_l.

| Function | What it does |
|---|---|
| `rlfpca_estimate()` | Robust longitudinal FPCA. The mean surface is estimated with Tukey bisquare weights for whole curves, outlying cells of the centred curves are winsorized, the covariance functions are estimated with the same curve weights, and the eigenfunctions and scores are computed. |
| `fit_rlfpcr()` | Fits the score regression by iteratively reweighted `mgcv::bam()` and reconstructs the coefficient functions through the Gram matrix of the estimated eigenfunctions. |
| `simulate_data()` | Generates the simulation designs of the paper, clean or with one of seven contamination mechanisms. |

`fit_rlfpcr()` returns a single coefficient function `beta_hat`, the
subject-level coefficient `beta_B`, its intercept and slope parts `beta_0` and
`beta_1`, the visit-level coefficient `beta_U`, the scalar effects `delta` and
the robust weights of the responses.

## Example

```r
library(RLFPCR)

# Case 1 of the paper: 100 subjects with 4 visits, 10% mixed contamination
dat <- simulate_data(case = "1", nI = 100, nV = 4,
                     contamination = "mixed", prop = 0.10, seed = 1)

lfpca <- rlfpca_estimate(Y = dat$X, subject = dat$subject, Time = dat$times,
                         L = 0.95, smooth = TRUE, tps = dat$tps)
fit <- fit_rlfpcr(y = dat$y, lfpca = lfpca, subject = dat$subject,
                  times = dat$times, scalar_cov = dat$scalar_cov)
```

The arguments that are not given keep their defaults, which are the settings
used in the paper: covariance basis dimension 8, Tukey constant 4.685, cells
winsorized at 4.5 MADs, at most four reweighting steps and Gram threshold
0.05. The two calls take about a second. `rlfpca_estimate()` may warn that
`ns()` moved interior knots inside the boundary. All first visits are at time
0, so a quantile knot falls on the boundary, and the warning can be ignored.

The simulated data come with the truth, so the fit can be checked:

```r
w <- RLFPCR:::trap_weights(dat$tps)   # trapezoidal weights on the grid
round(c(MSE = mean((fit$fitted - dat$mu_true)^2),
        MSE_delta = mean((fit$delta - dat$beta_scalar_true)^2),
        ISE = sum(w * (fit$beta_hat - dat$truth$beta_true)^2)), 3)
#>       MSE MSE_delta       ISE
#>     0.209     0.027     0.081
```

All 40 contaminated responses receive weights below 0.5, against 3 of the 360
clean ones:

```r
outlying <- seq_along(dat$y) %in% dat$idx_y_out
table(outlying, downweighted = fit$weights < 0.5)
#>         downweighted
#> outlying FALSE TRUE
#>    FALSE   357    3
#>    TRUE      0   40
```

```r
plot(dat$tps, dat$truth$beta_true, type = "l", lwd = 2,
     xlab = "s", ylab = "coefficient function")
lines(dat$tps, fit$beta_hat, col = "blue", lwd = 2)
```

In Case 2 the subject-level and visit-level components have their own
coefficient functions:

```r
dat2 <- simulate_data(case = "2", nI = 100, nV = 4,
                      contamination = "mixed", prop = 0.10, seed = 1)
lfpca2 <- rlfpca_estimate(Y = dat2$X, subject = dat2$subject,
                          Time = dat2$times, L = 0.95, smooth = TRUE,
                          tps = dat2$tps)
fit2 <- fit_rlfpcr(y = dat2$y, lfpca = lfpca2, subject = dat2$subject,
                   times = dat2$times, scalar_cov = dat2$scalar_cov)
round(c(ISE_B = sum(w * (fit2$beta_B - dat2$truth$beta_B_true)^2),
        ISE_U = sum(w * (fit2$beta_U - dat2$truth$beta_U_true)^2)), 3)
#> ISE_B ISE_U
#> 0.432 0.038
```

## Your own data

`rlfpca_estimate()` needs the curves as a matrix `Y` with one row per visit
and one column per grid point, together with the subject identifier and the
visit time of each row. The curves must be observed on a common dense grid
without missing values, and the rows must be sorted by subject (we also sort
the visits of a subject by time). Give the visit times on their original
scale, for example years since baseline, and the grid in `tps`. The
coefficient functions are on the scale of the integral over `tps`; with the
default grid 1, ..., `ncol(Y)` this is a sum over the grid points.
`fit_rlfpcr()` then needs the outcome of each row and, optionally, a matrix of
scalar covariates with one row per visit.

Useful parts of the results:

- `fit$beta_hat`, `fit$beta_B`, `fit$beta_0`, `fit$beta_1` and `fit$beta_U`:
  the coefficient functions on the grid;
- `fit$delta`, and `summary(fit$fit)` for the standard errors and p-values of
  the scalar covariates;
- `lfpca$curve_w` and `fit$weights`: the robust weights of the curves and of
  the responses; small values point to atypical scans or outcomes;
- `lfpca$N.X` and `lfpca$N.U`: the numbers of retained components, chosen by
  the explained-variance threshold `L` unless they are given.

The weights protect against sporadic atypical curves and outcomes. They do not
remove systematic differences between sites, scanners or batches, which should
be harmonized or modelled before the analysis.

## Reproducing the paper

The folder `replication/` holds the scripts behind the simulation study, the
sensitivity analyses and the ADNI application of the paper.
`replication/README.md` explains how to run them and which output gives each
table and figure. The scripts use the package code in `R/` of this repository.
The ADNI data cannot be redistributed and are not included; researchers can
apply for access at https://adni.loni.usc.edu.

## References

Gertheiss, J., Goldsmith, J., Crainiceanu, C. and Greven, S. (2013).
Longitudinal scalar-on-functions regression with application to tractography
data. Biostatistics, 14(3), 447-461.

Greven, S., Crainiceanu, C., Caffo, B. and Reich, D. (2010). Longitudinal
functional principal component analysis. Electronic Journal of Statistics, 4,
1022-1054.

## License

GPL-3. Questions and bug reports are welcome as GitHub issues or by email to
Ufuk Beyaztas (ufuk.beyaztas@marmara.edu.tr).
