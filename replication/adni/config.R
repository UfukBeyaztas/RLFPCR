## Settings shared by the ADNI scripts. Every script sources this file and
## must be run from the adni folder (setwd() to it, or call Rscript there).

ADNI_DIR  <- normalizePath(".")
PKG_DIR   <- normalizePath("../../R")          # R code of the RLFPCR package
SIM_DIR   <- normalizePath("../simulation")    # classical LFPCA and LFPCR

## The ADNI data cannot be redistributed, so this file is not part of the
## repository. README.md describes its layout.
DATA_FILE <- file.path(ADNI_DIR, "data", "longitudinal_only.xlsx")

OUT_DIR   <- file.path(ADNI_DIR, "results")
RAW_DIR   <- file.path(ADNI_DIR, "raw")        # one file per LOSO or CV fit
TAB_DIR   <- file.path(OUT_DIR, "tables")
FIG_DIR   <- file.path(OUT_DIR, "figures")
STOP_FILE <- file.path(ADNI_DIR, "STOP")       # create it to stop the workers

## Parallel workers for 02_loso.R and 03_cv.R; each needs about 2 GB.
N_WORKERS <- 12

## Grid of the 1000 ranked PET features on [0, 1].
TPS <- seq(0, 1, length.out = 1000)

## The fits use integrals over [0, 1]. Tables and figures report the
## coefficient functions multiplied by the grid spacing 1/999, i.e. on the
## scale of a sum over the 1000 features, which is the scale of refund::lpfr.
## Correlations and normalized distances do not depend on this choice.
COEF_SCALE <- 1 / (length(TPS) - 1)

## Relative eigenvalue threshold of the Gram correction.
GRAM_TOL <- 0.05

## Model settings. "main" is used throughout Section 5 of the paper. The
## classical mean surface uses k_time = 5 in time, because the visit times
## take only seven distinct values (0, 0.5, 1, 1.5, 2, 3 and 4 years).
## "first_version" (L = 0.90 and four basis functions for LFPCR) is the
## first-version setting of Table S11 in Additional file 1.
CONFIGS <- list(
  main = list(
    lfpcr  = list(L = 0.95, bf = 8, k_time = 5),
    rlfpcr = list(L = 0.95, bf = 8, df_time_mean = 5, mean_iter = 3,
                  clean_cells = TRUE, clean_c = 4.5, score_iter = 8,
                  max_iter = 4, psi = "bisquare", c = 4.685)
  ),
  first_version = list(
    lfpcr  = list(L = 0.90, bf = 4, k_time = 4),
    rlfpcr = list(L = 0.95, bf = 8, df_time_mean = 5, mean_iter = 3,
                  clean_cells = TRUE, clean_c = 4.5, score_iter = 8,
                  max_iter = 4, psi = "bisquare", c = 4.685)
  )
)
MAIN_CONFIG <- "main"

## Settings with a leave-one-subject-out analysis. Each one takes about
## 75 minutes with 12 workers.
LOSO_CONFIGS <- c("main")

## Subject-level cross-validation.
CV_FOLDS   <- 10
CV_REPEATS <- 3
CV_SEED    <- 20261113

## Curves and visits with a weight below this value count as downweighted.
DOWNWEIGHT_CUT <- 0.5

## Visit-time quantiles at which beta_B(s, t) is evaluated.
TIME_PROBS <- c(0.10, 0.50, 0.90)
