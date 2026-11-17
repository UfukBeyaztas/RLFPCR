## Step 1: full-data fits of LPFR, LFPCR and RLFPCR, the sensitivity fits of
## Table S11 (Additional file 1), and the downweighting and age diagnostics.
## Takes about 20 minutes and writes results/full_results.rds.

source("config.R")
source("adni_functions.R")
load_rlfpcr_code()
make_dirs()

d <- load_adni()
cat("Visits:", d$n, " subjects:", length(unique(d$subject)),
    " sites:", length(unique(d$site)), "\n")

## Sample description (Table 4)
first <- !duplicated(d$subject)
J <- as.numeric(table(d$subject))
age0 <- d$Z[first, "Age"]
edu0 <- d$Z[first, "Education"]
sample_summary <- list(
  n_subjects = sum(first),
  n_visits = d$n,
  visits_mean = mean(J), visits_range = range(J),
  male_pct = 100 * mean(d$gender[first] == "Male"),
  right_pct = 100 * mean(d$hand[first] == "Right"),
  age0_mean = mean(age0), age0_sd = sd(age0), age0_range = range(age0),
  edu_mean = mean(edu0), edu_sd = sd(edu0), edu_range = range(edu0),
  followup_max = max(d$times),
  n_sites = length(unique(d$site)),
  dup_subject_time = sum(duplicated(cbind(d$subject, d$times)))
)
str(sample_summary)

## Main and first-version settings
fits <- list()
cat("LPFR ...\n")
fits$LPFR <- fit_lpfr_adni(d)
heavy_main <- NULL
for (cn in names(CONFIGS)) {
  cat("Configuration", cn, ": LFPCR ...\n")
  f1 <- fit_lfpcr_adni(d, CONFIGS[[cn]]$lfpcr)
  cat("Configuration", cn, ": RLFPCR ...\n")
  f2 <- fit_rlfpcr_adni(d, CONFIGS[[cn]]$rlfpcr)
  cat(sprintf("  LFPCR  K_X=%d K_U=%d MSE=%.3f\n", f1$KX, f1$KU, mean((d$y - f1$fitted)^2)))
  cat(sprintf("  RLFPCR K_X=%d K_U=%d MSE=%.3f\n", f2$KX, f2$KU, mean((d$y - f2$fitted)^2)))
  if (cn == MAIN_CONFIG) heavy_main <- list(LFPCR = f1, RLFPCR = f2)
  fits[[cn]] <- list(LFPCR = light(f1), RLFPCR = light(f2))
}

## Check of the prediction code for new subjects. Applied to the training
## data, it reproduces the population-level fitted values, exactly for LFPCR
## and closely for RLFPCR, whose curve weights are recomputed.
chk_l <- predict_lfpcr_new(heavy_main$LFPCR, d)
chk_r <- predict_rlfpcr_new(heavy_main$RLFPCR, d)
pred_check <- c(
  LFPCR = max(abs(chk_l - heavy_main$LFPCR$fitted_pop)),
  RLFPCR = max(abs(chk_r - heavy_main$RLFPCR$fitted_pop)),
  RLFPCR_rmse = sqrt(mean((chk_r - heavy_main$RLFPCR$fitted_pop)^2)))
print(pred_check)
rm(chk_l, chk_r)

## Equal numbers of components. When the two methods retain different
## numbers of components, LFPCR is refitted with the numbers of RLFPCR and
## RLFPCR with those of LFPCR.
fits_eqK <- list()
for (cn in names(CONFIGS)) {
  kl <- c(fits[[cn]]$LFPCR$KX, fits[[cn]]$LFPCR$KU)
  kr <- c(fits[[cn]]$RLFPCR$KX, fits[[cn]]$RLFPCR$KU)
  if (all(kl == kr)) {
    cat("Configuration", cn, ": both methods retain K_X =", kl[1], "and K_U =", kl[2], "\n")
    next
  }
  cat("Configuration", cn, ": equal-component refits ...\n")
  fits_eqK[[cn]] <- list(
    LFPCR = light(fit_lfpcr_adni(d, CONFIGS[[cn]]$lfpcr, N.X = kr[1], N.U = kr[2])),
    RLFPCR = light(fit_rlfpcr_adni(d, CONFIGS[[cn]]$rlfpcr, N.X = kl[1], N.U = kl[2])))
}

## Site fixed effects added to the scalar covariates of all three methods.
cfg <- CONFIGS[[MAIN_CONFIG]]
Zsite <- cbind(d$Z, site_dummies(d$site))
cat("Site fixed effects ...\n")
fits_site <- list(
  LPFR = fit_lpfr_adni(d, Z = Zsite),
  LFPCR = light(fit_lfpcr_adni(d, cfg$lfpcr, Z = Zsite)),
  RLFPCR = light(fit_rlfpcr_adni(d, cfg$rlfpcr, Z = Zsite)))
fits_site$LPFR$pred_parts <- NULL

## Downweighted curves and visits and their distribution over the sites.
rob <- fits[[MAIN_CONFIG]]$RLFPCR
down_curve <- rob$curve_w < DOWNWEIGHT_CUT
down_resp <- rob$resp_w < DOWNWEIGHT_CUT
by_site <- data.frame(
  site = sort(unique(d$site)),
  subjects = as.numeric(tapply(d$subject, d$site, function(z) length(unique(z)))),
  visits = as.numeric(table(d$site)),
  curves_down = as.numeric(tapply(down_curve, d$site, sum)),
  visits_down = as.numeric(tapply(down_resp, d$site, sum)))
site_test <- function(flag) {
  if (sum(flag) == 0 || all(flag)) return(NA)
  set.seed(1)
  fisher.test(table(d$site, flag), simulate.p.value = TRUE, B = 1e5)$p.value
}
p_curve <- site_test(down_curve)
p_resp <- site_test(down_resp)
top_share <- function(k) {
  o <- order(-k)[1:3]
  c(share_down = sum(k[o]) / max(sum(k), 1), share_visits = sum(by_site$visits[o]) / d$n)
}
downweight <- list(
  cut = DOWNWEIGHT_CUT,
  curves = c(n = sum(down_curve), pct = 100 * mean(down_curve),
             subjects = length(unique(d$subject[down_curve])),
             sites = length(unique(d$site[down_curve])), p_site = p_curve,
             top3 = top_share(by_site$curves_down)),
  visits = c(n = sum(down_resp), pct = 100 * mean(down_resp),
             subjects = length(unique(d$subject[down_resp])),
             sites = length(unique(d$site[down_resp])), p_site = p_resp,
             top3 = top_share(by_site$visits_down)),
  both = sum(down_curve & down_resp),
  by_site = by_site)
print(downweight[c("curves", "visits", "both")])

## Age and the subject-level scores.
age_scores <- list()
for (m in c("LFPCR", "RLFPCR")) {
  f <- fits[[MAIN_CONFIG]][[m]]
  xi <- as.matrix(f$xi)
  idx <- match(d$subject, sort(unique(d$subject)))
  tstd <- (d$times - f$time_mean) / f$time_sd
  S <- cbind(xi[idx, , drop = FALSE], xi[idx, , drop = FALSE] * tstd, f$zeta)
  r2_scores <- summary(lm(d$Z[, "Age"] ~ S))$r.squared
  r2_all <- summary(lm(d$Z[, "Age"] ~ S + d$Z[, -1] + d$times))$r.squared
  age_scores[[m]] <- list(
    cor_baseline_age_xi = as.numeric(cor(age0, xi)),
    r2_age_on_scores = r2_scores,
    vif_age = 1 / (1 - r2_all))
}
age_scores$cor_age_time <- cor(d$Z[, "Age"], d$times)
print(age_scores)

## Save
full_results <- list(
  created = Sys.time(),
  sample_summary = sample_summary,
  times = d$times, y = d$y, subject = d$subject, site = d$site, tps = d$tps,
  t_quant = quantile(d$times, TIME_PROBS),
  configs = CONFIGS, main_config = MAIN_CONFIG, gram_tol = GRAM_TOL,
  fits = fits, fits_eqK = fits_eqK, fits_site = fits_site,
  downweight = downweight,
  age_scores = age_scores, pred_check = pred_check,
  session = sessionInfo())
saveRDS(full_results, file.path(OUT_DIR, "full_results.rds"))
cat("Saved", file.path(OUT_DIR, "full_results.rds"), "\n")
