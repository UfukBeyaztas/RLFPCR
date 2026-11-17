## Step 5: tables and figures of Section 5 and Table S11. Needs
## results/full_results.rds from step 1. The parts that use steps 2-4 are
## added when their files exist, so the script can be run at any time.
## Writes results/tables, results/figures and results/ADNI_key_numbers.txt,
## which collects the numbers quoted in the text.

source("config.R")
source("adni_functions.R")
load_rlfpcr_code()
make_dirs()
suppressPackageStartupMessages(library(ggplot2))

fr <- readRDS(file.path(OUT_DIR, "full_results.rds"))
tps <- fr$tps

## Report all coefficient functions on the scale COEF_SCALE (see config.R).
scale_betas <- function(f) {
  for (nm in intersect(c("beta_0", "beta_1", "beta_U", "beta_B", "beta"), names(f))) {
    f[[nm]] <- COEF_SCALE * f[[nm]]
  }
  f
}
fr$fits <- lapply(fr$fits, function(g) if (!is.null(g$method)) scale_betas(g) else lapply(g, scale_betas))
fr$fits_eqK <- lapply(fr$fits_eqK, function(g) lapply(g, scale_betas))
fr$fits_site <- lapply(fr$fits_site, scale_betas)
mc <- fr$main_config
LF <- fr$fits[[mc]]$LFPCR
RF <- fr$fits[[mc]]$RLFPCR
LP <- fr$fits$LPFR
y <- fr$y
t_med <- median(fr$times)
tq <- fr$t_quant
names(tq) <- c("Early visit", "Median visit", "Late visit")

key <- character(0)
note <- function(...) {
  txt <- paste0(...)
  key <<- c(key, txt)
  cat(txt, "\n")
}

theme_adni <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 8.5),
        axis.title = element_text(size = 10), axis.text = element_text(size = 8),
        strip.text = element_text(size = 9),
        strip.background = element_rect(fill = "grey90"))
save_fig <- function(p, name, width = 12, height = 10) {
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), plot = p,
         width = width, height = height, units = "cm")
}

## Table 4: sample summary
ss <- fr$sample_summary
write_tex_table(c(
  "\\begin{tabular}{ll}", "\\toprule", "Characteristic & Summary \\\\", "\\midrule",
  sprintf("Number of subjects & $%d$ \\\\", ss$n_subjects),
  sprintf("Number of visits & $%d$ \\\\", ss$n_visits),
  sprintf("Average number of visits per subject & $%.1f$ (range $%d$--$%d$) \\\\",
          ss$visits_mean, ss$visits_range[1], ss$visits_range[2]),
  sprintf("Number of ADNI sites & $%d$ \\\\", ss$n_sites),
  sprintf("Male sex & $%.1f\\%%$ \\\\", ss$male_pct),
  sprintf("Right-handed & $%.1f\\%%$ \\\\", ss$right_pct),
  sprintf("Baseline age, years & mean $%.1f$, SD $%.1f$, range $%.0f$--$%.0f$ \\\\",
          ss$age0_mean, ss$age0_sd, ss$age0_range[1], ss$age0_range[2]),
  sprintf("Years of education & mean $%.1f$, SD $%.1f$, range $%.0f$--$%.0f$ \\\\",
          ss$edu_mean, ss$edu_sd, ss$edu_range[1], ss$edu_range[2]),
  sprintf("Maximum follow-up, years & $%.1f$ \\\\", ss$followup_max),
  "\\bottomrule", "\\end{tabular}"), "tab_adni_sample_summary.tex")
note("Sample: ", ss$n_subjects, " subjects, ", ss$n_visits, " visits, ",
     sprintf("%.2f", ss$visits_mean), " visits/subject, ", ss$n_sites, " sites, ",
     ss$dup_subject_time, " duplicated (subject, time) pairs")

## Table 5: scalar effects
sc_rows <- list(LPFR = LP$scalar, LFPCR = LF$scalar, RLFPCR = RF$scalar)
sc_csv <- do.call(rbind, lapply(names(sc_rows), function(m) cbind(method = m, sc_rows[[m]])))
write.csv(sc_csv, file.path(TAB_DIR, "adni_scalar_effects.csv"), row.names = FALSE)
lines <- c("\\begin{tabular}{lrrrrrrrr}", "\\toprule",
           "& \\multicolumn{2}{c}{Age} & \\multicolumn{2}{c}{Gender} & \\multicolumn{2}{c}{Handedness} & \\multicolumn{2}{c}{Education} \\\\",
           "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}",
           "Method & Est. & $p$ & Est. & $p$ & Est. & $p$ & Est. & $p$ \\\\", "\\midrule")
for (m in names(sc_rows)) {
  s <- sc_rows[[m]][1:4, ]
  lines <- c(lines, paste0(m, " & ", paste(paste(fmt(s$estimate), fmt_p(s$p), sep = " & "),
                                           collapse = " & "), " \\\\"))
}
write_tex_table(c(lines, "\\bottomrule", "\\end{tabular}"), "tab_adni_scalar_effects.tex")

## Cross-validation summary (step 3)
cv_files <- list.files(file.path(RAW_DIR, "cv", mc), pattern = "^rep.*rds$", full.names = TRUE)
cv_sum <- NULL
if (length(cv_files) > 0) {
  cvp <- do.call(rbind, lapply(cv_files, function(f) readRDS(f)$pred))
  meths <- c("LPFR", "LFPCR", "RLFPCR")
  per_rep <- do.call(rbind, lapply(split(cvp, cvp$rep), function(z) {
    do.call(rbind, lapply(meths, function(m) {
      e <- z$y - z[[m]]
      data.frame(rep = z$rep[1], method = m, n = sum(!is.na(e)),
                 complete = nrow(z) == length(y),
                 MSE = mean(e^2, na.rm = TRUE), MAE = mean(abs(e), na.rm = TRUE),
                 MedAE = median(abs(e), na.rm = TRUE))
    }))
  }))
  write.csv(per_rep, file.path(TAB_DIR, "adni_cv_by_repeat.csv"), row.names = FALSE)
  full_reps <- per_rep[per_rep$complete, ]
  n_complete <- length(unique(full_reps$rep))
  if (n_complete == 0) full_reps <- per_rep
  cv_sum <- aggregate(cbind(MSE, MAE, MedAE) ~ method, data = full_reps, FUN = mean)
  cv_sd <- aggregate(MSE ~ method, data = full_reps, FUN = sd)
  cv_sum$MSE_sd <- cv_sd$MSE[match(cv_sum$method, cv_sd$method)]
  cv_sum$repeats <- n_complete
  note("CV (", CV_FOLDS, "-fold by subject): ", n_complete, " complete repeats out of ",
       CV_REPEATS, "; ", length(cv_files), " fold files",
       if (n_complete == 0) " (partial results)" else "")
  for (i in seq_len(nrow(cv_sum))) {
    note(sprintf("  %-7s CV MSE %.3f (sd over repeats %.3f), MAE %.3f, MedAE %.3f",
                 cv_sum$method[i], cv_sum$MSE[i], cv_sum$MSE_sd[i], cv_sum$MAE[i],
                 cv_sum$MedAE[i]))
  }
} else note("CV: no files yet (run 03_cv.R)")

## Table 6: in-sample MSE and cross-validated prediction errors
fit_tab <- data.frame(
  method = c("LPFR", "LFPCR", "RLFPCR"),
  KX = c(NA, LF$KX, RF$KX), KU = c(NA, LF$KU, RF$KU),
  MSE_in = c(mean((y - LP$fitted)^2), mean((y - LF$fitted)^2), mean((y - RF$fitted)^2)))
if (!is.null(cv_sum)) {
  fit_tab <- cbind(fit_tab, cv_sum[match(fit_tab$method, cv_sum$method),
                                   c("MSE", "MSE_sd", "MAE", "MedAE")])
  names(fit_tab)[5:8] <- c("CV_MSE", "CV_MSE_sd", "CV_MAE", "CV_MedAE")
}
write.csv(fit_tab, file.path(TAB_DIR, "adni_fit_summary.csv"), row.names = FALSE)
cvcell <- function(m, col) {
  if (is.null(cv_sum)) return("--")
  fmt(fit_tab[fit_tab$method == m, col])
}
lines <- c("\\begin{tabular}{lcccccc}", "\\toprule",
           "Method & $K_\\X$ & $K_U$ & In-sample MSE & CV MSE & CV MAE & CV median AE \\\\",
           "\\midrule")
for (m in fit_tab$method) {
  r <- fit_tab[fit_tab$method == m, ]
  lines <- c(lines, sprintf("%s & %s & %s & %s & %s & %s & %s \\\\", m,
                            ifelse(is.na(r$KX), "--", r$KX), ifelse(is.na(r$KU), "--", r$KU),
                            fmt(r$MSE_in), cvcell(m, "CV_MSE"), cvcell(m, "CV_MAE"),
                            cvcell(m, "CV_MedAE")))
}
write_tex_table(c(lines, "\\bottomrule", "\\end{tabular}"), "tab_adni_fit_summary.tex")
for (i in seq_len(nrow(fit_tab))) {
  note(sprintf("%-7s K_X=%s K_U=%s in-sample MSE %.3f", fit_tab$method[i],
               fit_tab$KX[i], fit_tab$KU[i], fit_tab$MSE_in[i]))
}

## Figure 4: coefficient functions (fig_3a, fig_3b and fig_3n)
comp_df <- do.call(rbind, lapply(list(LFPCR = LF, RLFPCR = RF), function(f) {
  rbind(data.frame(s = tps, beta = f$beta_0, Component = "beta[0](s)"),
        data.frame(s = tps, beta = f$beta_1, Component = "beta[1](s)"),
        data.frame(s = tps, beta = f$beta_U, Component = "beta[U](s)"))
}))
comp_df$Method <- rep(c("LFPCR", "RLFPCR"), each = 3 * length(tps))
comp_df$Component <- factor(comp_df$Component,
                            levels = c("beta[0](s)", "beta[1](s)", "beta[U](s)"))
p3a <- ggplot(comp_df, aes(x = s, y = beta, linetype = Method)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ Component, scales = "free_y", labeller = label_parsed, nrow = 1) +
  scale_linetype_manual(values = c(LFPCR = "solid", RLFPCR = "dashed")) +
  labs(x = "Functional domain", y = "Estimated coefficient",
       title = "Component-specific PET coefficient functions") +
  theme_adni
save_fig(p3a, "fig_3a")

## Subject-level coefficient beta_B(s): reconstructed from the joint
## subject-level dictionary {phi^(0)} U {phi^(1)} (coefficient of B_i0 + t B_i1
## when beta_0 = beta_1).
jb_df <- rbind(data.frame(s = tps, beta = LF$beta_B, Method = "LFPCR"),
               data.frame(s = tps, beta = RF$beta_B, Method = "RLFPCR"))
p3b <- ggplot(jb_df, aes(x = s, y = beta, linetype = Method)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_line(linewidth = 0.6) +
  scale_linetype_manual(values = c(LFPCR = "solid", RLFPCR = "dashed")) +
  labs(x = "Functional domain", y = expression(hat(beta)[B](s)),
       title = "Subject-level PET coefficient") +
  theme_adni
save_fig(p3b, "fig_3b")

## Single coefficient beta(s) from the full dictionary, compared with LPFR.
overall_df <- rbind(
  data.frame(s = tps, beta = LP$beta, Method = "LPFR"),
  data.frame(s = tps, beta = LF$beta, Method = "LFPCR"),
  data.frame(s = tps, beta = RF$beta, Method = "RLFPCR"))
overall_df$Method <- factor(overall_df$Method, levels = c("LPFR", "LFPCR", "RLFPCR"))
p3n <- ggplot(overall_df, aes(x = s, y = beta, linetype = Method)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_line(linewidth = 0.6) +
  scale_linetype_manual(values = c(LPFR = "solid", LFPCR = "dotted", RLFPCR = "dashed")) +
  labs(x = "Functional domain", y = expression(hat(beta)(s)),
       title = "Single-coefficient PET effects") +
  guides(linetype = guide_legend(nrow = 1)) +
  theme_adni
save_fig(p3n, "fig_3n")

nrm <- function(f) trapz_int(f^2, tps)
note(sprintf("Squared norms, LFPCR: beta_0 %.4g, beta_1 %.4g, beta_U %.4g, beta_B(t_med) %.4g",
             nrm(LF$beta_0), nrm(LF$beta_1), nrm(LF$beta_U), nrm(beta_B_at(LF, t_med))))
note(sprintf("Squared norms, RLFPCR: beta_0 %.4g, beta_1 %.4g, beta_U %.4g, beta_B(t_med) %.4g",
             nrm(RF$beta_0), nrm(RF$beta_1), nrm(RF$beta_U), nrm(beta_B_at(RF, t_med))))
note(sprintf("Squared norms, joint beta_B: LFPCR %.4g, RLFPCR %.4g; full beta: LFPCR %.4g, RLFPCR %.4g",
             nrm(LF$beta_B), nrm(RF$beta_B), nrm(LF$beta), nrm(RF$beta)))
note(sprintf("Squared norm, LPFR beta: %.4g", nrm(LP$beta)))
note(sprintf("Correlation of single coefficients: LPFR-LFPCR %.2f, LPFR-RLFPCR %.2f, LFPCR-RLFPCR %.2f",
             cor(LP$beta, LF$beta), cor(LP$beta, RF$beta), cor(LF$beta, RF$beta)))

## Figure 3: ADAS-Cog 13 trajectories
dat_y <- data.frame(y = y, subject = factor(fr$subject), time = fr$times)
p_adas <- ggplot(dat_y, aes(x = time, y = y, group = subject)) +
  geom_line(alpha = 0.20) +
  geom_point(alpha = 0.35, size = 0.8) +
  geom_smooth(aes(group = 1), method = "loess", formula = y ~ x, se = TRUE, linewidth = 0.9) +
  labs(x = "Visit time (years)", y = "ADAS-Cog 13") +
  theme_bw(base_size = 12)
save_fig(p_adas, "Fig_adas")

## Table S11: sensitivity of the coefficient functions
corr <- function(f, g) if (sd(f) > 0 && sd(g) > 0) cor(f, g) else NA
rel_isd <- function(f, g) isd(f, g, tps) / max(nrm(g), 1e-12)
setting_label <- function(cn) {
  switch(cn, main = "Main settings", first_version = "First-version settings", cn)
}
alt <- list()
for (cn in setdiff(names(fr$fits), c("LPFR", mc))) {
  alt[[setting_label(cn)]] <- fr$fits[[cn]]
}
for (cn in names(fr$fits_eqK)) {
  alt[[paste0(setting_label(cn), ", equal components")]] <- fr$fits_eqK[[cn]]
}
alt[["Site fixed effects"]] <- fr$fits_site[c("LFPCR", "RLFPCR")]
se_file <- file.path(OUT_DIR, "site_exclusion.rds")
if (file.exists(se_file)) {
  se <- readRDS(se_file)
  alt[[paste0("Without site ", se$site_removed)]] <- lapply(se$fits, scale_betas)
}

sens <- NULL
curves_sens <- NULL
add_curves <- function(f, lab, m) {
  rbind(data.frame(s = tps, beta = f$beta_0, Component = "beta[0](s)", Fit = lab, Method = m),
        data.frame(s = tps, beta = f$beta_1, Component = "beta[1](s)", Fit = lab, Method = m),
        data.frame(s = tps, beta = f$beta_U, Component = "beta[U](s)", Fit = lab, Method = m),
        data.frame(s = tps, beta = f$beta_B, Component = "beta[B](s)", Fit = lab, Method = m))
}
for (m in c("LFPCR", "RLFPCR")) {
  ref <- fr$fits[[mc]][[m]]
  curves_sens <- rbind(curves_sens, add_curves(ref, "Main", m))
  for (lab in names(alt)) {
    f <- alt[[lab]][[m]]
    curves_sens <- rbind(curves_sens, add_curves(f, lab, m))
    sens <- rbind(sens, data.frame(
      method = m, fit = lab, KX = f$KX, KU = f$KU,
      MSE_in = if (length(f$fitted) == length(y)) mean((y - f$fitted)^2) else NA,
      cor_betaB = corr(f$beta_B, ref$beta_B), risd_betaB = rel_isd(f$beta_B, ref$beta_B),
      cor_beta0 = corr(f$beta_0, ref$beta_0), cor_beta1 = corr(f$beta_1, ref$beta_1),
      cor_betaU = corr(f$beta_U, ref$beta_U),
      risd_beta0 = rel_isd(f$beta_0, ref$beta_0), risd_beta1 = rel_isd(f$beta_1, ref$beta_1),
      risd_betaU = rel_isd(f$beta_U, ref$beta_U),
      age = f$scalar$estimate[1], age_p = f$scalar$p[1],
      gender = f$scalar$estimate[2], gender_p = f$scalar$p[2]))
  }
}
between <- data.frame(
  method = "LFPCR vs RLFPCR", fit = "Main", KX = NA, KU = NA, MSE_in = NA,
  cor_betaB = corr(LF$beta_B, RF$beta_B), risd_betaB = rel_isd(LF$beta_B, RF$beta_B),
  cor_beta0 = corr(LF$beta_0, RF$beta_0), cor_beta1 = corr(LF$beta_1, RF$beta_1),
  cor_betaU = corr(LF$beta_U, RF$beta_U),
  risd_beta0 = rel_isd(LF$beta_0, RF$beta_0), risd_beta1 = rel_isd(LF$beta_1, RF$beta_1),
  risd_betaU = rel_isd(LF$beta_U, RF$beta_U),
  age = NA, age_p = NA, gender = NA, gender_p = NA)
sens <- rbind(sens, between)
write.csv(sens, file.path(TAB_DIR, "adni_sensitivity.csv"), row.names = FALSE)
lines <- c("\\begin{tabular}{llcccccccc}", "\\toprule",
           "Method & Fit & $K_\\X$ & $K_U$ & In-sample MSE & $\\rho(\\widehat\\beta_0)$ & $\\rho(\\widehat\\beta_1)$ & $\\rho(\\widehat\\beta_U)$ & $\\rho(\\widehat\\beta_B)$ & Age ($p$) \\\\",
           "\\midrule")
for (i in seq_len(nrow(sens))) {
  r <- sens[i, ]
  lines <- c(lines, sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
                            r$method, r$fit, ifelse(is.na(r$KX), "--", r$KX),
                            ifelse(is.na(r$KU), "--", r$KU),
                            ifelse(is.na(r$MSE_in), "--", fmt(r$MSE_in)),
                            fmt(r$cor_beta0, 2), fmt(r$cor_beta1, 2), fmt(r$cor_betaU, 2),
                            fmt(r$cor_betaB, 2),
                            ifelse(is.na(r$age), "--",
                                   paste0(fmt(r$age), " (", gsub("\\$", "", fmt_p(r$age_p)), ")"))))
}
write_tex_table(c(lines, "\\bottomrule", "\\end{tabular}"), "tab_adni_sensitivity.tex")
for (i in seq_len(nrow(sens))) {
  r <- sens[i, ]
  note(sprintf("Sensitivity %-7s %-28s K=(%s,%s) cor b0 %.2f b1 %.2f bU %.2f bB %.2f | rISD b0 %.3f b1 %.3f bU %.3f bB %.3f",
               r$method, r$fit, r$KX, r$KU, r$cor_beta0, r$cor_beta1, r$cor_betaU, r$cor_betaB,
               r$risd_beta0, r$risd_beta1, r$risd_betaU, r$risd_betaB))
}

curves_sens$Component <- factor(curves_sens$Component,
                                levels = c("beta[0](s)", "beta[1](s)", "beta[U](s)", "beta[B](s)"))
curves_sens$Fit <- factor(curves_sens$Fit, levels = c("Main", names(alt)))
p_sens <- ggplot(curves_sens, aes(x = s, y = beta, colour = Fit, linetype = Fit)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_line(linewidth = 0.5) +
  facet_grid(Method ~ Component, scales = "free_y", labeller = labeller(Component = label_parsed)) +
  labs(x = "Functional domain", y = "Estimated coefficient") +
  guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
  theme_adni
save_fig(p_sens, "fig_adni_sensitivity", width = 22, height = 12)

## Scalar effects with site fixed effects (CSV)
site_sc <- do.call(rbind, lapply(c("LPFR", "LFPCR", "RLFPCR"), function(m) {
  cbind(method = m, adjustment = "Site fixed effects", fr$fits_site[[m]]$scalar[1:4, ])
}))
write.csv(site_sc, file.path(TAB_DIR, "adni_scalar_effects_site.csv"), row.names = FALSE)

## Downweighting and age diagnostics
dw <- fr$downweight
write.csv(dw$by_site, file.path(TAB_DIR, "adni_downweight_by_site.csv"), row.names = FALSE)
dw_tab <- rbind(curves = dw$curves, visits = dw$visits)
write.csv(cbind(type = rownames(dw_tab), as.data.frame(dw_tab)),
          file.path(TAB_DIR, "adni_downweight_summary.csv"), row.names = FALSE)
note(sprintf("Downweighted curves (w < %.1f): %d (%.1f%%), %d subjects, %d sites; site test p = %.3f; top-3 sites hold %.0f%% of them and %.0f%% of visits",
             dw$cut, dw$curves["n"], dw$curves["pct"], dw$curves["subjects"], dw$curves["sites"],
             dw$curves["p_site"], 100 * dw$curves["top3.share_down"], 100 * dw$curves["top3.share_visits"]))
note(sprintf("Downweighted responses (w < %.1f): %d (%.1f%%), %d subjects, %d sites; site test p = %.3f; top-3 sites hold %.0f%% of them and %.0f%% of visits",
             dw$cut, dw$visits["n"], dw$visits["pct"], dw$visits["subjects"], dw$visits["sites"],
             dw$visits["p_site"], 100 * dw$visits["top3.share_down"], 100 * dw$visits["top3.share_visits"]))
note("Visits downweighted in both curve and response: ", dw$both)

as_ <- fr$age_scores
for (m in c("LFPCR", "RLFPCR")) {
  note(sprintf("%s: cor(baseline age, xi_k) = %s; R2(age ~ scores) = %.3f; VIF(age) = %.2f", m,
               paste(fmt(as_[[m]]$cor_baseline_age_xi, 3), collapse = ", "),
               as_[[m]]$r2_age_on_scores, as_[[m]]$vif_age))
}
note(sprintf("cor(age at visit, visit time) = %.3f", as_$cor_age_time))
note(sprintf("Prediction-code check (max abs difference to fitted values): LFPCR %.2e, RLFPCR %.2e (RMSE %.2e)",
             fr$pred_check["LFPCR"], fr$pred_check["RLFPCR"], fr$pred_check["RLFPCR_rmse"]))

## Table 7 and Figures 5 and 6: leave-one-subject-out stability (step 2)
cosd <- function(f, g) {
  den <- sqrt(nrm(f) * nrm(g))
  if (den == 0) return(NA)
  1 - trapz_int(f * g, tps) / den
}
loso_table <- function(cn) {
  files <- list.files(file.path(RAW_DIR, "loso", cn), pattern = "^subj_.*rds$", full.names = TRUE)
  if (length(files) == 0) return(NULL)
  full <- fr$fits[[cn]]
  out <- NULL
  for (f in files) {
    x <- readRDS(f)
    for (m in c("LFPCR", "RLFPCR")) {
      r <- scale_betas(x[[m]])
      g <- full[[m]]
      if (!isTRUE(r$success)) {
        out <- rbind(out, data.frame(subject_deleted = x$subject_deleted, method = m,
                                     success = FALSE, KX = NA, KU = NA, D_B = NA, D_U = NA,
                                     D_0 = NA, D_1 = NA, D_Bjoint = NA, D_full = NA, cos_B = NA,
                                     cos_U = NA, sse_pred = NA, n_pred = 0))
        next
      }
      bB <- beta_B_at(r, t_med)
      gB <- beta_B_at(g, t_med)
      out <- rbind(out, data.frame(
        subject_deleted = x$subject_deleted, method = m, success = TRUE,
        KX = r$KX, KU = r$KU,
        D_B = isd(bB, gB, tps), D_U = isd(r$beta_U, g$beta_U, tps),
        D_0 = isd(r$beta_0, g$beta_0, tps), D_1 = isd(r$beta_1, g$beta_1, tps),
        D_Bjoint = isd(r$beta_B, g$beta_B, tps),
        D_full = isd(r$beta, g$beta, tps),
        cos_B = cosd(bB, gB), cos_U = cosd(r$beta_U, g$beta_U),
        sse_pred = sum((x$y_test - r$pred_test)^2), n_pred = length(x$y_test)))
    }
  }
  norms <- sapply(c("LFPCR", "RLFPCR"), function(m) {
    g <- full[[m]]
    c(B = nrm(beta_B_at(g, t_med)), U = nrm(g$beta_U), b0 = nrm(g$beta_0),
      b1 = nrm(g$beta_1), Bjoint = nrm(g$beta_B), full = nrm(g$beta))
  })
  out$R_B <- out$D_B / norms["B", out$method]
  out$R_U <- out$D_U / norms["U", out$method]
  out$R_0 <- out$D_0 / norms["b0", out$method]
  out$R_1 <- out$D_1 / norms["b1", out$method]
  out$R_Bjoint <- out$D_Bjoint / norms["Bjoint", out$method]
  out$R_full <- out$D_full / norms["full", out$method]
  list(data = out, norms = norms, full_K = sapply(full, function(f) c(f$KX, f$KU)))
}

for (cn in LOSO_CONFIGS) {
  lt <- loso_table(cn)
  if (is.null(lt)) {
    note("LOSO (", cn, "): no files yet (run 02_loso.R)")
    next
  }
  dd <- lt$data
  write.csv(dd, file.path(TAB_DIR, paste0("adni_loso_", cn, "_by_subject.csv")), row.names = FALSE)
  ok_ids <- intersect(dd$subject_deleted[dd$method == "LFPCR" & dd$success],
                      dd$subject_deleted[dd$method == "RLFPCR" & dd$success])
  dv <- dd[dd$subject_deleted %in% ok_ids, ]
  wide <- function(col) {
    a <- dv[dv$method == "LFPCR", c("subject_deleted", col)]
    b <- dv[dv$method == "RLFPCR", c("subject_deleted", col)]
    merge(a, b, by = "subject_deleted", suffixes = c("_LFPCR", "_RLFPCR"))
  }
  summ <- NULL
  comps <- c(B = "beta_B(s,t_med)", U = "beta_U(s)", `0` = "beta_0(s)", `1` = "beta_1(s)",
             Bjoint = "beta_B(s) joint", full = "beta(s)")
  for (k in names(comps)) {
    w <- wide(paste0("R_", k))
    rl <- w[[2]]
    rr <- w[[3]]
    wr <- wide(paste0("D_", k))
    nm_row <- switch(k, B = "B", U = "U", `0` = "b0", `1` = "b1", Bjoint = "Bjoint", full = "full")
    summ <- rbind(summ, data.frame(
      component = comps[[k]], n = nrow(w),
      norm_LFPCR = lt$norms[nm_row, "LFPCR"], norm_RLFPCR = lt$norms[nm_row, "RLFPCR"],
      meanR_LFPCR = mean(rl), meanR_RLFPCR = mean(rr),
      medianR_LFPCR = median(rl), medianR_RLFPCR = median(rr),
      reduction_mean_pct = 100 * (1 - mean(rr) / mean(rl)),
      win_rate = mean(rr < rl),
      meanD_LFPCR = mean(wr[[2]]), meanD_RLFPCR = mean(wr[[3]])))
  }
  summ$meanCos_LFPCR <- NA
  summ$meanCos_RLFPCR <- NA
  for (k in c("B", "U")) {
    w <- wide(paste0("cos_", k))
    summ$meanCos_LFPCR[summ$component == comps[[k]]] <- mean(w[[2]], na.rm = TRUE)
    summ$meanCos_RLFPCR[summ$component == comps[[k]]] <- mean(w[[3]], na.rm = TRUE)
  }
  write.csv(summ, file.path(TAB_DIR, paste0("adni_loso_", cn, "_summary.csv")), row.names = FALSE)

  sameK <- tapply(seq_len(nrow(dv)), dv$method, function(ix) {
    m <- dv$method[ix[1]]
    mean(dv$KX[ix] == lt$full_K[1, m] & dv$KU[ix] == lt$full_K[2, m])
  })
  pred_mse <- tapply(seq_len(nrow(dv)), dv$method, function(ix) sum(dv$sse_pred[ix]) / sum(dv$n_pred[ix]))
  note("LOSO (", cn, "): ", length(ok_ids), " of ", length(unique(dd$subject_deleted)),
       " deletions with both fits successful")
  for (i in seq_len(nrow(summ))) {
    s <- summ[i, ]
    note(sprintf("  %-16s norm L %.4g R %.4g | mean R (x100) L %.3f R %.3f | median (x100) L %.3f R %.3f | reduction %.1f%% | win rate %.1f%%",
                 s$component, s$norm_LFPCR, s$norm_RLFPCR, 100 * s$meanR_LFPCR, 100 * s$meanR_RLFPCR,
                 100 * s$medianR_LFPCR, 100 * s$medianR_RLFPCR, s$reduction_mean_pct, 100 * s$win_rate))
  }
  note(sprintf("  share of LOSO fits with the full-data (K_X, K_U): LFPCR %.2f, RLFPCR %.2f",
               sameK["LFPCR"], sameK["RLFPCR"]))
  note(sprintf("  LOSO prediction MSE for the deleted subjects: LFPCR %.3f, RLFPCR %.3f",
               pred_mse["LFPCR"], pred_mse["RLFPCR"]))

  ## LaTeX table (relative distances x 100)
  tex_rows <- c(`0` = "$\\widehat\\beta_0(s)$", `1` = "$\\widehat\\beta_1(s)$",
                U = "$\\widehat\\beta_U(s)$", Bjoint = "$\\widehat\\beta_B(s)$",
                full = "$\\widehat\\beta(s)$")
  lines <- c("\\begin{tabular}{lrrrrrrrr}", "\\toprule",
             "& \\multicolumn{2}{c}{Full-data norm} & \\multicolumn{2}{c}{Mean} & \\multicolumn{2}{c}{Median} & & \\\\",
             "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}",
             "Coefficient & LFPCR & RLFPCR & LFPCR & RLFPCR & LFPCR & RLFPCR & Reduction & Win rate \\\\",
             "\\midrule")
  for (k in names(tex_rows)) {
    s <- summ[summ$component == comps[[k]], ]
    lines <- c(lines, sprintf("%s & %s & %s & %s & %s & %s & %s & %.1f\\%% & %.1f\\%% \\\\",
                              tex_rows[[k]], fmt(s$norm_LFPCR, 4), fmt(s$norm_RLFPCR, 4),
                              fmt(100 * s$meanR_LFPCR, 2), fmt(100 * s$meanR_RLFPCR, 2),
                              fmt(100 * s$medianR_LFPCR, 2), fmt(100 * s$medianR_RLFPCR, 2),
                              s$reduction_mean_pct, 100 * s$win_rate))
  }
  write_tex_table(c(lines, "\\bottomrule", "\\end{tabular}"),
                  paste0("tab_adni_loso_", cn, ".tex"))

  if (cn == mc) {
    long <- rbind(
      data.frame(subject_deleted = dv$subject_deleted, Method = dv$method,
                 Component = "Intercept coefficient", R = dv$R_0),
      data.frame(subject_deleted = dv$subject_deleted, Method = dv$method,
                 Component = "Visit-level coefficient", R = dv$R_U))
    p4a <- ggplot(long, aes(x = Method, y = log10(R + 1e-8))) +
      geom_boxplot(outlier.shape = NA, width = 0.55) +
      geom_jitter(width = 0.12, alpha = 0.45, size = 1.0) +
      facet_wrap(~ Component, scales = "free_y") +
      labs(x = NULL, y = expression(log[10]("normalized ISD")),
           title = "Leave-one-subject-out relative stability on the log scale") +
      theme_adni + theme(legend.position = "none")
    save_fig(p4a, "Fig_4a")

    wb <- wide("R_0")
    wu <- wide("R_U")
    sc_df <- rbind(data.frame(Component = "Intercept coefficient",
                              LFPCR = wb[[2]], RLFPCR = wb[[3]]),
                   data.frame(Component = "Visit-level coefficient",
                              LFPCR = wu[[2]], RLFPCR = wu[[3]]))
    p4b <- ggplot(sc_df, aes(x = LFPCR, y = RLFPCR)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      geom_point(alpha = 0.7, size = 1.4) +
      facet_wrap(~ Component, scales = "free") +
      labs(x = "LFPCR normalized distance", y = "RLFPCR normalized distance",
           title = "Subject-level comparison of relative perturbations") +
      theme_adni
    save_fig(p4b, "Fig_4b")

    tot <- data.frame(subject_deleted = wb$subject_deleted,
                      LFPCR = wb[[2]] + wu[[2]][match(wb$subject_deleted, wu$subject_deleted)],
                      RLFPCR = wb[[3]] + wu[[3]][match(wb$subject_deleted, wu$subject_deleted)])
    tot <- tot[order(-abs(tot$LFPCR - tot$RLFPCR)), ][1:min(20, nrow(tot)), ]
    top_long <- rbind(data.frame(subject_deleted = tot$subject_deleted, Method = "LFPCR", T = tot$LFPCR),
                      data.frame(subject_deleted = tot$subject_deleted, Method = "RLFPCR", T = tot$RLFPCR))
    top_long$subject_deleted <- factor(top_long$subject_deleted, levels = tot$subject_deleted)
    p5 <- ggplot(top_long, aes(x = subject_deleted, y = T, fill = Method)) +
      geom_col(position = "dodge") +
      scale_fill_grey(start = 0.3, end = 0.7) +
      labs(x = "Deleted subject", y = "Total normalized distance",
           title = "Largest subject-deletion effects") +
      theme_adni + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
    save_fig(p5, "Fig_5")
  }
}

writeLines(key, file.path(OUT_DIR, "ADNI_key_numbers.txt"))
cat("\nTables in", TAB_DIR, "\nFigures in", FIG_DIR, "\n")
