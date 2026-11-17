## Step 4: LFPCR and RLFPCR refitted without the site that holds most of the
## downweighted PET curves, and compared with the full-data fits (Section 5.6
## and Table S11). Needs results/full_results.rds from step 1. Takes about
## 5 minutes and writes results/site_exclusion.rds and
## results/tables/adni_site_exclusion.csv.

source("config.R")
source("adni_functions.R")
load_rlfpcr_code()
make_dirs()

fr <- readRDS(file.path(OUT_DIR, "full_results.rds"))
d <- load_adni()
bs <- fr$downweight$by_site
drop_site <- bs$site[which.max(bs$curves_down)]
keep <- which(d$site != drop_site)
cat("Site", drop_site, ":", sum(d$site == drop_site), "visits,",
    length(unique(d$subject[d$site == drop_site])), "subjects,",
    bs$curves_down[bs$site == drop_site], "of", sum(bs$curves_down),
    "downweighted curves\n")

dk <- subset_adni(d, keep)
cfg <- CONFIGS[[MAIN_CONFIG]]
fits <- list(LFPCR = light(fit_lfpcr_adni(dk, cfg$lfpcr)),
             RLFPCR = light(fit_rlfpcr_adni(dk, cfg$rlfpcr)))

t_med <- median(fr$times)
corr <- function(f, g) cor(f, g)
rel_isd <- function(f, g) isd(f, g, TPS) / trapz_int(g^2, TPS)
out <- NULL
for (m in c("LFPCR", "RLFPCR")) {
  f <- fits[[m]]
  g <- fr$fits[[MAIN_CONFIG]][[m]]
  out <- rbind(out, data.frame(
    method = m, site_removed = drop_site, KX = f$KX, KU = f$KU,
    cor_beta0 = corr(f$beta_0, g$beta_0), cor_beta1 = corr(f$beta_1, g$beta_1),
    cor_betaU = corr(f$beta_U, g$beta_U),
    cor_betaB_med = corr(beta_B_at(f, t_med), beta_B_at(g, t_med)),
    risd_beta0 = rel_isd(f$beta_0, g$beta_0), risd_beta1 = rel_isd(f$beta_1, g$beta_1),
    risd_betaU = rel_isd(f$beta_U, g$beta_U),
    risd_betaB_med = rel_isd(beta_B_at(f, t_med), beta_B_at(g, t_med)),
    curves_down = if (m == "RLFPCR") sum(f$curve_w < DOWNWEIGHT_CUT) else NA,
    age = f$scalar$estimate[1], age_p = f$scalar$p[1],
    gender = f$scalar$estimate[2], gender_p = f$scalar$p[2]))
}
print(out, digits = 3)
write.csv(out, file.path(TAB_DIR, "adni_site_exclusion.csv"), row.names = FALSE)
saveRDS(list(site_removed = drop_site, fits = fits, summary = out),
        file.path(OUT_DIR, "site_exclusion.rds"))
