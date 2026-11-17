## Step 2: leave-one-subject-out analysis (Section 5.5 of the paper). For
## every subject, LFPCR and RLFPCR are refitted without that subject, and the
## coefficient functions and the predictions for the deleted subject are
## saved in raw/loso/<setting>/subj_XXX.rds.
##
## Finished subjects are skipped when the script is run again. An empty file
## named STOP in this folder stops the workers after their current subject
## (delete it before starting again). Takes about 75 minutes with 12 workers
## for each setting in LOSO_CONFIGS.

source("config.R")
source("adni_functions.R")
load_rlfpcr_code()
make_dirs()

d <- load_adni()
subjects <- sort(unique(d$subject))

tasks <- list()
for (cn in LOSO_CONFIGS) {
  done <- file.exists(file.path(RAW_DIR, "loso", cn, sprintf("subj_%03d.rds", subjects)))
  for (i in subjects[!done]) tasks[[length(tasks) + 1]] <- list(i_subj = i, cfg_name = cn)
}
cat(length(tasks), "LOSO fits to run with", N_WORKERS, "workers\n")

t0 <- Sys.time()
status <- run_parallel(tasks, "run_loso_task")
print(table(status))
print(Sys.time() - t0)

for (cn in LOSO_CONFIGS) {
  n_done <- sum(file.exists(file.path(RAW_DIR, "loso", cn, sprintf("subj_%03d.rds", subjects))))
  cat("Configuration", cn, ":", n_done, "of", length(subjects), "subjects finished\n")
}
