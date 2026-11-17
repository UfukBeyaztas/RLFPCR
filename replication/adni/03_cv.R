## Step 3: subject-level cross-validation (Table 6). The folds are formed by
## subjects, so all visits of a subject are predicted by a model that has
## not seen that subject: the subject effect is set to zero and the PET
## scores are computed from the training eigenfunctions. Each (repeat, fold)
## task fits LPFR, LFPCR and RLFPCR and is saved in
## raw/cv/main/repRR_foldKK.rds.
##
## Restarting and stopping work as in 02_loso.R. The 3 x 10 tasks take about
## 20 minutes with 12 workers.

source("config.R")
source("adni_functions.R")
load_rlfpcr_code()
make_dirs()

grid <- expand.grid(k = seq_len(CV_FOLDS), r = seq_len(CV_REPEATS))
done <- file.exists(file.path(RAW_DIR, "cv", MAIN_CONFIG,
                              sprintf("rep%02d_fold%02d.rds", grid$r, grid$k)))
tasks <- lapply(which(!done), function(j) list(r = grid$r[j], k = grid$k[j]))
cat(length(tasks), "CV tasks to run with", N_WORKERS, "workers\n")

t0 <- Sys.time()
status <- run_parallel(tasks, "run_cv_task")
print(table(status))
print(Sys.time() - t0)
