## Sensitivity of RLFPCR to its tuning constants and its behaviour under
## increasing contamination (Tables S7-S10 of Additional file 1).
##
## Part A refits RLFPCR with the Tukey constant c in {3.5, 4.685, 6} (curves
## and responses) and the cellwise threshold in {3.5, 4.5, 6, none} on the
## clean, mixed 10% and local-spike data of the main study (same seeds), for
## I = 100 and 250 with m = 4.
## Part B fits LFPCR, Bag-LFPCR and RLFPCR under mixed contamination with
## proportion p from 0 to 0.4, for I = 100 and m = 4.
##
## Run from this folder with Rscript run_sensitivity.R (or setwd() to this
## folder and source() it). Finished tasks are skipped when the script is run
## again, and an empty file named STOP in RAW_DIR stops the workers.

PKG_DIR   <- normalizePath("../../R")          # R code of the RLFPCR package
SIM_DIR   <- normalizePath("../simulation")    # sim_functions.R
SENS_DIR  <- normalizePath(".")
RAW_DIR   <- file.path(SENS_DIR, "raw")
OUT_DIR   <- file.path(SENS_DIR, "results")
N_WORKERS <- 28
N_REP     <- 100

## RLFPCR variants of part A and contamination proportions of part B.
VARIANTS <- list(
  "c = 3.5"       = list(c = 3.5,   clean_cells = TRUE,  clean_c = 4.5),
  "c = 4.685"     = list(c = 4.685, clean_cells = TRUE,  clean_c = 4.5),
  "c = 6"         = list(c = 6,     clean_cells = TRUE,  clean_c = 4.5),
  "cells 3.5 MAD" = list(c = 4.685, clean_cells = TRUE,  clean_c = 3.5),
  "cells 6 MAD"   = list(c = 4.685, clean_cells = TRUE,  clean_c = 6),
  "no cell step"  = list(c = 4.685, clean_cells = FALSE, clean_c = 4.5)
)
P_GRID <- c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40)

load_sim_code <- function(pkg_dir, sim_dir) {
  suppressPackageStartupMessages({
    library(mgcv)
    library(splines)
    library(refund)
    library(rainbow)
  })
  for (f in list.files(pkg_dir, pattern = "[.]R$", full.names = TRUE)) source(f)
  source(file.path(sim_dir, "sim_functions.R"))
  invisible(TRUE)
}

## One task per data set; each task fits all methods of its part.
taskA <- expand.grid(case = 1:2, design = c(1, 3),
                     scenario = c("clean", "mixed10", "spike"),
                     rep = seq_len(N_REP), stringsAsFactors = FALSE)
taskA$part <- "A"
taskA$p <- NA
taskB <- expand.grid(case = 1:2, design = 1, p = P_GRID, rep = seq_len(N_REP))
taskB$part <- "B"
taskB$scenario <- NA
tasks <- rbind(taskA[, c("part", "case", "design", "scenario", "p", "rep")],
               taskB[, c("part", "case", "design", "scenario", "p", "rep")])
tasks$file <- with(tasks, file.path(RAW_DIR, ifelse(
  part == "A",
  sprintf("A_case%d_d%d_%s_rep%04d.rds", case, design, scenario, rep),
  sprintf("B_case%d_p%03d_rep%04d.rds", case, round(100 * p), rep))))

run_task <- function(k) {
  tk <- tasks[k, ]
  if (file.exists(file.path(RAW_DIR, "STOP"))) return("stopped")
  if (file.exists(tk$file)) return("exists")
  dt <- design_table()[tk$design, ]
  if (tk$part == "A") {
    sc <- scenario_table()
    sc <- sc[sc$scenario == tk$scenario, ]
    d <- simulate_data(case = as.character(tk$case), nI = dt$nI, nV = dt$nV,
                       contamination = sc$contamination, prop = max(sc$prop, 0.01),
                       seed = sim_seed(tk$case, tk$design, tk$rep))
    fits <- c(list(LFPCR = function(d, w) m_lfpcr(d, w)),
              lapply(VARIANTS, function(v) {
                force(v)
                function(d, w) m_rlfpcr(d, w, c = v$c, clean_cells = v$clean_cells,
                                        clean_c = v$clean_c)
              }))
  } else {
    d <- simulate_data(case = as.character(tk$case), nI = dt$nI, nV = dt$nV,
                       contamination = if (tk$p == 0) "none" else "mixed",
                       prop = max(tk$p, 0.01),
                       seed = sim_seed(tk$case, tk$design, tk$rep))
    fits <- list(LFPCR = function(d, w) m_lfpcr(d, w),
                 "Bag-LFPCR" = function(d, w) m_screen_lfpcr(d, w, "bagplot"),
                 RLFPCR = function(d, w) m_rlfpcr(d, w))
  }
  w <- trap_weights(d$tps)
  n <- length(d$y)
  frac_curves <- length(unique(c(d$idx_subject_out, d$idx_visit_out))) / n
  frac_resp <- length(d$idx_y_out) / n
  rows <- lapply(names(fits), function(mn) {
    out <- tryCatch(fits[[mn]](d, w), error = function(e) e)
    base <- data.frame(part = tk$part, case = tk$case, nI = dt$nI, nV = dt$nV,
                       scenario = tk$scenario, p = tk$p, rep = tk$rep, method = mn,
                       frac_curves = frac_curves, frac_resp = frac_resp,
                       stringsAsFactors = FALSE)
    if (inherits(out, "error")) {
      return(cbind(base, mse = NA, mse_delta = NA, ise = NA, ise_B = NA, ise_U = NA,
                   n_curve_drop = NA, n_resp_drop = NA, error = conditionMessage(out)))
    }
    met <- evaluate_fit(out, d, w)
    cbind(base, mse = met[["mse"]], mse_delta = met[["mse_delta"]], ise = met[["ise"]],
          ise_B = met[["ise_B"]], ise_U = met[["ise_U"]],
          n_curve_drop = if (is.null(out$n_curve_drop)) NA else out$n_curve_drop,
          n_resp_drop = if (is.null(out$n_resp_drop)) NA else out$n_resp_drop,
          error = "")
  })
  res <- do.call(rbind, rows)
  tmp <- paste0(tk$file, ".tmp", Sys.getpid())
  saveRDS(res, tmp)
  file.rename(tmp, tk$file)
  "ok"
}

## Run the missing tasks in parallel.
dir.create(RAW_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
load_sim_code(PKG_DIR, SIM_DIR)
todo <- which(!file.exists(tasks$file))
cat(sprintf("%s: %d of %d tasks done, %d to run.\n", format(Sys.time()),
            nrow(tasks) - length(todo), nrow(tasks), length(todo)))
if (length(todo) > 0) {
  library(parallel)
  cl <- makePSOCKcluster(min(N_WORKERS, length(todo)))
  status <- tryCatch({
    clusterExport(cl, c("tasks", "RAW_DIR", "VARIANTS", "PKG_DIR", "SIM_DIR",
                        "load_sim_code", "run_task"))
    invisible(clusterEvalQ(cl, {
      Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
      load_sim_code(PKG_DIR, SIM_DIR)
    }))
    unlist(parLapplyLB(cl, todo, run_task, chunk.size = 1))
  }, finally = stopCluster(cl))
  print(table(status))
}

## Collect the task files, summarize them and write the tables.
files <- tasks$file[file.exists(tasks$file)]
all <- do.call(rbind, lapply(files, readRDS))
saveRDS(all, file.path(OUT_DIR, "sensitivity_metrics_all.rds"))
cat("Failed fits:", sum(all$error != ""), "\n")
all$scenario[is.na(all$scenario)] <- "-"
all$p[is.na(all$p)] <- -1
grp <- split(all, list(all$part, all$case, all$nI, all$nV, all$scenario, all$p, all$method),
             drop = TRUE)
mn <- function(x) mean(x, na.rm = TRUE)
md <- function(x) median(x, na.rm = TRUE)
agg <- do.call(rbind, lapply(grp, function(z) {
  data.frame(z[1, c("part", "case", "nI", "nV", "scenario", "p", "method")],
             n = sum(!is.na(z$mse)),
             mse = mn(z$mse), mse_delta = mn(z$mse_delta), ise = mn(z$ise),
             ise_B = mn(z$ise_B), ise_U = mn(z$ise_U),
             mse_median = md(z$mse), ise_median = md(z$ise),
             ise_B_median = md(z$ise_B), ise_U_median = md(z$ise_U),
             n_curve_drop = mn(z$n_curve_drop), n_resp_drop = mn(z$n_resp_drop),
             frac_curves = mn(z$frac_curves), frac_resp = mn(z$frac_resp),
             stringsAsFactors = FALSE)
}))
rownames(agg) <- NULL
write.csv(agg, file.path(OUT_DIR, "sensitivity_summary.csv"), row.names = FALSE)
source(file.path(SENS_DIR, "make_sensitivity_tables.R"))
make_sensitivity_tables(agg, OUT_DIR)
