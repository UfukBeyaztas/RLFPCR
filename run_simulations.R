## Simulation study of Section 4 of the paper.
##
## Run from this folder, either with
##   Rscript run_simulations.R
## or in R after setwd() to this folder, with source("run_simulations.R").
##
## Each replication is saved to its own file in RAW_DIR as soon as it
## finishes. If a run is interrupted, start it again with the same command;
## finished replications are skipped. To stop the workers after their current
## replication, create an empty file named STOP in RAW_DIR (delete it before
## starting again).

PKG_DIR <- normalizePath("../../R")    # R code of the RLFPCR package
SIM_DIR <- normalizePath(".")
RAW_DIR <- file.path(SIM_DIR, "raw")   # 20,000 small files

N_REP     <- 250
N_WORKERS <- 28      # we used 28 workers on a 32-core PC
CASES     <- c(1, 2)
DESIGNS   <- 1:4     # 1: I = 100, m = 4; 2: I = 100, m = 8;
                     # 3: I = 250, m = 4; 4: I = 250, m = 8
SCENARIOS <- c("clean", "mixed05", "mixed10", "spike", "drift", "cluster",
               "covariate", "batch20", "batch30", "instrument")
METHODS   <- c("LPFR", "LFPCR", "MAD-LFPCR", "Bag-LFPCR", "KFPCA-LS", "rLMM",
               "RLFPCR")
CURVE_DESIGNS <- 3   # estimated coefficient curves are kept for these designs

load_sim_code <- function(pkg_dir, sim_dir) {
  suppressPackageStartupMessages({
    library(mgcv)
    library(splines)
    library(refund)
    library(robustlmm)
    library(rainbow)
  })
  for (f in list.files(pkg_dir, pattern = "[.]R$", full.names = TRUE)) source(f)
  source(file.path(sim_dir, "sim_functions.R"))
  invisible(TRUE)
}

load_sim_code(PKG_DIR, SIM_DIR)
dir.create(RAW_DIR, showWarnings = FALSE, recursive = TRUE)

cell_name <- function(case, design, scenario) {
  dt <- design_table()[design, ]
  sprintf("case%d_I%d_m%d_%s", case, dt$nI, dt$nV, scenario)
}

## One task per replication. Replication 1 of every cell comes first, then
## replication 2 and so on, so that a partial run gives all cells about the
## same number of replications.
tasks <- expand.grid(scenario = SCENARIOS, design = DESIGNS, case = CASES,
                     rep = seq_len(N_REP), stringsAsFactors = FALSE)
tasks <- tasks[order(tasks$rep, tasks$case, tasks$design,
                     match(tasks$scenario, SCENARIOS)), ]
tasks$cell <- mapply(cell_name, tasks$case, tasks$design, tasks$scenario)
tasks$file <- file.path(RAW_DIR, tasks$cell, sprintf("rep_%04d.rds", tasks$rep))
for (cl_dir in unique(file.path(RAW_DIR, tasks$cell))) {
  dir.create(cl_dir, showWarnings = FALSE)
}
todo <- tasks[!file.exists(tasks$file), ]
cat(sprintf("%s: %d of %d replications done, %d to run.\n",
            format(Sys.time()), nrow(tasks) - nrow(todo), nrow(tasks), nrow(todo)))

run_task <- function(k) {
  tk <- todo[k, ]
  if (file.exists(file.path(RAW_DIR, "STOP"))) return("stopped")
  if (file.exists(tk$file)) return("exists")
  t0 <- proc.time()[3]
  out <- tryCatch(run_one(tk$case, tk$design, tk$scenario, tk$rep,
                          methods = METHODS,
                          keep_curves = tk$design %in% CURVE_DESIGNS),
                  error = function(e) e)
  if (inherits(out, "error")) {
    cat(sprintf("%s ERROR %s rep %d: %s\n", format(Sys.time()), tk$cell, tk$rep,
                conditionMessage(out)),
        file = file.path(RAW_DIR, "errors.log"), append = TRUE)
    return("error")
  }
  tmp <- paste0(tk$file, ".tmp", Sys.getpid())
  saveRDS(out, tmp)
  file.rename(tmp, tk$file)
  cat(sprintf("%s %s rep %d done in %.0f s\n", format(Sys.time()), tk$cell,
              tk$rep, proc.time()[3] - t0),
      file = file.path(RAW_DIR, "progress.log"), append = TRUE)
  "ok"
}

## Run the missing replications in parallel, then combine all results.
if (nrow(todo) > 0) {
  library(parallel)
  cl <- makePSOCKcluster(min(N_WORKERS, nrow(todo)))
  status <- tryCatch({
    clusterExport(cl, c("todo", "RAW_DIR", "METHODS", "CURVE_DESIGNS", "PKG_DIR",
                        "SIM_DIR", "load_sim_code", "run_task"))
    invisible(clusterEvalQ(cl, {
      Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
                 MKL_NUM_THREADS = "1")
      load_sim_code(PKG_DIR, SIM_DIR)
    }))
    unlist(parLapplyLB(cl, seq_len(nrow(todo)), function(k) run_task(k),
                       chunk.size = 1))
  }, finally = stopCluster(cl))
  print(table(status))
}

source(file.path(SIM_DIR, "collect_results.R"))
collect_results(RAW_DIR, file.path(SIM_DIR, "results"))
