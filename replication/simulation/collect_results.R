## Combines the replication files of run_simulations.R into
##   results/sim_metrics_all.rds  one row per replication and method
##   results/sim_curves.rds       estimated coefficient curves (CURVE_DESIGNS)
## It also works while the simulation is still running, e.g.
##   source("collect_results.R"); collect_results("raw", "results")
## check_progress() prints the number of finished replications in each cell.

collect_results <- function(raw_dir, out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  cells <- list.dirs(raw_dir, full.names = TRUE, recursive = FALSE)
  metrics <- list()
  curves <- list()
  for (cd in cells) {
    files <- list.files(cd, pattern = "^rep_[0-9]+[.]rds$", full.names = TRUE)
    if (length(files) == 0) next
    res <- lapply(files, readRDS)
    metrics[[basename(cd)]] <- do.call(rbind, lapply(res, `[[`, "metrics"))
    if (length(res[[1]]$curves) > 0) {
      reps <- as.integer(sub("^rep_([0-9]+)[.]rds$", "\\1", basename(files)))
      curves[[basename(cd)]] <- list(reps = reps,
                                     curves = lapply(res, `[[`, "curves"),
                                     truth = res[[1]]$truth, tps = res[[1]]$tps)
    }
  }
  m <- do.call(rbind, metrics)
  rownames(m) <- NULL
  saveRDS(m, file.path(out_dir, "sim_metrics_all.rds"))
  saveRDS(curves, file.path(out_dir, "sim_curves.rds"))
  done <- aggregate(rep ~ case + nI + nV + scenario, data = unique(m[, c("case", "nI", "nV", "scenario", "rep")]),
                    FUN = length)
  names(done)[5] <- "replications"
  cat(sprintf("Collected %d replications in %d cells into %s\n",
              sum(done$replications), nrow(done), out_dir))
  invisible(done)
}

check_progress <- function(raw_dir = "raw") {
  cells <- list.dirs(raw_dir, full.names = FALSE, recursive = FALSE)
  n <- sapply(cells, function(cd) length(list.files(file.path(raw_dir, cd),
                                                    pattern = "^rep_[0-9]+[.]rds$")))
  print(data.frame(cell = cells, done = n, row.names = NULL))
  cat("Total:", sum(n), "\n")
  invisible(n)
}
