## LaTeX tables of the sensitivity study (Tables S7-S10 of Additional file 1).
## run_sensitivity.R calls make_sensitivity_tables() at the end. To rebuild
## the tables from the saved summary, run from this folder
##   source("make_sensitivity_tables.R")
##   make_sensitivity_tables(read.csv("results/sensitivity_summary.csv"),
##                           "results")

make_sensitivity_tables <- function(agg, out_dir) {
  f3 <- function(x) ifelse(is.na(x), "--", formatC(x, format = "f", digits = 3))
  pick <- function(z, method, scen, col) {
    v <- z[z$method == method & z$scenario == scen, col]
    if (length(v) == 0) NA else v
  }

  ## Part A (Tables S7 and S8): one table per case, rows are the variants
  ## within each design.
  methA <- c("LFPCR", "c = 3.5", "c = 4.685", "c = 6", "cells 3.5 MAD",
             "cells 6 MAD", "no cell step")
  labA <- c("LFPCR", "$c=3.5$", "$c=4.685$ (default)", "$c=6$",
            "cells at 3.5 MAD", "cells at 6 MAD", "no cellwise step")
  scen <- c("clean", "mixed10", "spike")
  for (cs in 1:2) {
    zA <- agg[agg$part == "A" & agg$case == cs, ]
    metrics <- if (cs == 1) c("mse", "ise") else c("mse", "ise_B", "ise_U")
    mlab <- c(mse = "MSE", ise = "ISE", ise_B = "ISE$_B$", ise_U = "ISE$_U$")
    ncol <- 3 * length(metrics)
    head1 <- paste(sprintf("\\multicolumn{3}{c}{%s}", mlab[metrics]), collapse = " & ")
    rules <- paste(sprintf("\\cmidrule(lr){%d-%d}", 2 + 3 * (seq_along(metrics) - 1),
                           4 + 3 * (seq_along(metrics) - 1)), collapse = "")
    head2 <- paste(rep("Clean & Mixed & Spike", length(metrics)), collapse = " & ")
    lines <- c(sprintf("\\begin{tabular}{l%s}", strrep("r", ncol)), "\\toprule",
               paste0(" & ", head1, " \\\\"), rules,
               paste0("Method & ", head2, " \\\\"), "\\midrule")
    for (nI in sort(unique(zA$nI))) {
      z <- zA[zA$nI == nI, ]
      lines <- c(lines, sprintf("\\multicolumn{%d}{l}{\\emph{$I=%d$, $m=4$}} \\\\",
                                ncol + 1, nI))
      for (k in seq_along(methA)) {
        vals <- unlist(lapply(metrics, function(mt)
          sapply(scen, function(sc) f3(pick(z, methA[k], sc, mt)))))
        lines <- c(lines, paste0(labA[k], " & ", paste(vals, collapse = " & "), " \\\\"))
      }
    }
    lines <- c(lines, "\\bottomrule", "\\end{tabular}")
    writeLines(lines, file.path(out_dir, sprintf("tab_tuning_case%d.tex", cs)))
  }

  ## Part B (Tables S9 and S10): one table per case, rows are the
  ## contamination proportions.
  methB <- c("LFPCR", "Bag-LFPCR", "RLFPCR")
  for (cs in 1:2) {
    zB <- agg[agg$part == "B" & agg$case == cs, ]
    metrics <- if (cs == 1) c("mse", "ise") else c("mse", "ise_B", "ise_U")
    mlab <- c(mse = "MSE", ise = "ISE", ise_B = "ISE$_B$", ise_U = "ISE$_U$")
    ncol <- 3 + 3 * length(metrics)
    head1 <- paste(sprintf("\\multicolumn{3}{c}{%s}", mlab[metrics]), collapse = " & ")
    rules <- paste(sprintf("\\cmidrule(lr){%d-%d}", 4 + 3 * (seq_along(metrics) - 1),
                           6 + 3 * (seq_along(metrics) - 1)), collapse = "")
    head2 <- paste(rep("LFPCR & Bag & RLFPCR", length(metrics)), collapse = " & ")
    lines <- c(sprintf("\\begin{tabular}{rrr%s}", strrep("r", 3 * length(metrics))),
               "\\toprule",
               paste0("& \\multicolumn{2}{c}{Contaminated (\\%)} & ", head1, " \\\\"),
               paste0("\\cmidrule(lr){2-3}", rules),
               paste0("$p$ & Curves & Responses & ", head2, " \\\\"), "\\midrule")
    for (pp in sort(unique(zB$p))) {
      z <- zB[zB$p == pp, ]
      fc <- 100 * z$frac_curves[z$method == "RLFPCR"]
      fr <- 100 * z$frac_resp[z$method == "RLFPCR"]
      if (pp == 0) {
        fc <- 0
        fr <- 0
      }
      vals <- unlist(lapply(metrics, function(mt)
        sapply(methB, function(m) f3(z[z$method == m, mt]))))
      lines <- c(lines, paste0(formatC(pp, format = "f", digits = 2), " & ",
                               formatC(fc, format = "f", digits = 0), " & ",
                               formatC(fr, format = "f", digits = 0), " & ",
                               paste(vals, collapse = " & "), " \\\\"))
    }
    lines <- c(lines, "\\bottomrule", "\\end{tabular}")
    writeLines(lines, file.path(out_dir, sprintf("tab_breakdown_case%d.tex", cs)))
  }
  invisible(TRUE)
}
