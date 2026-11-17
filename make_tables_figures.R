## Tables and figures of the simulation study, built from
## results/sim_metrics_all.rds and results/sim_curves.rds (see
## collect_results.R). No model is refitted, so the script also works on a
## partial set of replications. Run it from the simulation folder; the output
## goes to results/tables and results/figures.

N_REP <- 250   # replications 1, ..., N_REP of every cell are used

suppressPackageStartupMessages(library(ggplot2))

res_dir <- "results"
tab_dir <- file.path(res_dir, "tables")
fig_dir <- file.path(res_dir, "figures")
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

m <- readRDS(file.path(res_dir, "sim_metrics_all.rds"))
curves <- readRDS(file.path(res_dir, "sim_curves.rds"))
m <- m[m$rep <= N_REP, ]
curves <- lapply(curves, function(cc) {
  keep <- cc$reps <= N_REP
  cc$reps <- cc$reps[keep]
  cc$curves <- cc$curves[keep]
  cc
})

method_order <- c("LPFR", "LFPCR", "MAD-LFPCR", "Bag-LFPCR", "KFPCA-LS", "rLMM", "RLFPCR")
scen_order <- c("clean", "mixed05", "mixed10", "spike", "drift", "cluster",
                "covariate", "batch20", "batch30", "instrument")
scen_label <- c(clean = "Clean", mixed05 = "Mixed 5\\%", mixed10 = "Mixed 10\\%",
                spike = "Local spikes", drift = "Subject drift", cluster = "Clustered visits",
                covariate = "Covariate outliers", batch20 = "Batch 20\\%",
                batch30 = "Batch 30\\%", instrument = "Instrument drift")
metrics <- c("mse", "mse_delta", "ise", "ise_proj", "ise_B", "ise_B_proj", "ise_U", "ise_U_proj")
m$method <- factor(m$method, levels = intersect(method_order, unique(m$method)))
m$scenario <- factor(m$scenario, levels = intersect(scen_order, unique(m$scenario)))
m$design <- paste0("I=", m$nI, ", m=", m$nV)
designs <- unique(m[order(m$nI, m$nV), c("nI", "nV", "design")])

## Only replications in which every method returned a fit are used.
m$ok <- is.finite(m$mse)
okrep <- aggregate(ok ~ case + nI + nV + scenario + rep, data = m, FUN = all)
names(okrep)[6] <- "all_ok"
m <- merge(m, okrep, by = c("case", "nI", "nV", "scenario", "rep"))
mc <- m[m$all_ok, ]

## Summaries (CSV)

summ <- do.call(rbind, lapply(split(mc, list(mc$case, mc$nI, mc$nV, mc$scenario, mc$method),
                                    drop = TRUE), function(x) {
  data.frame(case = x$case[1], nI = x$nI[1], nV = x$nV[1], scenario = x$scenario[1],
             method = x$method[1], n = nrow(x),
             t(sapply(metrics, function(v) mean(x[[v]]))),
             t(setNames(sapply(metrics, function(v) sd(x[[v]])), paste0(metrics, "_sd"))),
             time = median(x$time), check.names = FALSE)
}))
summ <- summ[order(summ$case, summ$nI, summ$nV, summ$scenario, summ$method), ]
write.csv(summ, file.path(tab_dir, "summary_all.csv"), row.names = FALSE)

fails <- aggregate(cbind(failed = !ok) ~ case + nI + nV + scenario + method, data = m, FUN = sum)
write.csv(fails[fails$failed > 0, ], file.path(tab_dir, "failures.csv"), row.names = FALSE)
nrep <- aggregate(all_ok ~ case + nI + nV + scenario, data = okrep, FUN = function(z) c(total = length(z), common = sum(z)))
write.csv(do.call(data.frame, nrep), file.path(tab_dir, "replications.csv"), row.names = FALSE)

## LaTeX helpers

fmt <- function(x, d = 3) ifelse(is.na(x), "--", formatC(x, digits = d, format = "f"))
cell_ms <- function(mu, s) ifelse(is.na(mu), "--", paste0(fmt(mu), " (", fmt(s), ")"))
get_val <- function(cs, nI, nV, sc, me, v) {
  r <- summ[summ$case == cs & summ$nI == nI & summ$nV == nV &
              summ$scenario == sc & summ$method == me, ]
  if (nrow(r) == 0) return(c(NA, NA))
  c(r[[v]], r[[paste0(v, "_sd")]])
}
header_block <- function(scens, vlab) {
  nv <- length(vlab)
  top <- paste0(" & ", paste0("\\multicolumn{", nv, "}{c}{", scen_label[scens], "}",
                               collapse = " & "), " \\\\")
  rules <- paste0(sapply(seq_along(scens), function(k) {
    a <- 2 + (k - 1) * nv
    sprintf("\\cmidrule(lr){%d-%d}", a, a + nv - 1)
  }), collapse = " ")
  sub <- paste0("Method & ", paste(rep(vlab, length(scens)), collapse = " & "), " \\\\")
  c(top, rules, sub)
}
write_table <- function(lines, file) writeLines(lines, file.path(tab_dir, file))

## Tables 1 and 2: rows are designs and methods, columns are settings and
## errors; entries are mean (SD).
main_table <- function(cs, scens, vars, vlab, caption, label, file) {
  meths <- levels(droplevels(mc$method[mc$case == cs]))
  L <- c("\\begin{sidewaystable}[p]", "\\centering",
         paste0("\\caption{", caption, "}"), paste0("\\label{", label, "}"),
         "\\scriptsize", "\\setlength{\\tabcolsep}{2pt}",
         paste0("\\begin{tabular}{l", strrep("c", length(scens) * length(vars)), "}"),
         "\\toprule", header_block(scens, vlab), "\\midrule")
  for (k in seq_len(nrow(designs))) {
    dz <- designs[k, ]
    L <- c(L, sprintf("\\multicolumn{%d}{l}{\\textit{$I=%d$, $m=%d$}} \\\\",
                      1 + length(scens) * length(vars), dz$nI, dz$nV))
    for (me in meths) {
      cells <- unlist(lapply(scens, function(sc) sapply(vars, function(v) {
        z <- get_val(cs, dz$nI, dz$nV, sc, me, v)
        cell_ms(z[1], z[2])
      })))
      L <- c(L, paste0(me, " & ", paste(cells, collapse = " & "), " \\\\"))
    }
    if (k < nrow(designs)) L <- c(L, "\\midrule")
  }
  L <- c(L, "\\bottomrule", "\\end{tabular}", "\\end{sidewaystable}")
  write_table(L, file)
}

## Tables S1-S4, the seven additional settings: rows are settings and
## methods, columns are designs and errors; entries are means.
new_table <- function(cs, scens, vars, vlab, caption, label, file) {
  meths <- levels(droplevels(mc$method[mc$case == cs]))
  nv <- length(vars)
  top <- paste0(" & ", paste0("\\multicolumn{", nv, "}{c}{$I=", designs$nI, "$, $m=",
                               designs$nV, "$}", collapse = " & "), " \\\\")
  rules <- paste0(sapply(seq_len(nrow(designs)), function(k) {
    a <- 2 + (k - 1) * nv
    sprintf("\\cmidrule(lr){%d-%d}", a, a + nv - 1)
  }), collapse = " ")
  sub <- paste0("Method & ", paste(rep(vlab, nrow(designs)), collapse = " & "), " \\\\")
  L <- c("\\begin{sidewaystable}[p]", "\\centering",
         paste0("\\caption{", caption, "}"), paste0("\\label{", label, "}"),
         "\\scriptsize", "\\setlength{\\tabcolsep}{3pt}",
         paste0("\\begin{tabular}{l", strrep("c", nrow(designs) * nv), "}"),
         "\\toprule", top, rules, sub, "\\midrule")
  for (j in seq_along(scens)) {
    sc <- scens[j]
    L <- c(L, sprintf("\\multicolumn{%d}{l}{\\textit{%s}} \\\\", 1 + nrow(designs) * nv,
                      scen_label[sc]))
    for (me in meths) {
      cells <- unlist(lapply(seq_len(nrow(designs)), function(k) sapply(vars, function(v) {
        fmt(get_val(cs, designs$nI[k], designs$nV[k], sc, me, v)[1])
      })))
      L <- c(L, paste0(me, " & ", paste(cells, collapse = " & "), " \\\\"))
    }
    if (j < length(scens)) L <- c(L, "\\midrule")
  }
  L <- c(L, "\\bottomrule", "\\end{tabular}", "\\end{sidewaystable}")
  write_table(L, file)
}

orig <- intersect(c("clean", "mixed05", "mixed10"), levels(mc$scenario))
newsc <- intersect(scen_order[-(1:3)], levels(mc$scenario))
mlab <- c("MSE", "MSE$_\\delta$")
n_note <- paste0("Entries are Monte Carlo means with standard deviations in parentheses over ",
                 "replications in which all methods returned a fit. ")

if (any(mc$case == 1)) {
  main_table(1, orig, c("mse", "mse_delta", "ise"), c(mlab, "ISE"),
             paste0("Monte Carlo results for Case~1. ", n_note,
                    "ISE is computed against the coefficient that generates the response, ",
                    "$1.40\\beta$, after Gram-corrected reconstruction."),
             "tab:case1", "tab_case1_main.tex")
  if (length(newsc) > 0) {
    half <- split(newsc, ceiling(seq_along(newsc) / 4))
    for (h in seq_along(half)) {
      new_table(1, half[[h]], c("mse", "mse_delta", "ise"), c(mlab, "ISE"),
                paste0("Case~1 under the additional contamination settings (part ", h,
                       "). Entries are Monte Carlo means."),
                paste0("tab:case1new", h), paste0("tab_case1_new_part", h, ".tex"))
    }
  }
}
if (any(mc$case == 2)) {
  main_table(2, orig, c("mse", "mse_delta", "ise_B", "ise_U"),
             c(mlab, "ISE$_B$", "ISE$_U$"),
             paste0("Monte Carlo results for Case~2. ", n_note,
                    "ISE$_B$ and ISE$_U$ are computed against $1.35\\beta_B$ and $1.10\\beta_U$. ",
                    "LPFR, KFPCA-LS and rLMM estimate a single coefficient and have no componentwise ISE."),
             "tab:case2", "tab_case2_main.tex")
  if (length(newsc) > 0) {
    half <- split(newsc, ceiling(seq_along(newsc) / 4))
    for (h in seq_along(half)) {
      new_table(2, half[[h]], c("mse", "mse_delta", "ise_B", "ise_U"),
                c(mlab, "ISE$_B$", "ISE$_U$"),
                paste0("Case~2 under the additional contamination settings (part ", h,
                       "). Entries are Monte Carlo means."),
                paste0("tab:case2new", h), paste0("tab_case2_new_part", h, ".tex"))
    }
  }
}

## Table 3 (RLFPCR relative to LFPCR) and Table S6 (best competitor)

paired <- function(cs, v) {
  a <- mc[mc$case == cs & mc$method == "RLFPCR", c("nI", "nV", "scenario", "rep", v)]
  b <- mc[mc$case == cs & mc$method == "LFPCR", c("nI", "nV", "scenario", "rep", v)]
  names(a)[5] <- "r"
  names(b)[5] <- "l"
  merge(a, b, by = c("nI", "nV", "scenario", "rep"))
}
ratio_rows <- function(cs, vars, vlab) {
  scens <- levels(droplevels(mc$scenario[mc$case == cs]))
  out <- NULL
  for (sc in scens) {
    vals <- c()
    for (k in seq_len(nrow(designs))) for (v in vars) {
      p <- paired(cs, v)
      p <- p[p$scenario == sc & p$nI == designs$nI[k] & p$nV == designs$nV[k], ]
      vals <- c(vals, if (nrow(p) == 0) "--" else
        paste0(fmt(mean(p$r) / mean(p$l), 2), " (", round(100 * mean(p$r < p$l)), ")"))
    }
    out <- c(out, paste0(scen_label[sc], " & ", paste(vals, collapse = " & "), " \\\\"))
  }
  out
}
ratio_table <- function() {
  L <- c("\\begin{sidewaystable}[p]", "\\centering",
         "\\caption{RLFPCR relative to LFPCR in all contamination settings. Each entry is the ratio of Monte Carlo means (RLFPCR/LFPCR); values below one favour RLFPCR. The percentage of replications in which RLFPCR has the smaller error is given in parentheses.}",
         "\\label{tab:ratio}", "\\scriptsize", "\\setlength{\\tabcolsep}{3pt}")
  for (cs in intersect(1:2, unique(mc$case))) {
    vars <- if (cs == 1) c("mse", "ise") else c("mse", "ise_B", "ise_U")
    vlab <- if (cs == 1) c("MSE", "ISE") else c("MSE", "ISE$_B$", "ISE$_U$")
    nv <- length(vars)
    top <- paste0(" & ", paste0("\\multicolumn{", nv, "}{c}{$I=", designs$nI, "$, $m=",
                                 designs$nV, "$}", collapse = " & "), " \\\\")
    rules <- paste0(sapply(seq_len(nrow(designs)), function(k) {
      a <- 2 + (k - 1) * nv
      sprintf("\\cmidrule(lr){%d-%d}", a, a + nv - 1)
    }), collapse = " ")
    L <- c(L, if (cs == 2) "\\vspace{0.3cm}" else NULL,
           paste0("\\begin{tabular}{l", strrep("c", nrow(designs) * nv), "}"),
           "\\toprule", sprintf("\\multicolumn{%d}{l}{\\textbf{Case %d}} \\\\", 1 + nv * nrow(designs), cs),
           top, rules,
           paste0("Setting & ", paste(rep(vlab, nrow(designs)), collapse = " & "), " \\\\"),
           "\\midrule", ratio_rows(cs, vars, vlab), "\\bottomrule", "\\end{tabular}")
    if (cs == 1) L <- c(L, "")
  }
  write_table(c(L, "\\end{sidewaystable}"), "tab_ratio_rlfpcr_lfpcr.tex")
}
ratio_table()

best_table <- function() {
  L <- c("\\begin{table}[!htbp]", "\\centering",
         "\\caption{Prediction MSE of RLFPCR divided by the smallest MSE among the six competitors; the competitor attaining it is given in parentheses. Values below one mean that RLFPCR is the most accurate method.}",
         "\\label{tab:best}", "\\scriptsize", "\\setlength{\\tabcolsep}{3pt}",
         paste0("\\begin{tabular}{ll", strrep("c", nrow(designs)), "}"), "\\toprule",
         paste0("Case & Setting & ", paste0("$I=", designs$nI, "$, $m=", designs$nV, "$",
                                             collapse = " & "), " \\\\"), "\\midrule")
  for (cs in intersect(1:2, unique(summ$case))) {
    for (sc in levels(droplevels(summ$scenario[summ$case == cs]))) {
      vals <- sapply(seq_len(nrow(designs)), function(k) {
        s <- summ[summ$case == cs & summ$scenario == sc & summ$nI == designs$nI[k] &
                    summ$nV == designs$nV[k], ]
        if (!"RLFPCR" %in% s$method) return("--")
        oth <- s[s$method != "RLFPCR", ]
        b <- oth[which.min(oth$mse), ]
        paste0(fmt(s$mse[s$method == "RLFPCR"] / b$mse, 2), " (", b$method, ")")
      })
      L <- c(L, paste0(cs, " & ", scen_label[sc], " & ", paste(vals, collapse = " & "), " \\\\"))
    }
    if (cs == 1 && 2 %in% summ$case) L <- c(L, "\\midrule")
  }
  write_table(c(L, "\\bottomrule", "\\end{tabular}", "\\end{table}"), "tab_best_competitor.tex")
}
best_table()

## Table S5: ISE against the true coefficient and against its projection onto
## the span retained by each method (clean and mixed settings).
proj_table <- function() {
  L <- c("\\begin{table}[!htbp]", "\\centering",
         "\\caption{Coefficient errors against two targets: the generating coefficient (ISE) and its projection onto the span retained by each method (ISE$^{\\mathrm{proj}}$). For the score-based methods the span is that of the estimated eigenfunctions; for rLMM it is the space of bin-wise constant functions; LPFR has no finite retained span.}",
         "\\label{tab:proj}", "\\scriptsize", "\\setlength{\\tabcolsep}{3pt}")
  meth_proj <- c("LFPCR", "MAD-LFPCR", "Bag-LFPCR", "KFPCA-LS", "rLMM", "RLFPCR")
  for (cs in intersect(1:2, unique(summ$case))) {
    vars <- if (cs == 1) c("ise", "ise_proj") else c("ise_B", "ise_B_proj", "ise_U", "ise_U_proj")
    vlab <- if (cs == 1) c("ISE", "ISE$^{\\mathrm{proj}}$") else
      c("ISE$_B$", "ISE$_B^{\\mathrm{proj}}$", "ISE$_U$", "ISE$_U^{\\mathrm{proj}}$")
    mp <- if (cs == 1) meth_proj else c("LFPCR", "MAD-LFPCR", "Bag-LFPCR", "RLFPCR")
    L <- c(L, paste0("\\begin{tabular}{lll", strrep("c", length(vars)), "}"), "\\toprule",
           sprintf("\\multicolumn{%d}{l}{\\textbf{Case %d}} \\\\", 3 + length(vars), cs),
           paste0("Design & Setting & Method & ", paste(vlab, collapse = " & "), " \\\\"), "\\midrule")
    for (k in seq_len(nrow(designs))) for (sc in orig) for (me in mp) {
      s <- summ[summ$case == cs & summ$scenario == sc & summ$nI == designs$nI[k] &
                  summ$nV == designs$nV[k] & summ$method == me, ]
      if (nrow(s) == 0) next
      L <- c(L, paste0("$I=", designs$nI[k], "$, $m=", designs$nV[k], "$ & ", scen_label[sc],
                       " & ", me, " & ", paste(fmt(unlist(s[vars])), collapse = " & "), " \\\\"))
    }
    L <- c(L, "\\bottomrule", "\\end{tabular}", "")
  }
  write_table(c(L, "\\end{table}"), "tab_ise_projection.tex")
}
proj_table()

## Screening counts and run times (CSV)
scr <- aggregate(cbind(n_curve_drop, n_resp_drop) ~ case + nI + nV + scenario + method,
                 data = mc[mc$method %in% c("MAD-LFPCR", "Bag-LFPCR", "RLFPCR"), ], FUN = mean)
write.csv(scr, file.path(tab_dir, "screening_counts.csv"), row.names = FALSE)
rt <- aggregate(time ~ nI + nV + method, data = mc, FUN = median)
write.csv(rt, file.path(tab_dir, "runtime_median_seconds.csv"), row.names = FALSE)

## Figures 1 and 2, and the log ratios of the prediction errors

gray_col <- rgb(0.7, 0.7, 0.7, 0.45)
curve_panel <- function(mat, truth, tps, main, ylim, ylab = NULL) {
  plot(tps, truth, type = "n", ylim = ylim, xlab = "", ylab = "", main = main,
       bty = "l", cex.main = 1.4, cex.axis = 1.2)
  for (i in seq_len(nrow(mat))) lines(tps, mat[i, ], col = gray_col, lwd = 1)
  lines(tps, truth, col = "black", lwd = 3)
  lines(tps, colMeans(mat), col = "blue", lwd = 3)
  if (!is.null(ylab)) mtext(ylab, side = 2, line = 2.3, cex = 1.2)
}
curve_matrix <- function(cc, me, comp) {
  rows <- lapply(cc$curves, function(z) z[[me]][[comp]])
  rows <- rows[!sapply(rows, is.null)]
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

for (cn in names(curves)) {
  cc <- curves[[cn]]
  if (!grepl("_(clean|mixed10)$", cn)) next
  cs <- as.integer(sub("^case([0-9]).*", "\\1", cn))
  if (cs == 1) {
    meths <- intersect(c("LPFR", "LFPCR", "MAD-LFPCR", "KFPCA-LS", "RLFPCR"),
                       names(cc$curves[[1]]))
    mats <- lapply(meths, curve_matrix, cc = cc, comp = "beta")
    yl <- range(cc$truth$beta_true, unlist(lapply(mats, function(x) apply(x, 2, quantile, c(0.01, 0.99)))))
    pdf(file.path(fig_dir, paste0("fig_", cn, ".pdf")), width = 3 * length(meths), height = 3.2)
    par(mfrow = c(1, length(meths)), mar = c(2.5, 3.5, 2, 0.5))
    for (j in seq_along(meths)) curve_panel(mats[[j]], cc$truth$beta_true, cc$tps, meths[j], yl,
                                            if (j == 1) expression(beta(s)) else NULL)
    dev.off()
  } else {
    meths <- intersect(c("LFPCR", "MAD-LFPCR", "Bag-LFPCR", "RLFPCR"), names(cc$curves[[1]]))
    pdf(file.path(fig_dir, paste0("fig_", cn, ".pdf")), width = 3 * length(meths), height = 6)
    par(mfrow = c(2, length(meths)), mar = c(2.5, 3.5, 2, 0.5))
    for (comp in c("beta_B", "beta_U")) {
      tr <- if (comp == "beta_B") cc$truth$beta_B_true else cc$truth$beta_U_true
      yname <- if (comp == "beta_B") expression(beta[B](s)) else expression(beta[U](s))
      mats <- lapply(meths, curve_matrix, cc = cc, comp = comp)
      yl <- range(tr, unlist(lapply(mats, function(x) apply(x, 2, quantile, c(0.01, 0.99)))))
      for (j in seq_along(meths)) curve_panel(mats[[j]], tr, cc$tps, meths[j], yl,
                                              if (j == 1) yname else NULL)
    }
    dev.off()
  }
}

## Paired log ratio of prediction MSE, RLFPCR versus LFPCR, in every setting.
pr <- do.call(rbind, lapply(intersect(1:2, unique(mc$case)), function(cs) {
  p <- paired(cs, "mse")
  p$case <- paste("Case", cs)
  p
}))
if (!is.null(pr) && nrow(pr) > 0) {
  pr$design <- factor(paste0("I=", pr$nI, ", m=", pr$nV),
                      levels = paste0("I=", designs$nI, ", m=", designs$nV))
  pr$setting <- factor(gsub("\\\\", "", scen_label[as.character(pr$scenario)]),
                       levels = rev(gsub("\\\\", "", scen_label[levels(mc$scenario)])))
  g <- ggplot(pr, aes(x = log(r / l), y = setting)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    geom_boxplot(outlier.size = 0.4, fill = "grey85") +
    facet_grid(case ~ design) +
    labs(x = "log(MSE of RLFPCR / MSE of LFPCR)", y = NULL) +
    theme_bw(base_size = 10)
  ggsave(file.path(fig_dir, "fig_logratio_mse_rlfpcr_lfpcr.pdf"), g, width = 10, height = 6)
}

cat("Tables written to", tab_dir, "\nFigures written to", fig_dir, "\n")
