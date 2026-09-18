# R2.7: the joint-regression baseline the reviewer asks for -- platform entered as a
# covariate in the model, with NO batch-corrected matrix produced first.
#
# Two arms, scored identically:
#   limma_2step : removeBatchEffect to build a corrected matrix, then the Wilcoxon
#                 rank-sum test on the pooled result  (the manuscript's limma arm)
#   limma_1step : lmFit on the pooled data with design ~ group + platform, eBayes,
#                 p-value and logFC taken from the group coefficient -- no corrected
#                 matrix at any point
#
# Both are scored with the SAME rules as the manuscript's scorer, reimplemented here
# because that scorer takes corrected matrices and runs Wilcoxon internally, which the
# one-step arm by definition does not produce:
#   calls      : p <= 0.05, split by the sign of the fold change
#   power      : correctly-directed calls among the 2*n_DE true DE genes
#   Type I     : false positives among the true null genes / number of true nulls
# Type I error is therefore only defined on the typeI grid, which has n_DE = 0.
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
source("/work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/code/sim_rerun_common.R")
suppressPackageStartupMessages({library(limma)})

grp_of_cn <- function(cn) ifelse(grepl("control", cn, ignore.case = TRUE), 0, 1)

# shared scoring, applied to (p, lfc) from either arm
score_pl <- function(p, lfc, genes, n_DE, alpha = 0.05) {
  G <- length(p)
  true_up   <- if (n_DE > 0) genes[1:n_DE] else character(0)
  true_down <- if (n_DE > 0) genes[(n_DE + 1):(2 * n_DE)] else character(0)
  true_de   <- c(true_up, true_down)
  true_null <- setdiff(genes, true_de)
  pred_up   <- genes[lfc > 0 & p <= alpha]
  pred_down <- genes[lfc < 0 & p <= alpha]
  pred_all  <- c(pred_up, pred_down)
  TP <- sum(pred_up %in% true_up) + sum(pred_down %in% true_down)
  FP <- sum(pred_all %in% true_null)
  c(typeI = FP / length(true_null),
    power = if (length(true_de)) TP / length(true_de) else NA_real_)
}

# ---- arm 1: correct, then Wilcoxon (the manuscript's limma) ------------------
arm_2step <- function(a, s, n_DE) {
  r  <- METHODS$limma(a, s)
  da <- as.matrix(r$array); ds <- as.matrix(r$seq)
  dimnames(da) <- dimnames(a); dimnames(ds) <- dimnames(s)
  d  <- cbind(da, ds); y <- grp_of_cn(colnames(d))
  ctl <- d[, y == 0, drop = FALSE]; cas <- d[, y == 1, drop = FALSE]
  p <- lfc <- numeric(nrow(d))
  for (i in seq_len(nrow(d))) {
    p[i]   <- suppressWarnings(wilcox.test(ctl[i, ], cas[i, ])$p.value)
    lfc[i] <- mean(cas[i, ]) - mean(ctl[i, ])
  }
  score_pl(p, lfc, rownames(d), n_DE)
}

# ---- arm 2: one model, platform as a covariate, no corrected matrix ----------
arm_1step <- function(a, s, n_DE) {
  q <- qn_within(a, s)                       # same preprocessing as the limma arm
  d <- as.matrix(cbind(q$array, q$seq))
  dimnames(d) <- list(rownames(a), c(colnames(a), colnames(s)))
  grp  <- factor(ifelse(grp_of_cn(colnames(d)) == 1, "case", "control"),
                 levels = c("control", "case"))
  plat <- factor(c(rep("array", ncol(a)), rep("seq", ncol(s))))
  des  <- model.matrix(~ grp + plat)
  fit  <- eBayes(lmFit(d, des))
  tt   <- topTable(fit, coef = "grpcase", number = Inf, sort.by = "none")
  score_pl(tt$P.Value, tt$logFC, rownames(tt), n_DE)
}

seed <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "101"))
out  <- file.path(ROOT, "revision_repo/results_limma_1step")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
NDE <- list(typeI_balance = 0, power_size = 1000, power_balance = 1000, power_effect = 1000)
rows <- list(); t0 <- Sys.time()
for (gname in names(GRIDS)) {
  g <- GRIDS[[gname]]
  for (lv in g$levels) {
    r <- g$gen(lv, seed)
    a <- as.matrix(r$sim_array); s <- as.matrix(r$sim_seq)
    for (arm in c("limma_2step", "limma_1step")) {
      v <- tryCatch(if (arm == "limma_2step") arm_2step(a, s, NDE[[gname]])
                    else                      arm_1step(a, s, NDE[[gname]]),
                    error = function(e) { cat("  ERR", arm, gname, lv, conditionMessage(e), "\n")
                                          c(typeI = NA_real_, power = NA_real_) })
      rows[[length(rows) + 1]] <- data.frame(grid = gname, level = as.character(lv),
        seed = seed, method = arm, Type_I_error = v[["typeI"]], Power = v[["power"]])
      write.csv(do.call(rbind, rows), file.path(out, sprintf("l1s_seed%d.csv", seed)), row.names = FALSE)
      cat(sprintf("[%5.1f min] %-14s %-15s %-12s tI=%-8s pw=%-8s\n",
          as.numeric(difftime(Sys.time(), t0, units = "mins")), gname, lv, arm,
          format(v[["typeI"]], digits = 3), format(v[["power"]], digits = 3))); flush.console()
    }
  }
}
cat("\nDONE seed", seed, "\n")
