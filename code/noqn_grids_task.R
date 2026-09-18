# R2.5: does the within-platform quantile normalisation, rather than the batch-correction
# algorithm, drive the Type I error and power results for ComBat and limma?
#
# Runs the manuscript's OWN four grids and OWN scorers, so Type I error keeps its usual
# definition (the typeI_balance grid has n_DE = 0) and power is the usual power. The only
# thing that changes between arms is whether qn_within() is applied before the correction.
#
# Output goes to results_noqn_grids/, NOT results/ -- sim_rerun_plots.R globs every CSV in
# results/, so writing there would silently add these arms to Figures 1 and 3.
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
source("/work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/code/sim_rerun_common.R")
suppressPackageStartupMessages({library(sva); library(limma)})

# identical to m_ComBat / m_limma except qn_within() -> log(s+1) only. The log is kept
# because both are additive models and would otherwise fit a log2-scale array against raw
# counts; it is the QUANTILE NORMALISATION that is being removed, not the transformation.
noqn <- function(a, s) list(array = as.matrix(a), seq = log(as.matrix(s) + 1))
METHODS$ComBat_noQN <- function(a, s) {
  q <- noqn(a, s); d <- as.matrix(cbind(q$array, q$seq))
  b <- factor(c(rep("array", ncol(q$array)), rep("seq", ncol(q$seq))))
  out <- sva::ComBat(dat = d, batch = b, mod = model.matrix(~ 1, data.frame(x = rep(1, ncol(d)))))
  list(array = out[, seq_len(ncol(q$array)), drop = FALSE],
       seq   = out[, ncol(q$array) + seq_len(ncol(q$seq)), drop = FALSE])
}
METHODS$limma_noQN <- function(a, s) {
  q <- noqn(a, s); d <- as.matrix(cbind(q$array, q$seq))
  b <- factor(c(rep("array", ncol(q$array)), rep("seq", ncol(q$seq))))
  out <- limma::removeBatchEffect(d, batch = b)
  list(array = out[, seq_len(ncol(q$array)), drop = FALSE],
       seq   = out[, ncol(q$array) + seq_len(ncol(q$seq)), drop = FALSE])
}

seed <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "101"))
out  <- file.path(ROOT, "revision_repo/results_noqn_grids")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
cat("seed:", seed, "\n")
run_seed(seed, c("ComBat", "ComBat_noQN", "limma", "limma_noQN"),
         file.path(out, sprintf("noqn_seed%d.csv", seed)))
cat("\nDONE seed", seed, "\n")
