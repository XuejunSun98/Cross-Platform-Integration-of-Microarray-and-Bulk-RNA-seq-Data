# Does the Type I inflation of ComBat and limma under directional imbalance come
# from the methods themselves, or from the within-platform quantile normalisation
# they are wrapped in?
#
# The diagnostic said QN: ComBat and limma produce IDENTICAL null-gene shifts
# (0.0155 / -0.0890 array, 0.0000 / -0.0734 seq at 1:1 / 1:0), which two different
# batch-removal algorithms cannot do by chance, and neither batch step can change a
# standardised within-batch shift (ComBat's is a per-gene affine map applied to all
# samples in a batch; removeBatchEffect subtracts a per-gene batch coefficient).
#
# This script tests that directly: same ratio sweep, same scorer, four arms --
# ComBat and limma as published (with qn_within) and with the QN step replaced by a
# plain log(s+1) so the two platforms are merely put on a comparable scale.
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
source("/work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/code/de_ratio_task_defs.R")

# no-QN counterparts: identical except qn_within() -> log(s+1) only
noqn <- function(a, s) list(array = as.matrix(a), seq = log(as.matrix(s) + 1))
m_ComBat_noQN <- function(a, s) {
  q <- noqn(a, s); d <- as.matrix(cbind(q$array, q$seq))
  b <- factor(c(rep("array", ncol(q$array)), rep("seq", ncol(q$seq))))
  out <- sva::ComBat(dat = d, batch = b, mod = model.matrix(~ 1, data.frame(x = rep(1, ncol(d)))))
  list(array = out[, seq_len(ncol(q$array)), drop = FALSE],
       seq   = out[, ncol(q$array) + seq_len(ncol(q$seq)), drop = FALSE])
}
m_limma_noQN <- function(a, s) {
  q <- noqn(a, s); d <- as.matrix(cbind(q$array, q$seq))
  b <- factor(c(rep("array", ncol(q$array)), rep("seq", ncol(q$seq))))
  out <- limma::removeBatchEffect(d, batch = b)
  list(array = out[, seq_len(ncol(q$array)), drop = FALSE],
       seq   = out[, ncol(q$array) + seq_len(ncol(q$seq)), drop = FALSE])
}
ARMS <- list(ComBat = SWEEP_METHODS$ComBat, ComBat_noQN = m_ComBat_noQN,
             limma  = SWEEP_METHODS$limma,  limma_noQN  = m_limma_noQN)

i   <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "101"))
out <- file.path(ROOT, "revision_repo/results_noqn")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
rows <- list(); t0 <- Sys.time()
for (rt in RATIOS) {
  k <- mk(rt)
  r <- sim_array_seq_deratio(n_per_group = 50, n_up = k$n_up, n_down = k$n_down, seed = i)
  a <- as.matrix(r$sim_array); s <- as.matrix(r$sim_seq)
  for (mn in names(ARMS)) {
    v <- tryCatch(score_sd(ARMS[[mn]](a, s), dimnames(r$sim_array),
                           dimnames(r$sim_seq), k$n_up, k$n_down),
                  error = function(e) { cat("  ", mn, rt, "ERR:", conditionMessage(e), "\n")
                                        rep(NA_real_, 6) })
    rows[[length(rows)+1]] <- data.frame(ratio = rt, n_up = k$n_up, n_down = k$n_down,
      seed = i, method = mn, t(v))
    write.csv(do.call(rbind, rows), file.path(out, sprintf("noqn_seed%d.csv", i)), row.names = FALSE)
    cat(sprintf("[%5.1f min] %-4s %-12s tI=%.4f pw=%.4f AUC=%.4f FDR=%.4f\n",
        as.numeric(difftime(Sys.time(), t0, units = "mins")), rt, mn,
        v[["typeI_raw"]], v[["power_raw"]], v[["AUC"]], v[["FDR_BH"]])); flush.console()
  }
}
cat("\nDONE seed", i, "\n")
