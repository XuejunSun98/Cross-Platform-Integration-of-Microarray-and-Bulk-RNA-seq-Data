# Up:down ratio sweep at FIXED total DE = 2000, evaluated as in
# Soneson & Delorenzi 2013 (Genome Biology), the benchmark compcodeR implements.
#
# Following that paper:
#  * the SAME threshold is applied to every method -- no per-method calibration;
#  * significance is BH-adjusted p < 0.05 (an FDR threshold), and we report the
#    TRUE FDR realised at it, alongside TPR;
#  * AUC is reported as the threshold-free metric -- it is what exposed the
#    degradation under all-upregulated DE in their Figure 1C vs 1D;
#  * the paper's own raw-p Type I error and power are kept alongside so the
#    numbers stay comparable with the rest of this manuscript.
#
# Total DE is fixed at 2000 so the ratio is the only factor varying. The extreme
# arm 1:0 (all up) matches their B01250 / B04000 scenarios.
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
source("/work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/code/de_sweep_common.R")
suppressPackageStartupMessages({library(MatchMixeR); library(limma); library(sva)
  library(preprocessCore); library(TDM); library(batchelor); library(data.table); library(pROC)})

TOTAL  <- 2000
RATIOS <- c("1:1","2:1","3:1","4:1","9:1","1:0")      # 1:0 = all up (extreme)
mk <- function(rt) {
  if (rt == "1:0") return(list(n_up = TOTAL, n_down = 0))
  r <- as.numeric(strsplit(rt, ":")[[1]]); u <- round(TOTAL * r[1]/sum(r))
  list(n_up = u, n_down = TOTAL - u)
}

gene_stats <- function(da, ds) {                       # per-gene Wilcoxon + log FC
  d <- cbind(da, ds); cn <- colnames(d)
  g <- ifelse(grepl("control", cn, ignore.case = TRUE), 0, 1)
  ctl <- as.matrix(d[, g == 0, drop = FALSE]); cas <- as.matrix(d[, g == 1, drop = FALSE])
  p <- lfc <- numeric(nrow(d))
  for (i in seq_len(nrow(d))) {
    xc <- ctl[i, ]; xk <- cas[i, ]
    # direction is taken from the MEAN DIFFERENCE, not the log ratio. MatchMixeR
    # returns fitted values that can be negative (min -0.075 here), so log(mean)
    # is undefined for 945 genes -- 147 of them truly DE. The paper's scorer
    # silently drops those; the sign of the mean difference agrees with the sign
    # of the log ratio on 100% of genes where the latter is defined, and is
    # well-defined everywhere.
    lfc[i] <- mean(xk) - mean(xc)
    p[i]   <- suppressWarnings(wilcox.test(xc, xk)$p.value)
  }
  list(p = p, lfc = lfc)
}

score_sd <- function(res, dn_a, dn_s, n_up, n_down, alpha = 0.05) {
  da <- as.data.frame(res$array); ds <- as.data.frame(res$seq)
  dimnames(da) <- dn_a; dimnames(ds) <- dn_s
  o <- gene_stats(da, ds); G <- length(o$p)
  up  <- if (n_up   > 0) 1:n_up else integer(0)
  dn  <- if (n_down > 0) (n_up + 1):(n_up + n_down) else integer(0)   # empty when 1:0
  nul <- (n_up + n_down + 1):G
  is_de <- rep(FALSE, G); is_de[c(up, dn)] <- TRUE
  padj  <- p.adjust(o$p, method = "BH")

  called  <- function(th, pv) pv <= th
  correct <- function(sel) sum(sel[up] & o$lfc[up] > 0) + sum(sel[dn] & o$lfc[dn] < 0)
  # --- S&D: BH-adjusted threshold, same for every method ---
  sel <- called(alpha, padj)
  n_called <- sum(sel)
  c(TPR_BH   = if (length(c(up,dn))) correct(sel)/(n_up + n_down) else NA_real_,
    FDR_BH   = if (n_called > 0) sum(sel[nul])/n_called else 0,
    n_called = n_called,
    AUC      = as.numeric(pROC::auc(pROC::roc(is_de, -o$p, quiet = TRUE, direction = "<"))),
    # --- this manuscript's own convention, kept for comparability ---
    typeI_raw = sum(o$p[nul] <= alpha)/length(nul),
    power_raw = correct(called(alpha, o$p))/(n_up + n_down))
}

i    <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "101"))
meth <- strsplit(Sys.getenv("METHODS",
         unset = "QN,Angel,TDM,MMR,ComBat,ComBat_seq,RNABC,limma,MNN"), ",")[[1]]
# namespace the output by array, otherwise the cheap and xpn arrays write the
# SAME filename for a given seed and the second to finish silently overwrites
# the first
TAG <- sub("^[a-z]+_", "", Sys.getenv("SLURM_JOB_NAME", unset = "manual"))
out  <- file.path(ROOT, "revision_repo/results_de_ratio")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
rows <- list(); t0 <- Sys.time()
for (rt in RATIOS) {
  k <- mk(rt)
  r <- sim_array_seq_deratio(n_per_group = 50, n_up = k$n_up, n_down = k$n_down, seed = i)
  a <- as.matrix(r$sim_array); s <- as.matrix(r$sim_seq)
  for (mn in meth) {
    v <- tryCatch(score_sd(SWEEP_METHODS[[mn]](a, s), dimnames(r$sim_array),
                           dimnames(r$sim_seq), k$n_up, k$n_down),
                  error = function(e) { cat("  ", mn, rt, "ERR:", conditionMessage(e), "\n")
                                        rep(NA_real_, 6) })
    rows[[length(rows)+1]] <- data.frame(ratio = rt, n_up = k$n_up, n_down = k$n_down,
      seed = i, method = mn, t(v))
    write.csv(do.call(rbind, rows), file.path(out, sprintf("%s_seed%d.csv", TAG, i)), row.names = FALSE)
    cat(sprintf("[%5.1f min] %-4s %-11s AUC=%.4f TPR=%.4f FDR=%.4f | tI=%.4f pw=%.4f\n",
        as.numeric(difftime(Sys.time(), t0, units = "mins")), rt, mn,
        v[["AUC"]], v[["TPR_BH"]], v[["FDR_BH"]], v[["typeI_raw"]], v[["power_raw"]]))
    flush.console()
  }
}
cat("\nDONE seed", i, "\n")
