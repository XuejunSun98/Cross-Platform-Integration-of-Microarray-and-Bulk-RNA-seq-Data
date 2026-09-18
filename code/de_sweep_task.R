# One SLURM array task = one seed, all 12 DE-sweep scenarios.
#   Rscript de_sweep_task.R <tag> <method[,method...]>
args <- commandArgs(trailingOnly = TRUE)
tag  <- args[1]; methods <- strsplit(args[2], ",")[[1]]
seed <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "101"))
source("/work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/code/de_sweep_common.R")
suppressPackageStartupMessages({
  library(MatchMixeR); library(limma)
  if (any(methods=="MNN")) library(batchelor)
  if (any(methods %in% c("ComBat","ComBat_seq","RNABC"))) library(sva)
  if (any(methods=="TDM")) library(TDM)})

out <- file.path(ROOT, "revision_repo/results_de_sweep")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
f <- file.path(out, sprintf("%s_seed%d.csv", tag, seed))
rows <- list(); t0 <- Sys.time()
for (lv in GRIDS[["de_sweep"]]$levels) {
  r <- GRIDS[["de_sweep"]]$gen(lv, seed)
  a <- as.matrix(r$sim_array); s <- as.matrix(r$sim_seq)
  dn_a <- dimnames(r$sim_array); dn_s <- dimnames(r$sim_seq)
  k <- SWEEP[SWEEP$level == lv, ]
  for (mn in methods) {
    t1 <- Sys.time()
    v <- tryCatch(score_sweep(SWEEP_METHODS[[mn]](a, s), dn_a, dn_s, lv),
                  error = function(e) { cat("  ", mn, lv, "ERR:", conditionMessage(e), "\n")
                    c(typeI=NA_real_, power=NA_real_, power_up=NA_real_, power_down=NA_real_) })
    rows[[length(rows)+1]] <- data.frame(level = lv, total = k$total, ratio = k$ratio,
      n_up = k$n_up, n_down = k$n_down, seed = seed, method = mn,
      Type_I_error = v[["typeI"]], Power = v[["power"]],
      Power_up = v[["power_up"]], Power_down = v[["power_down"]],
      secs = round(as.numeric(difftime(Sys.time(), t1, units="secs")), 1))
    write.csv(do.call(rbind, rows), f, row.names = FALSE)
    cat(sprintf("[%5.1f min] %-10s %-11s tI=%-8s pw=%-8s up=%-7s dn=%-7s %5.1fs\n",
        as.numeric(difftime(Sys.time(), t0, units="mins")), lv, mn,
        format(v[["typeI"]],digits=3), format(v[["power"]],digits=3),
        format(v[["power_up"]],digits=3), format(v[["power_down"]],digits=3),
        rows[[length(rows)]]$secs)); flush.console()
  }
}
cat("\nDONE", tag, "seed", seed, "\n")
