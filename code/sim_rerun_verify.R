# ---------------------------------------------------------------
# VERIFICATION GATE -- run this before submitting anything.
#
# The stage-1 .rda files are gone, so every grid is regenerated in-script. That
# is only safe if the regenerated data is identical to what produced the shipped
# results. Proven for power_size (the MMR pilot matched published values
# exactly); NOT yet proven for typeI_balance, power_balance or power_effect,
# which use different generators.
#
# Test: run a PUBLISHED method (QN -- fast, deterministic, in every d_plot) on
# regenerated data and compare with the shipped value for the same (level, seed).
# ---------------------------------------------------------------
source("/work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/code/sim_rerun_common.R")
suppressPackageStartupMessages(library(preprocessCore))
setwd(SIM)

qn_pair <- function(a, s) {          # QN as the method blocks do it
  tgt <- as.matrix(a)[, 1]
  list(array = normalize.quantiles.use.target(as.matrix(a), target = tgt),
       seq   = normalize.quantiles.use.target(as.matrix(s), target = tgt))
}
shipped <- function(f) { e <- new.env(); load(f, envir = e); e$d_plot }

checks <- list(
  list(grid="power_size",    lv=50,         file="d_plot_Power_Sample_Size.rda", col="Power",
       key=function(d,lv,sd) d$sample_size==lv & d$seed==sd),
  list(grid="power_balance", lv="same_imbalance", file="d_plot_Power_Imbalance.rda", col="Power",
       key=function(d,lv,sd) d$balance==lv & d$seed==sd),
  list(grid="power_effect",  lv=0.5,        file="d_plot_Power_Effect_Size.rda", col="Power",
       key=function(d,lv,sd) d$effect_size==lv & d$seed==sd),
  list(grid="typeI_balance", lv="same_imbalance", file="d_plot_Type_I_error_small.rda", col="Type_I_error",
       key=function(d,lv,sd) d$balance==lv & d$seed==sd))

cat(sprintf("%-16s %-16s %10s %10s %10s  %s\n","grid","level","regenerated","shipped","diff","verdict"))
for (ck in checks) {
  ok <- TRUE
  for (sd in c(101, 102)) {
    r <- GRIDS[[ck$grid]]$gen(ck$lv, sd)
    v <- score_pair(qn_pair(r$sim_array, r$sim_seq),
                    dimnames(r$sim_array), dimnames(r$sim_seq), ck$grid)
    mine <- if (ck$col == "Power") v[["power"]] else v[["typeI"]]
    d    <- shipped(ck$file)
    pubv <- d[ck$key(d, ck$lv, sd) & d$method == "QN", ck$col]
    pubv <- if (length(pubv)) pubv[1] else NA_real_
    dif  <- abs(mine - pubv)
    good <- isTRUE(dif < 1e-8); ok <- ok && good
    cat(sprintf("%-16s %-16s %10.5f %10.5f %10.2e  %s\n", ck$grid,
                paste0(ck$lv," s",sd), mine, pubv, dif,
                if (good) "MATCH" else "*** MISMATCH ***"))
  }
  if (!ok) cat("   -> regeneration does NOT reproduce this grid; do not submit it.\n")
}
cat("\nAll four must MATCH before submitting. A mismatch means the generator or its\n",
    "arguments differ from what produced the shipped results.\n")
