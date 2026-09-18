# GATE: at n_up = n_down = 1000 the generalised generator and scorer must be
# IDENTICAL to the published ones. Scenario 2000_1:1 is the published
# power-by-sample-size point at m=50, so QN must reproduce the shipped d_plot
# value exactly. Run this before submitting anything.
source("/work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/code/de_sweep_common.R")
suppressPackageStartupMessages(library(preprocessCore))
qn_pair <- function(a, s) { tgt <- as.matrix(a)[,1]
  list(array = normalize.quantiles.use.target(as.matrix(a), target = tgt),
       seq   = normalize.quantiles.use.target(as.matrix(s), target = tgt)) }
e <- new.env(); load(file.path(SIM, "d_plot_Power_Sample_Size.rda"), envir = e); pub <- e$d_plot

cat(sprintf("%-8s %12s %12s %10s  %s\n", "seed", "regenerated", "shipped", "diff", "verdict"))
ok_all <- TRUE
for (sd in c(101, 102, 103)) {
  # 1. identical DATA: new generator at n_up=n_down=1000 vs published at n_DE=1000
  r1 <- sim_array_seq_deratio(n_per_group = 50, n_up = 1000, n_down = 1000, seed = sd)
  r2 <- sim_array_seq_balanced(n_per_group = 50, n_DE = 1000, seed = sd)
  same <- identical(as.matrix(r1$sim_array), as.matrix(r2$sim_array)) &&
          identical(as.matrix(r1$sim_seq),   as.matrix(r2$sim_seq))
  # 2. identical SCORE against the shipped value
  v <- score_sweep(qn_pair(r1$sim_array, r1$sim_seq),
                   dimnames(r1$sim_array), dimnames(r1$sim_seq), "2000_1:1")
  p <- pub$Power[pub$method == "QN" & pub$sample_size == 50 & pub$seed == sd]
  d <- abs(v[["power"]] - p); good <- isTRUE(d < 1e-12) && same
  ok_all <- ok_all && good
  cat(sprintf("%-8d %12.5f %12.5f %10.2e  %s%s\n", sd, v[["power"]], p, d,
              if (good) "MATCH" else "*** MISMATCH ***",
              if (!same) "  (DATA DIFFERS)" else ""))
}
cat(if (ok_all) "\nGATE PASSED - generalisation changes nothing at 1:1.\n"
    else "\nGATE FAILED - do not submit.\n")
