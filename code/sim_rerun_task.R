# One SLURM array task = one seed. Usage:
#   Rscript sim_rerun_task.R <tag> <method[,method...]>
# with the seed taken from SLURM_ARRAY_TASK_ID (101..200).
args <- commandArgs(trailingOnly = TRUE)
tag  <- args[1]; methods <- strsplit(args[2], ",")[[1]]
seed <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "101"))

source("/work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/code/sim_rerun_common.R")
suppressPackageStartupMessages({
  library(MatchMixeR)                                   # provides MM(), xpn() + MatchMixeR.so
  if (any(methods == "MNN"))        library(batchelor)
  if (any(methods == "ComBat_seq")) library(sva)
})

dir.create(file.path(ROOT, "revision_repo/results"), showWarnings = FALSE, recursive = TRUE)
out <- file.path(ROOT, sprintf("revision_repo/results/%s_seed%d.csv", tag, seed))
cat("tag:", tag, " seed:", seed, " methods:", paste(methods, collapse=","), "\n")
cat("out:", out, "\n\n")
run_seed(seed, methods, out)
cat("\nDONE", tag, "seed", seed, "\n")
