# Shambhala2 for the two power grids it was never completed on:
#   power_balance  400/400 missing in the shipped d_plot
#   power_effect   280/400 missing (only 30 of 100 seeds per level)
# It is slow -- which is why those arms were abandoned -- so one array task
# handles ONE (grid, level, seed) combination. 2 grids x 4 levels x 100 seeds = 800.
#
# Each task runs in its OWN directory: Shambhala2 writes fixed filenames
# (Input.csv, P0.csv, Q0.csv, P_prim.txt, args.txt, Cu_bis.txt), so concurrent
# tasks sharing a directory would overwrite each other's inputs.
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
source("/work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/code/sim_rerun_common.R")
source(file.path(ROOT, "longleaf_external_tools.R"))

# upstream Shambhala2.R defines the wrapper AND calls it at the bottom -- take
# the function only, then swap its hardcoded MATLAB invocation
#   system("matlab -nodesktop -nosplash -nodisplay -r \"run('Shambhala2.m');exit;\"")
# for shambhala2_matlab(), which copies Shambhala2.m / CuBlock.m /
# readExpressionData.m into the working directory first and calls
# `matlab -batch`. Without the copy the run dies with "cannot open the
# connection", because the .m files are not in the per-task directory.
src <- grab(readLines(file.path(ROOT, "Shambhala2_upstream/Shambhala2.R"), warn = FALSE),
            "^Shambhala2 <- function")
src <- sub('system\\("matlab[^\n]*\\)', 'shambhala2_matlab()', src)
stopifnot(grepl("shambhala2_matlab\\(\\)", src))
# grab() takes only the function body, dropping upstream's line-1 library(matrixStats);
# rowSds() is then missing at the very last step, AFTER every sample is harmonised.
suppressPackageStartupMessages(library(matrixStats))
eval(parse(text = src), envir = globalenv())

i     <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
combos <- expand.grid(seed = 101:200, level = 1:4, grid = c("power_balance", "power_effect"),
                      stringsAsFactors = FALSE)
stopifnot(i >= 1, i <= nrow(combos))
cb   <- combos[i, ]
gname <- cb$grid
lv    <- GRIDS[[gname]]$levels[cb$level]
seed  <- cb$seed

wd <- file.path(Sys.getenv("TMPDIR", unset = tempdir()), sprintf("sbl_%s_%s_%d", gname, lv, seed))
dir.create(wd, recursive = TRUE, showWarnings = FALSE)
old <- setwd(wd); on.exit(setwd(old))
cat("task", i, ":", gname, lv, "seed", seed, "\n  wd:", wd, "\n"); flush.console()

r  <- GRIDS[[gname]]$gen(lv, seed)
# Scale handling follows the original scripts (unbalanced_power_sample_size_sbl.R and
# friends): the RNA-seq is passed RAW, and Shambhala2's output is used AS IS. Shambhala2
# quantile-maps Input onto Q0's distribution, and Q0 is the log2-scale microarray, so the
# harmonized matrix is already on the array scale. Logging the input and the output (as an
# earlier version of this script did) compresses the result and destroys the comparison.
a  <- as.matrix(r$sim_array); s <- as.matrix(r$sim_seq)
# Q = microarray (the definitive dataset that sets the target shape);
# P = Input = the RNA-seq. P = Input is a known deviation, kept because a
# two-platform comparison provides no third dataset -- see RUNLIST.
write.csv(data.frame(SYMBOL = rownames(r$sim_array), a, check.names = FALSE), "Q0.csv", row.names = FALSE)
write.csv(data.frame(SYMBOL = rownames(r$sim_seq),   s, check.names = FALSE), "P0.csv", row.names = FALSE)
file.copy("P0.csv", "Input.csv", overwrite = TRUE)

t0  <- Sys.time()
res <- tryCatch({
  H  <- Shambhala2("Input.csv", "P0.csv", "Q0.csv", delete_buffer_files = TRUE, k = 5)
  ds <- as.matrix(apply(H[, -1, drop = FALSE], 2, as.numeric))
  rownames(ds) <- as.character(H[, 1])
  ds <- ds[rownames(r$sim_seq), , drop = FALSE]   # already on the array's scale
  v  <- score_pair(list(array = a, seq = ds), dimnames(r$sim_array), dimnames(r$sim_seq), gname)
  data.frame(typeI = v[["typeI"]], power = v[["power"]], status = "OK")
}, error = function(e) data.frame(typeI = NA_real_, power = NA_real_,
                                  status = substr(conditionMessage(e), 1, 80)))

out <- file.path(ROOT, "revision_repo/results_sbl")
dir.create(out, showWarnings = FALSE)
write.csv(data.frame(grid = gname, level = as.character(lv), seed = seed, method = "Shambhala2",
                     Type_I_error = res$typeI, Power = res$power, status = res$status,
                     mins = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)),
          file.path(out, sprintf("sbl_%s_%s_%d.csv", gname, lv, seed)), row.names = FALSE)
cat(sprintf("  %s  power=%s  %.1f min\n", res$status, format(res$power, digits = 4),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
unlink(wd, recursive = TRUE)
