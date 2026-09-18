# Shambhala2 on METSIM, CHUNKED.
#
# Measured rate on the single-process run (job 1210276): 12 of 330 samples in 40.7 min
# = 3.39 min/sample, i.e. ~18.7 h for all 330 -- over the 12 h wall limit it was given.
#
# Shambhala2's `Input` argument is the set of samples to harmonise (NH = ncol(MAS)-1);
# `P` is the calibration pool it is merged against and `Q` the target shape. Splitting
# Input by columns while keeping P0 and Q0 complete therefore leaves the calibration
# identical -- each chunk's samples are harmonised against exactly the same pool as in
# the single-process run.
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
source("/work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/code/sim_rerun_common.R")
source(file.path(ROOT, "revision_repo/code/longleaf_external_tools.R"))
# upstream Shambhala2.R has library(matrixStats) on line 1, but grab() extracts ONLY the
# function definition, so that call is left behind and rowSds() is missing at the very
# last step -- after all 30 samples have been harmonised. Load it explicitly.
suppressPackageStartupMessages(library(matrixStats))
src <- grab(readLines(file.path(ROOT, "Shambhala2_upstream/Shambhala2.R"), warn = FALSE),
            "^Shambhala2 <- function")
src <- sub('system\\("matlab[^\n]*\\)', 'shambhala2_matlab()', src)
stopifnot(grepl("shambhala2_matlab\\(\\)", src))
eval(parse(text = src), envir = globalenv())

load(file.path(ROOT, "real_data_METSIM/METSIM_analysis.rda"))
a <- as.matrix(d_array); s <- as.matrix(d_seq)   # RAW: Shambhala2 maps onto Q0's (log2 array) scale

NCHUNK <- 11
ci  <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
idx <- split(seq_len(ncol(s)), cut(seq_len(ncol(s)), NCHUNK, labels = FALSE))[[ci]]
cat(sprintf("chunk %d/%d: %d samples (%d..%d)\n", ci, NCHUNK, length(idx), min(idx), max(idx)))

wd <- file.path(Sys.getenv("TMPDIR", unset = tempdir()), sprintf("sbl_metsim_%d", ci))
dir.create(wd, recursive = TRUE, showWarnings = FALSE)
old <- setwd(wd); on.exit(setwd(old))

# Q0 = microarray (target shape); P0 = FULL seq (calibration pool); Input = this chunk
write.csv(data.frame(SYMBOL = rownames(a), a, check.names = FALSE), "Q0.csv", row.names = FALSE)
write.csv(data.frame(SYMBOL = rownames(s), s, check.names = FALSE), "P0.csv", row.names = FALSE)
write.csv(data.frame(SYMBOL = rownames(s), s[, idx, drop = FALSE], check.names = FALSE),
          "Input.csv", row.names = FALSE)

t0 <- Sys.time()
H  <- Shambhala2("Input.csv", "P0.csv", "Q0.csv", delete_buffer_files = TRUE, k = 5)
ds <- as.matrix(apply(H[, -1, drop = FALSE], 2, as.numeric))
rownames(ds) <- as.character(H[, 1])
cat(sprintf("chunk %d done in %.1f min, dim %d x %d\n", ci,
    as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(ds), ncol(ds)))

out <- file.path(ROOT, "revision_repo/results_metsim_sbl")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(chunk = ci, idx = idx, cols = colnames(s)[idx], ds = ds),
        file.path(out, sprintf("sbl_chunk%02d.rds", ci)))
