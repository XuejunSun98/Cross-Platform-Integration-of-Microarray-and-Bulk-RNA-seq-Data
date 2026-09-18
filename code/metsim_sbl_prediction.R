# Shambhala2 on METSIM, then the same cross-platform prediction as metsim_prediction.R.
# Kept separate because Shambhala2 shells out to MATLAB and needs its own working
# directory (it writes fixed filenames: Input.csv / P0.csv / Q0.csv -> Cu_bis.txt).
#
# Orientation follows the validated sbl_rerun_task.R:
#   Q0 = microarray (the definitive dataset that sets the target shape)
#   P0 = Input     = the RNA-seq.  P = Input is a known deviation from the published
#                    3-dataset design, kept because a two-platform comparison offers
#                    no third dataset (see RUNLIST).
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
source("/work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/code/sim_rerun_common.R")
source(file.path(ROOT, "revision_repo/code/longleaf_external_tools.R"))
suppressPackageStartupMessages({library(glmnet); library(pROC)})

# upstream Shambhala2.R defines the wrapper AND calls it at the bottom; take only the
# function, and swap its hardcoded MATLAB system() call for shambhala2_matlab().
src <- grab(readLines(file.path(ROOT, "Shambhala2_upstream/Shambhala2.R"), warn = FALSE),
            "^Shambhala2 <- function")
src <- sub('system\\("matlab[^\n]*\\)', 'shambhala2_matlab()', src)
stopifnot(grepl("shambhala2_matlab\\(\\)", src))
eval(parse(text = src), envir = globalenv())

load(file.path(ROOT, "real_data_METSIM/METSIM_analysis.rda"))
a <- as.matrix(d_array); s <- as.matrix(d_seq)
y <- ifelse(grepl("^case", colnames(a)), 1, 0)

wd <- file.path(Sys.getenv("TMPDIR", unset = tempdir()), "sbl_metsim")
dir.create(wd, recursive = TRUE, showWarnings = FALSE)
old <- setwd(wd); on.exit(setwd(old))
cat("wd:", wd, "\n"); flush.console()

write.csv(data.frame(SYMBOL = rownames(a), a, check.names = FALSE), "Q0.csv", row.names = FALSE)
write.csv(data.frame(SYMBOL = rownames(s), log(s + 1), check.names = FALSE), "P0.csv", row.names = FALSE)
file.copy("P0.csv", "Input.csv", overwrite = TRUE)

t0 <- Sys.time()
H  <- Shambhala2("Input.csv", "P0.csv", "Q0.csv", delete_buffer_files = TRUE, k = 5)
ds <- as.matrix(apply(H[, -1, drop = FALSE], 2, as.numeric))
rownames(ds) <- as.character(H[, 1])
ds <- log(ds + 1)[rownames(s), , drop = FALSE]
colnames(ds) <- colnames(s)
cat(sprintf("Shambhala2 done in %.1f min | array range %.2f-%.2f  corrected seq range %.2f-%.2f\n",
    as.numeric(difftime(Sys.time(), t0, units = "mins")),
    min(a), max(a), min(ds, na.rm = TRUE), max(ds, na.rm = TRUE))); flush.console()

fit_eval <- function(xtr, ytr, xte, yte, alpha) {
  keep <- apply(xtr, 2, sd, na.rm = TRUE) > 0
  xtr <- xtr[, keep, drop = FALSE]; xte <- xte[, keep, drop = FALSE]
  cv  <- cv.glmnet(xtr, ytr, family = "binomial", alpha = alpha, nfolds = 10)
  fit <- glmnet(xtr, ytr, family = "binomial", alpha = alpha, lambda = cv$lambda.1se)
  pp  <- as.numeric(predict(fit, newx = xte, s = cv$lambda.1se, type = "response"))
  c(AUC_prob  = as.numeric(auc(roc(yte, pp, quiet = TRUE))),
    AUC_class = as.numeric(auc(roc(yte, as.numeric(pp > 0.5), quiet = TRUE))),
    nonzero   = sum(as.matrix(coef(fit))[-1, 1] != 0))
}

ok <- complete.cases(a) & complete.cases(ds)
ca <- a[ok, , drop = FALSE]; cs <- ds[ok, , drop = FALSE]
set.seed(98); half <- sample(ncol(ca), ncol(ca) %/% 2)
ALPHA <- c(lasso = 1, elastic_net = 0.5, ridge = 0); out <- list()
for (dir in c("array_to_seq", "seq_to_array")) {
  tr <- if (dir == "array_to_seq") ca else cs
  te <- if (dir == "array_to_seq") cs else ca
  for (pn in names(ALPHA)) {
    pub <- fit_eval(t(tr), y, t(te), y, ALPHA[[pn]])
    spl <- fit_eval(t(tr[, half]), y[half], t(te[, -half]), y[-half], ALPHA[[pn]])
    out[[length(out) + 1]] <- data.frame(method = "Shambhala2", direction = dir, penalty = pn,
      paired_AUC_prob = pub[["AUC_prob"]], paired_AUC_class = pub[["AUC_class"]],
      split_AUC_prob = spl[["AUC_prob"]], split_AUC_class = spl[["AUC_class"]],
      nonzero = pub[["nonzero"]], genes = nrow(ca))
  }
}
res <- do.call(rbind, out)
write.csv(res, file.path(ROOT, "revision_repo/results/metsim_prediction_shambhala2.csv"), row.names = FALSE)
print(res, row.names = FALSE)
