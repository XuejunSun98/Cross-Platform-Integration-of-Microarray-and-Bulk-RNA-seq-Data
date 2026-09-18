# Reassemble the 11 Shambhala2 chunks (job array 1212298) and run the same
# cross-platform prediction as metsim_prediction.R, so the Shambhala2 row can sit
# alongside the other methods in the expanded Table 2.
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
suppressPackageStartupMessages({library(glmnet); library(pROC)})
ROOT <- "/work/users/x/u/xuejun1/Integration_paper_Sim"
load(file.path(ROOT, "real_data_METSIM/METSIM_analysis.rda"))
a <- as.matrix(d_array); s <- as.matrix(d_seq)
y <- ifelse(grepl("^case", colnames(a)), 1, 0)

f <- sort(list.files(file.path(ROOT, "revision_repo/results_metsim_sbl"),
                     "^sbl_chunk[0-9]+\\.rds$", full.names = TRUE))
stopifnot(length(f) == 11)
ch <- lapply(f, readRDS)

# rebuild in ORIGINAL column order; refuse to guess if anything is missing or misaligned
idx <- unlist(lapply(ch, `[[`, "idx"))
stopifnot(!anyDuplicated(idx), setequal(idx, seq_len(ncol(s))))
genes <- Reduce(intersect, lapply(ch, function(z) rownames(z$ds)))
ds <- matrix(NA_real_, nrow = length(genes), ncol = ncol(s),
             dimnames = list(genes, colnames(s)))
for (z in ch) {
  stopifnot(identical(colnames(z$ds), NULL) || ncol(z$ds) == length(z$idx))
  ds[, z$idx] <- z$ds[genes, , drop = FALSE]
}
stopifnot(!anyNA(ds))
# no transform: Shambhala2's output is already on the array's log2 scale
cat(sprintf("reassembled: %d genes x %d samples from %d chunks\n", nrow(ds), ncol(ds), length(ch)))
cat(sprintf("array range %.2f-%.2f | corrected seq range %.2f-%.2f\n",
            min(a), max(a), min(ds), max(ds)))

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

g  <- intersect(rownames(a), rownames(ds))
ca <- a[g, , drop = FALSE]; cs <- ds[g, , drop = FALSE]
ok <- complete.cases(ca) & complete.cases(cs)
ca <- ca[ok, , drop = FALSE]; cs <- cs[ok, , drop = FALSE]
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
