# METSIM cross-platform prediction -- expansion of manuscript Table 2 to three penalties.
#
# Mirrors the original Real_data_prediction script: set.seed(98), cv.glmnet with
# nfolds = 10, family = "binomial", refit at lambda.1se, train on one platform and
# test on the other, 330 paired subjects, outcome = Matsuda Index (R = insulin
# resistant = case if <= 4, N = control otherwise).
#
# TWO DEVIATIONS FROM THE ORIGINAL ARE RECORDED RATHER THAN SILENTLY FIXED, and both
# are reported side by side so the published column can be reproduced and compared:
#
#  (1) AUC definition. The original computes AUC for 8 of the 9 methods from the
#      PREDICTED CLASS LABELS (predict(type="class") recoded to 0/1), which makes AUC
#      equal to balanced accuracy, (sens+spec)/2. Only the limma block uses predicted
#      PROBABILITIES. We report AUC_class (published convention) and AUC_prob (the
#      standard definition) for every method, so limma is no longer scored on a
#      different metric from its competitors.
#  (2) Subject overlap. The original trains on all 330 subjects' array profiles and
#      tests on THE SAME 330 subjects' RNA-seq profiles, so every test subject's label
#      was used in training. We report the published "paired" design and a "split"
#      design in which subjects are partitioned 50/50, training on one platform in the
#      training subjects and testing on the other platform in the held-out subjects.
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
source("/work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/code/sim_rerun_common.R")
suppressPackageStartupMessages({
  library(glmnet); library(pROC); library(preprocessCore); library(limma)
  library(sva); library(MatchMixeR); library(TDM); library(data.table)})
load("/work/users/x/u/xuejun1/Integration_paper_Sim/real_data_METSIM/METSIM_analysis.rda")

a <- as.matrix(d_array); s <- as.matrix(d_seq)
y <- ifelse(grepl("^case", colnames(a)), 1, 0)          # case = Matsuda < 4 = "R"
cat(sprintf("METSIM: %d genes x %d paired subjects | cases %d controls %d\n",
            nrow(a), ncol(a), sum(y == 1), sum(y == 0)))

M <- list(
  No_correction = function(a,s) list(array=a, seq=log2(s+1)),
  QN = function(a,s) { tgt <- as.matrix(a)[,1]
    list(array=normalize.quantiles.use.target(as.matrix(a), target=tgt),
         seq  =normalize.quantiles.use.target(log2(as.matrix(s)+1), target=tgt)) },
  TDM = function(a,s) { r <- tdm_transform(
      ref_data    = data.table(cbind(gene=rownames(a), as.data.frame(a))),
      target_data = data.table(cbind(gene=rownames(s), as.data.frame(round(s)))))
    m <- as.matrix(r[,-1]); rownames(m) <- r$gene; list(array=a, seq=m[rownames(a),,drop=FALSE]) },
  Angel = function(a,s) { d <- cbind(a, log2(s+1))
    d <- apply(d, 2, function(x) rank(x)/length(x))
    list(array=d[,1:ncol(a),drop=FALSE], seq=d[,ncol(a)+seq_len(ncol(s)),drop=FALSE]) },
  MMR = METHODS$MMR, ComBat = METHODS$ComBat, RNABC = METHODS$RNABC, limma = METHODS$limma,
  # added in this revision (R1.2 / R3.m1): XPN, MNN and ComBat-seq
  XPN = METHODS$XPN, MNN = METHODS$MNN, ComBat_seq = METHODS$ComBat_seq)
# COCONUT and Rank-In are excluded throughout the prediction analysis because they use
# the outcome during correction. ComBat-seq expects integer counts; METSIM RNA-seq is
# distributed as TPM, so the TPM matrix is rounded to serve as pseudo-counts -- an
# approximation that is disclosed rather than hidden.

ALPHA <- c(lasso = 1, elastic_net = 0.5, ridge = 0)

fit_eval <- function(xtr, ytr, xte, yte, alpha) {
  keep <- apply(xtr, 2, sd, na.rm = TRUE) > 0
  xtr <- xtr[, keep, drop=FALSE]; xte <- xte[, keep, drop=FALSE]
  cv  <- cv.glmnet(xtr, ytr, family="binomial", alpha=alpha, nfolds=10)
  fit <- glmnet(xtr, ytr, family="binomial", alpha=alpha, lambda=cv$lambda.1se)
  pp  <- as.numeric(predict(fit, newx=xte, s=cv$lambda.1se, type="response"))
  pc  <- as.numeric(pp > 0.5)
  c(AUC_prob  = as.numeric(auc(roc(yte, pp, quiet=TRUE))),
    AUC_class = as.numeric(auc(roc(yte, pc, quiet=TRUE))),
    nonzero   = sum(as.matrix(coef(fit))[-1,1] != 0))
}

set.seed(98)
half <- sample(ncol(a), ncol(a) %/% 2)                  # subject split for the "split" design
out <- list()
for (mn in names(M)) {
  t0 <- Sys.time()
  r <- tryCatch(M[[mn]](a, s), error=function(e){cat("  ",mn,"ERR:",conditionMessage(e),"\n"); NULL})
  if (is.null(r)) { cat(sprintf("%-14s SKIPPED\n", mn)); next }
  ca <- as.matrix(r$array); cs <- as.matrix(r$seq)
  dimnames(ca) <- dimnames(a); dimnames(cs) <- dimnames(s)
  ok <- complete.cases(ca) & complete.cases(cs)
  ca <- ca[ok,,drop=FALSE]; cs <- cs[ok,,drop=FALSE]
  for (dir in c("array_to_seq","seq_to_array")) {
    tr <- if (dir=="array_to_seq") ca else cs
    te <- if (dir=="array_to_seq") cs else ca
    for (pn in names(ALPHA)) {
      pub <- fit_eval(t(tr), y, t(te), y, ALPHA[[pn]])                 # published design
      spl <- fit_eval(t(tr[,half]), y[half], t(te[,-half]), y[-half], ALPHA[[pn]])
      out[[length(out)+1]] <- data.frame(method=mn, direction=dir, penalty=pn,
        paired_AUC_prob=pub[["AUC_prob"]], paired_AUC_class=pub[["AUC_class"]],
        split_AUC_prob=spl[["AUC_prob"]], split_AUC_class=spl[["AUC_class"]],
        nonzero=pub[["nonzero"]], genes=nrow(ca))
    }
  }
  cat(sprintf("%-14s done (%ds)\n", mn, round(as.numeric(difftime(Sys.time(),t0,units="secs")))))
  flush.console()
}
res <- do.call(rbind, out)
write.csv(res, "/work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/results/metsim_prediction.csv", row.names=FALSE)
cat("\nwrote results/metsim_prediction.csv\n")
print(res[res$penalty=="lasso", c("method","direction","paired_AUC_class","paired_AUC_prob")], row.names=FALSE)
