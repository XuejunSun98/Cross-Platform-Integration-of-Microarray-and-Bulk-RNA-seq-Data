# METSIM prediction under REPEATED STRATIFIED SUBJECT SPLITS.
#
# Why: METSIM is fully paired, so the published design trains on 330 subjects' arrays and
# tests on THE SAME 330 subjects' RNA-seq -- every test subject's label was used in
# training. The tell is that "no correction" is mid-table in both directions (published
# Table 2: 0.785 / 0.672 inside a 0.77-0.90 field). Splitting SUBJECTS removes that.
#
# A single 50/50 split is too noisy at 70 cases / 330 subjects to separate methods that
# differ by 0.02-0.10, so we repeat NREP stratified splits and report mean and SD.
#
# The correction is fitted ONCE on all 330 subjects (transductive), as in the manuscript
# and the LUSC analysis. That is standard in this literature but is a real caveat: the
# test samples' EXPRESSION (never their labels) informs the correction. It is disclosed,
# not hidden. One array task per method.
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
source("/work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/code/sim_rerun_common.R")
suppressPackageStartupMessages({
  library(glmnet); library(pROC); library(preprocessCore); library(limma)
  library(sva); library(MatchMixeR); library(TDM); library(data.table)})
load(file.path(ROOT, "real_data_METSIM/METSIM_analysis.rda"))
a <- as.matrix(d_array); s <- as.matrix(d_seq)
y <- ifelse(grepl("^case", colnames(a)), 1, 0)

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
  XPN = METHODS$XPN, MNN = METHODS$MNN, ComBat_seq = METHODS$ComBat_seq)

ti <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
mn <- names(M)[ti]
NREP <- 50; TRAIN_FRAC <- 0.7
cat(sprintf("method %s (%d/%d), %d reps, train frac %.2f\n", mn, ti, length(M), NREP, TRAIN_FRAC))

r  <- M[[mn]](a, s)
ca <- as.matrix(r$array); cs <- as.matrix(r$seq)
dimnames(ca) <- dimnames(a); dimnames(cs) <- dimnames(s)
ok <- complete.cases(ca) & complete.cases(cs)
ca <- ca[ok,,drop=FALSE]; cs <- cs[ok,,drop=FALSE]

auc_of <- function(xtr, ytr, xte, yte, alpha) {
  keep <- apply(xtr, 2, sd, na.rm = TRUE) > 0
  xtr <- xtr[, keep, drop=FALSE]; xte <- xte[, keep, drop=FALSE]
  cv  <- cv.glmnet(xtr, ytr, family="binomial", alpha=alpha, nfolds=10)
  fit <- glmnet(xtr, ytr, family="binomial", alpha=alpha, lambda=cv$lambda.1se)
  pp  <- as.numeric(predict(fit, newx=xte, s=cv$lambda.1se, type="response"))
  as.numeric(auc(roc(yte, pp, quiet=TRUE)))
}

set.seed(2026)
ALPHA <- c(lasso = 1, elastic_net = 0.5, ridge = 0)
rows <- list(); t0 <- Sys.time()
for (rep in seq_len(NREP)) {
  # stratify so every split keeps the 21% case prevalence
  tr_i <- c(sample(which(y==1), round(TRAIN_FRAC*sum(y==1))),
            sample(which(y==0), round(TRAIN_FRAC*sum(y==0))))
  te_i <- setdiff(seq_along(y), tr_i)
  for (dir in c("array_to_seq","seq_to_array")) {
    TR <- if (dir=="array_to_seq") ca else cs
    TE <- if (dir=="array_to_seq") cs else ca
    for (pn in names(ALPHA)) {
      v <- tryCatch(auc_of(t(TR[, tr_i]), y[tr_i], t(TE[, te_i]), y[te_i], ALPHA[[pn]]),
                    error = function(e) NA_real_)
      rows[[length(rows)+1]] <- data.frame(method=mn, rep=rep, direction=dir, penalty=pn, AUC=v)
    }
  }
  if (rep %% 5 == 0) { cat(sprintf("  rep %d/%d (%.1f min)\n", rep, NREP,
      as.numeric(difftime(Sys.time(), t0, units="mins")))); flush.console() }
}
out <- file.path(ROOT, "revision_repo/results_metsim_split")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
write.csv(do.call(rbind, rows), file.path(out, sprintf("%s.csv", mn)), row.names = FALSE)
cat("DONE", mn, "\n")
