# Cross-platform stage prediction on TCGA-LUSC (R1.4 + R2.3).
#
# Train on one platform, test on the other, for each correction method and each
# of three penalties (lasso / elastic net / ridge).
#
# Design notes:
#  * Endpoint is pathologic stage I vs II+. Within-platform 5-fold CV gives
#    AUC 0.619 on the RNA-seq (n=497) and 0.537 on the microarray (n=132), so
#    this is a WEAK signal and the array side is close to chance. Absolute AUCs
#    should be read with that ceiling in mind; the comparison between methods is
#    still like-for-like because every method sees the same data and folds.
#  * ONLY the RNA-seq -> microarray direction is REPORTED. The array platform has
#    132 samples and a within-platform CV AUC of 0.537, so it cannot support model
#    training: in the array -> seq direction lasso and ridge selected zero genes at
#    lambda.1se for most methods, giving AUC exactly 0.500 (no model, not a tie).
#    Both directions are still computed and written to the CSV as a record.
#  * 129 patients appear on BOTH platforms. Testing on them would score the model
#    on patients it trained on, with identical labels, so they are removed from
#    the TEST set. The overlapping version is reported alongside as a sensitivity.
#  * Supervised methods (COCONUT, Rank-In) are excluded, as in the manuscript's
#    own prediction analysis, because they use outcome labels during correction.
#  * Correction is fitted on train and test together (transductive), which is
#    standard in this literature but is a real caveat.
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
source("/work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/code/sim_rerun_common.R")
suppressPackageStartupMessages({
  library(glmnet); library(pROC); library(preprocessCore); library(limma)
  library(sva); library(MatchMixeR); library(batchelor); library(TDM); library(data.table)})
load(file.path(ROOT, "real_data_TCGA/LUSC_analysis.rda"))

a <- as.matrix(d_array); s <- as.matrix(d_seq)
pid <- function(x) sub("^[a-z]+_", "", x)
lab <- function(x) factor(ifelse(grepl("^control", x), "I", "II+"), levels = c("I","II+"))

# ---- the methods, each returning list(array=, seq=) on a common scale --------
M <- list(
  No_correction = function(a,s) list(array=a, seq=log2(s+1)),
  QN = function(a,s) { tgt <- as.matrix(a)[,1]
    list(array=normalize.quantiles.use.target(as.matrix(a), target=tgt),
         seq  =normalize.quantiles.use.target(log2(as.matrix(s)+1), target=tgt)) },
  Angel = function(a,s) { d <- cbind(a, log2(s+1))
    d <- apply(d, 2, function(x) rank(x)/length(x))
    list(array=d[,1:ncol(a),drop=FALSE], seq=d[,ncol(a)+seq_len(ncol(s)),drop=FALSE]) },
  TDM = function(a,s) { r <- tdm_transform(
      ref_data    = data.table(cbind(gene=rownames(a), as.data.frame(a))),
      target_data = data.table(cbind(gene=rownames(s), as.data.frame(s))))
    m <- as.matrix(r[,-1]); rownames(m) <- r$gene; list(array=a, seq=m[rownames(a),,drop=FALSE]) },
  MMR = METHODS$MMR, ComBat = METHODS$ComBat, ComBat_seq = METHODS$ComBat_seq,
  RNABC = METHODS$RNABC, limma = METHODS$limma, MNN = METHODS$MNN, XPN = METHODS$XPN)

ALPHA <- c(lasso = 1, elastic_net = 0.5, ridge = 0)
out <- list(); set.seed(1)

for (mn in names(M)) {
  t0 <- Sys.time()
  # MatchMixeR needs equal, paired sample counts (132 vs 497 here), so it errors
  # with "non-conformable arrays". Report the skip loudly -- an earlier version
  # skipped before logging and the method vanished from the results silently.
  r <- tryCatch(M[[mn]](a, s), error = function(e) { cat("  ", mn, "ERR:", conditionMessage(e), "\n"); NULL })
  if (is.null(r)) { cat(sprintf("%-14s SKIPPED\n", mn)); flush.console(); next }
  ca <- as.matrix(r$array); cs <- as.matrix(r$seq)
  dimnames(ca) <- dimnames(a); dimnames(cs) <- dimnames(s)
  ca <- ca[complete.cases(ca) & complete.cases(cs), , drop=FALSE]
  cs <- cs[rownames(ca), , drop=FALSE]
  corr_s <- round(as.numeric(difftime(Sys.time(), t0, units="secs")))

  for (dir in c("array_to_seq","seq_to_array")) {
    # 129 of the 132 array patients are also in the seq cohort. Dropping the
    # overlap from the TEST set leaves only 3 array patients, so instead the
    # overlapping patients are assigned to ONE side: they remain in whichever
    # platform is the test set and are removed from training. No patient is ever
    # in both, and no sample is wasted.
    tr <- if (dir=="array_to_seq") ca else cs
    te <- if (dir=="array_to_seq") cs else ca
    if (dir == "seq_to_array")
      tr <- tr[, !(pid(colnames(tr)) %in% pid(colnames(te))), drop=FALSE]
    ytr <- lab(colnames(tr)); yte <- lab(colnames(te))
    keep <- if (dir=="array_to_seq") !(pid(colnames(te)) %in% pid(colnames(tr)))
            else rep(TRUE, ncol(te))
    for (an in names(ALPHA)) {
      nz <- NA_integer_
      v <- tryCatch({
        fit <- cv.glmnet(t(tr), ytr, family="binomial", alpha=ALPHA[an], nfolds=5)
        nz <<- sum(as.matrix(coef(fit, s="lambda.min"))[-1,1] != 0)
        p   <- predict(fit, t(te), s="lambda.min", type="response")[,1]
        c(all = as.numeric(auc(roc(yte, p, quiet=TRUE, direction="<"))),
          nov = as.numeric(auc(roc(yte[keep], p[keep], quiet=TRUE, direction="<"))))
      }, error=function(e) c(all=NA_real_, nov=NA_real_))
      out[[length(out)+1]] <- data.frame(method=mn, direction=dir, penalty=an,
        AUC_no_overlap=v[["nov"]], AUC_all=v[["all"]],
        n_train=ncol(tr), n_test=sum(keep), genes=nrow(ca), corr_secs=corr_s,
        nonzero=nz)
    }
  }
  cat(sprintf("%-14s corrected in %4ds | a->s %.3f/%.3f/%.3f | s->a %.3f/%.3f/%.3f\n", mn, corr_s,
      out[[length(out)-5]]$AUC_no_overlap, out[[length(out)-4]]$AUC_no_overlap, out[[length(out)-3]]$AUC_no_overlap,
      out[[length(out)-2]]$AUC_no_overlap, out[[length(out)-1]]$AUC_no_overlap, out[[length(out)]]$AUC_no_overlap))
  flush.console()
  write.csv(do.call(rbind, out), file.path(ROOT,"revision_repo/results/lusc_prediction.csv"), row.names=FALSE)
}
tab <- do.call(rbind, out)
cat("\n=============== AUC, test patients not seen in training ===============\n")
print(reshape(tab[,c("method","direction","penalty","AUC_no_overlap")],
      idvar=c("method","direction"), timevar="penalty", direction="wide"), row.names=FALSE, digits=3)
cat("\nWithin-platform 5-fold CV ceiling: seq 0.619, array 0.537 -- read the absolute\nvalues against that, not against 1.0.\n")
