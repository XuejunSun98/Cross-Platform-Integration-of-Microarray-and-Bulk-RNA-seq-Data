# R2.8: quantitative metrics for the real-data evaluation, replacing inspection of
# histograms and clustering tiles.
#
# Two properties are scored on separate axes, because a method can succeed at one by
# destroying the other:
#   platform silhouette   -- near 0 means the two platforms are fully intermixed (good)
#   biological silhouette -- larger means samples of the same group stay together (good)
# plus the adjusted Rand index of hierarchical clustering against the biological labels.
# All on Euclidean distance over z-scored genes (see the note in score()).
#
# Datasets: SEQC, CCLE, METSIM and TCGA-LUSC.
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
source("/work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/code/sim_rerun_common.R")
suppressPackageStartupMessages({
  library(cluster); library(mclust); library(preprocessCore); library(limma)
  library(sva); library(MatchMixeR); library(TDM); library(batchelor); library(data.table)})

## ---- datasets: each gives array, seq (raw), and a biological label ----------
# SEQC, CCLE and METSIM come from real_all2.rda, the authors' own processed matrices
# (METSIM there is the 16,160-gene version behind the published results, not the GEO
# rebuild). LUSC was assembled for this revision.
RA <- new.env(); load(file.path(ROOT, "real_all2.rda"), envir = RA)

load_seqc <- function() {
  a <- as.matrix(RA$d_affy_SEQC); s <- as.matrix(RA$d_seq_SEQC)
  list(array = a, seq = s,
       g_array = ifelse(grepl("^A", colnames(a)), "A", "B"),
       g_seq   = ifelse(grepl("^A", colnames(s)), "A", "B"))
}
load_ccle <- function() {
  a <- as.matrix(RA$d_affy_CCLE); s <- as.matrix(RA$d_seq_CCLE)
  list(array = a, seq = s,
       g_array = sub("_array_[0-9]+$", "", colnames(a)),
       g_seq   = sub("_seq_[0-9]+$",   "", colnames(s)))
}
load_metsim <- function() {
  a <- as.matrix(RA$d_affy_BMI); s <- as.matrix(RA$d_seq_BMI)
  # Matsuda Index: <= 4 insulin resistant (case), > 4 control -- as in the original code.
  # One subject has no Matsuda value and is dropped from both platforms.
  m  <- suppressWarnings(as.numeric(RA$meta_BMI[["matsuda:ch1"]]))
  ok <- !is.na(m)
  g  <- ifelse(m[ok] <= 4, "case", "control")
  list(array = a[, ok, drop = FALSE], seq = s[, ok, drop = FALSE], g_array = g, g_seq = g)
}
load_lusc <- function() {
  e <- new.env(); load(file.path(ROOT, "real_data_TCGA/LUSC_analysis.rda"), envir = e)
  a <- as.matrix(e$d_array); s <- as.matrix(e$d_seq)
  list(array = a, seq = s,
       g_array = ifelse(grepl("^control", colnames(a)), "I", "II+"),
       g_seq   = ifelse(grepl("^control", colnames(s)), "I", "II+"))
}
DATASETS <- list(SEQC = load_seqc, CCLE = load_ccle, METSIM = load_metsim, `TCGA-LUSC` = load_lusc)

## ---- methods ---------------------------------------------------------------
M <- list(
  No_correction = function(a,s) list(array=a, seq=log2(s+1)),
  QN = function(a,s) { tgt <- as.matrix(a)[,1]
    list(array=normalize.quantiles.use.target(as.matrix(a), target=tgt),
         seq  =normalize.quantiles.use.target(log2(as.matrix(s)+1), target=tgt)) },
  Angel = function(a,s) { d <- cbind(a, log2(s+1)); d <- apply(d, 2, function(x) rank(x)/length(x))
    list(array=d[,1:ncol(a),drop=FALSE], seq=d[,ncol(a)+seq_len(ncol(s)),drop=FALSE]) },
  TDM = function(a,s) { r <- tdm_transform(
      ref_data    = data.table(cbind(gene=rownames(a), as.data.frame(a))),
      target_data = data.table(cbind(gene=rownames(s), as.data.frame(round(s)))))
    m <- as.matrix(r[,-1]); rownames(m) <- r$gene; list(array=a, seq=m[rownames(a),,drop=FALSE]) },
  MMR = METHODS$MMR, ComBat = METHODS$ComBat, RNABC = METHODS$RNABC, limma = METHODS$limma,
  XPN = METHODS$XPN, MNN = METHODS$MNN, ComBat_seq = METHODS$ComBat_seq)

## ---- metrics ---------------------------------------------------------------
score <- function(ca, cs, ga, gs) {
  d <- cbind(ca, cs)
  ok <- complete.cases(d) & apply(d, 1, function(x) sd(x) > 0)
  d <- d[ok, , drop = FALSE]
  # Euclidean on standardised genes, NOT Spearman-correlation distance. Spearman is
  # invariant to monotone transformation, so QN, Angel and TDM -- which are exactly that --
  # cannot change it: they returned silhouettes identical to the uncorrected data to three
  # decimals in all four datasets. Genes are z-scored so no single high-variance gene
  # dominates the distance, but the metric remains sensitive to the values each method
  # actually alters. This also matches the Methods, which specify Euclidean distance.
  z <- t(scale(t(d)))
  z <- z[complete.cases(z), , drop = FALSE]
  D <- dist(t(z))
  plat <- factor(c(rep("array", ncol(ca)), rep("seq", ncol(cs))))
  bio  <- factor(c(ga, gs))
  sil <- function(f) if (nlevels(f) < 2) NA_real_ else
    summary(silhouette(as.integer(f), D))$avg.width
  hc <- cutree(hclust(D, method = "average"), k = nlevels(bio))
  c(sil_platform = sil(plat), sil_biological = sil(bio),
    ARI = mclust::adjustedRandIndex(hc, bio))
}

out <- list()
for (dn in names(DATASETS)) {
  D <- DATASETS[[dn]](); a <- D$array; s <- D$seq
  cat("\n===", dn, ":", nrow(a), "genes |", ncol(a), "array /", ncol(s), "seq ===\n"); flush.console()
  for (mn in names(M)) {
    t0 <- Sys.time()
    v <- tryCatch({
      r  <- M[[mn]](a, s)
      ca <- as.matrix(r$array); cs <- as.matrix(r$seq)
      dimnames(ca) <- dimnames(a); dimnames(cs) <- dimnames(s)
      score(ca, cs, D$g_array, D$g_seq)
    }, error = function(e) { cat("   ", mn, "ERR:", substr(conditionMessage(e),1,60), "\n")
                             c(sil_platform=NA_real_, sil_biological=NA_real_, ARI=NA_real_) })
    out[[length(out)+1]] <- data.frame(dataset = dn, method = mn, t(v),
      secs = round(as.numeric(difftime(Sys.time(), t0, units="secs")), 1))
    cat(sprintf("  %-14s platform %7s  biological %7s  ARI %7s\n", mn,
        format(v[["sil_platform"]], digits=3), format(v[["sil_biological"]], digits=3),
        format(v[["ARI"]], digits=3))); flush.console()
  }
}
res <- do.call(rbind, out)
write.csv(res, file.path(ROOT, "revision_repo/results/real_silhouette.csv"), row.names = FALSE)
cat("\nwrote results/real_silhouette.csv\n")
