# Distributional alignment across all four real datasets and all evaluated methods.
# Replaces hist_all_real_3datasets_sorted.png, which covered three datasets and the
# original method set only.
#
# For each dataset x method the corrected matrices are reduced to a density over the
# pooled values of each platform, so microarray and RNA-seq can be overlaid. Values are
# subsampled (200k per platform) purely to keep the density estimation tractable; the
# shapes are unchanged.
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
source("/work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/code/sim_rerun_common.R")
suppressPackageStartupMessages({
  library(preprocessCore); library(limma); library(sva); library(MatchMixeR)
  library(TDM); library(batchelor); library(data.table); library(COCONUT)})

RA <- new.env(); load(file.path(ROOT, "real_all2.rda"), envir = RA)
load_seqc <- function() list(array = as.matrix(RA$d_affy_SEQC), seq = as.matrix(RA$d_seq_SEQC))
load_ccle <- function() list(array = as.matrix(RA$d_affy_CCLE), seq = as.matrix(RA$d_seq_CCLE))
load_metsim <- function() {
  m <- suppressWarnings(as.numeric(RA$meta_BMI[["matsuda:ch1"]])); ok <- !is.na(m)
  list(array = as.matrix(RA$d_affy_BMI)[, ok, drop = FALSE],
       seq   = as.matrix(RA$d_seq_BMI)[, ok, drop = FALSE], grp = ifelse(m[ok] <= 4, "case", "control")) }
load_lusc <- function() { e <- new.env(); load(file.path(ROOT, "real_data_TCGA/LUSC_analysis.rda"), envir = e)
  list(array = as.matrix(e$d_array), seq = as.matrix(e$d_seq)) }
DATASETS <- list(SEQC = load_seqc, CCLE = load_ccle, METSIM = load_metsim, `TCGA-LUSC` = load_lusc)

M <- list(
  `No correction` = function(a,s) list(array=a, seq=log2(s+1)),
  QN = function(a,s) { tgt <- as.matrix(a)[,1]
    list(array=normalize.quantiles.use.target(as.matrix(a), target=tgt),
         seq  =normalize.quantiles.use.target(log2(as.matrix(s)+1), target=tgt)) },
  Angel = function(a,s) { d <- cbind(a, log2(s+1)); d <- apply(d, 2, function(x) rank(x)/length(x))
    list(array=d[,1:ncol(a),drop=FALSE], seq=d[,ncol(a)+seq_len(ncol(s)),drop=FALSE]) },
  TDM = function(a,s) { r <- tdm_transform(
      ref_data    = data.table(cbind(gene=rownames(a), as.data.frame(a))),
      target_data = data.table(cbind(gene=rownames(s), as.data.frame(round(s)))))
    m <- as.matrix(r[,-1]); rownames(m) <- r$gene; list(array=a, seq=m[rownames(a),,drop=FALSE]) },
  MMR = METHODS$MMR, ComBat = METHODS$ComBat, `ComBat-seq` = METHODS$ComBat_seq,
  RNABC = METHODS$RNABC, limma = METHODS$limma, XPN = METHODS$XPN, MNN = METHODS$MNN)

SUB <- 2e5
dens <- function(x) { x <- x[is.finite(x)]
  if (length(x) > SUB) x <- sample(x, SUB)
  d <- density(x, n = 512); data.frame(x = d$x, y = d$y) }

out <- list()
for (dn in names(DATASETS)) {
  D <- DATASETS[[dn]](); a <- D$array; s <- D$seq
  cat("\n===", dn, ":", nrow(a), "genes |", ncol(a), "array /", ncol(s), "seq ===\n"); flush.console()
  for (mn in names(M)) {
    t0 <- Sys.time()
    r <- tryCatch(M[[mn]](a, s), error = function(e) { cat("   ", mn, "ERR:", substr(conditionMessage(e),1,50), "\n"); NULL })
    if (is.null(r)) next
    ca <- as.matrix(r$array); cs <- as.matrix(r$seq)
    out[[length(out)+1]] <- rbind(
      cbind(dens(as.vector(ca)), platform = "Microarray", dataset = dn, method = mn),
      cbind(dens(as.vector(cs)), platform = "RNA-seq",    dataset = dn, method = mn))
    cat(sprintf("  %-14s ok (%ds)\n", mn, round(as.numeric(difftime(Sys.time(), t0, units="secs"))))); flush.console()
  }
}
res <- do.call(rbind, out)
saveRDS(res, file.path(ROOT, "revision_repo/results/real_hist_density.rds"))
cat("\nwrote results/real_hist_density.rds:", nrow(res), "rows,",
    length(unique(res$method)), "methods x", length(unique(res$dataset)), "datasets\n")
