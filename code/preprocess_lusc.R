# Preprocess TCGA-LUSC for the cross-platform arm.
#
# Gene mapping follows the project's documented Affymetrix pipeline
# (xuejunsun98.github.io/HR-VILAGE-3K3M/docs/Affymetrix_microarray_process.html):
#   * take the first symbol when a probe maps to several  (sub("///.*", "", .))
#   * standardise with HGNChelper::checkGeneSymbols()$Suggested.Symbol, which
#     resolves aliases and outdated symbols
#   * drop NA symbols
#   * where several rows collapse to one symbol, keep the row with the HIGHEST
#     ROW MEDIAN (slice_max(RowMedian)), not the mean
#   * quantile-normalise the microarray after RMA
#
# Platform notes:
#   * Firehose "gene_rma" is HT-HG-U133A RMA, log2, already summarised to gene
#     symbols by the Broad -- so there are no probes left to collapse, but the
#     symbols still need standardising to match the RNA-seq side.
#   * RSEM "raw_count" are EXPECTED counts and not integers (e.g. 941.95);
#     rounded so ComBat-seq and DESeq2 can be used. Unlike METSIM, which offered
#     only TPM, these are genuine counts.
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
suppressPackageStartupMessages({library(HGNChelper); library(preprocessCore)})
setwd("/work/users/x/u/xuejun1/Integration_paper_Sim/real_data_TCGA")
load("LUSC_raw_matrices.rda")

# ---- standardise symbols, collapse by highest row median --------------------
standardise <- function(m, what) {
  sym <- sub("///.*", "", rownames(m))
  sug <- suppressWarnings(checkGeneSymbols(sym, species = "human"))$Suggested.Symbol
  sug <- sub(" ///.*", "", sug)                     # checkGeneSymbols can return multiples
  keep <- !is.na(sug) & sug != ""
  m <- m[keep, , drop = FALSE]; sug <- sug[keep]
  rmed <- matrixStats::rowMedians(m)
  ord  <- order(sug, -rmed)                          # highest median first per symbol
  m    <- m[ord, , drop = FALSE]; sug <- sug[ord]
  dup  <- duplicated(sug)
  cat(sprintf("  %-6s %5d rows -> %5d symbols (%d unmapped, %d duplicate collapsed)\n",
              what, nrow(m) + sum(!keep), sum(!dup), sum(!keep), sum(dup)))
  m <- m[!dup, , drop = FALSE]; rownames(m) <- sug[!dup]
  m
}
cat("standardising gene symbols (HGNChelper):\n")
a <- standardise(arr, "array")
s <- standardise(seq, "seq")

g <- sort(intersect(rownames(a), rownames(s)))
cat("\noverlap genes:", length(g), "\n")
a <- a[g, , drop = FALSE]; s <- s[g, , drop = FALSE]

drop <- rowSums(s) == 0
cat("dropped (zero in every seq sample):", sum(drop), "  -> genes kept:", sum(!drop), "\n")
a <- a[!drop, , drop = FALSE]; s <- round(s[!drop, , drop = FALSE])
storage.mode(s) <- "integer"

# ---- quantile-normalise the microarray, per the pipeline --------------------
qa <- normalize.quantiles(as.matrix(a)); dimnames(qa) <- dimnames(a); a <- qa

# ---- outcome: stage I vs II+ ------------------------------------------------
grp <- setNames(clt$stage_grp, rownames(clt))
lab <- function(ids) ifelse(grp[ids] == "I", "control", "case")
a <- a[, !is.na(grp[colnames(a)]), drop = FALSE]
s <- s[, !is.na(grp[colnames(s)]), drop = FALSE]
d_array <- a; d_seq <- s
colnames(d_array) <- paste0(lab(colnames(a)), "_", colnames(a))
colnames(d_seq)   <- paste0(lab(colnames(s)), "_", colnames(s))
surv   <- clt[, c("os_time","os_event","stage_grp")]
paired <- intersect(colnames(a), colnames(s))

cat("\n--- analysis matrices ---\n")
cat("d_array:", dim(d_array), " log2 RMA + QN, range", round(range(d_array),2), "\n")
cat("d_seq  :", dim(d_seq),   " integer counts, range", range(d_seq), "\n")
cat("array stage I/II+:", sum(grepl("^control",colnames(d_array))), "/", sum(grepl("^case",colnames(d_array))), "\n")
cat("seq   stage I/II+:", sum(grepl("^control",colnames(d_seq))),   "/", sum(grepl("^case",colnames(d_seq))), "\n")
cat("paired patients:", length(paired), "\n")
save(d_array, d_seq, surv, paired, file = "LUSC_analysis.rda")
cat("\nsaved LUSC_analysis.rda\n")
