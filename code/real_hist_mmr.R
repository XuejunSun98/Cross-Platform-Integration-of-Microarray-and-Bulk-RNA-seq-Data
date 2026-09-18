# MMR densities for the distributional-alignment figure on the datasets where the fast
# script could not run it.
#
# MM() regresses one platform on the other sample by sample and so needs equal, matched
# columns. real_hist.R passes the full matrices, which are unequal for CCLE (107 array /
# 105 seq) and TCGA-LUSC (132 / 497) -- both failed with "non-conformable arrays". Here
# each of those is first restricted to the samples measured on both platforms, matched by
# cell line and by TCGA barcode respectively, and kept in the same order. All other methods
# keep the full matrices, so this restriction is noted in the figure caption.
#
# SEQC has no shared sample identifiers either, but its two platforms measure the same two
# MAQC reference RNA pools (A = UHRR, B = HBRR) with 4 microarray replicates against 8
# sequencing lanes per pool. Every replicate within a pool is the same RNA, so the platforms
# are matched by truncating each pool to the smaller count; which replicate pairs with which
# is immaterial.
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
source("/work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/code/sim_rerun_common.R")

RA <- new.env(); load(file.path(ROOT, "real_all2.rda"), envir = RA)
load_seqc <- function() list(array = as.matrix(RA$d_affy_SEQC), seq = as.matrix(RA$d_seq_SEQC),
  g_array = ifelse(grepl("^A", colnames(RA$d_affy_SEQC)), "A", "B"),
  g_seq   = ifelse(grepl("^A", colnames(RA$d_seq_SEQC)),  "A", "B"))
load_ccle <- function() list(array = as.matrix(RA$d_affy_CCLE), seq = as.matrix(RA$d_seq_CCLE))
load_lusc <- function() { e <- new.env(); load(file.path(ROOT, "real_data_TCGA/LUSC_analysis.rda"), envir = e)
  list(array = as.matrix(e$d_array), seq = as.matrix(e$d_seq)) }
DATASETS <- list(SEQC = load_seqc, CCLE = load_ccle, `TCGA-LUSC` = load_lusc)

pkey <- function(x) sub("_(array|seq)_", "_", x)
SUB  <- 2e5
dens <- function(x) { x <- x[is.finite(x)]; if (length(x) > SUB) x <- sample(x, SUB)
  d <- density(x, n = 512); data.frame(x = d$x, y = d$y) }

set.seed(101)
out <- file.path(ROOT, "revision_repo/results_hist_sup"); dir.create(out, showWarnings = FALSE)
for (dn in names(DATASETS)) {
  D <- DATASETS[[dn]](); a <- D$array; s <- D$seq
  ka <- pkey(colnames(a)); ks <- pkey(colnames(s)); common <- intersect(ka, ks)
  if (length(common) >= 2) {
    a <- a[, match(common, ka), drop = FALSE]; s <- s[, match(common, ks), drop = FALSE]
  } else {
    ba <- split(seq_along(D$g_array), D$g_array); bs <- split(seq_along(D$g_seq), D$g_seq)
    gg <- intersect(names(ba), names(bs)); n <- pmin(lengths(ba[gg]), lengths(bs[gg]))
    a <- a[, unlist(Map(function(ix,k) sort(ix)[seq_len(k)], ba[gg], n), use.names=FALSE), drop=FALSE]
    s <- s[, unlist(Map(function(ix,k) sort(ix)[seq_len(k)], bs[gg], n), use.names=FALSE), drop=FALSE]
  }
  cat("===", dn, ":", nrow(a), "genes |", ncol(a), "matched pairs ===\n"); flush.console()
  r  <- METHODS$MMR(a, s)
  res <- rbind(cbind(dens(as.vector(as.matrix(r$array))), platform = "Microarray",
                     dataset = dn, method = "MMR"),
               cbind(dens(as.vector(as.matrix(r$seq))),   platform = "RNA-seq",
                     dataset = dn, method = "MMR"))
  saveRDS(res, file.path(out, sprintf("hist_%s_MMR.rds", gsub("[^A-Za-z0-9]", "", dn))))
  cat("  wrote", dn, "\n"); flush.console()
}
cat("DONE\n")
