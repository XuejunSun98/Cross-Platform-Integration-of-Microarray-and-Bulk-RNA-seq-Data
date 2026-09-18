# METSIM: rebuild the paired microarray / RNA-seq dataset behind manuscript Table 2.
#
# Sources (public GEO):
#   GSE70353  microarray, Affymetrix HG-U219 (GPL13667), 770 samples, RMA by submitters.
#             NOTE: the manuscript states platform GPL20609 and the `primeviewhsentrezg`
#             custom CDF. That is wrong -- GPL20609/PrimeView CDFs cannot be applied to
#             U219 CELs. The U219 equivalent would be hgu219hsentrezg. We use the GEO
#             series-matrix RMA values and map probes with hgu219.db.
#   GSE135134 RNA-seq TPM, 434 samples (matches the manuscript: "RNA-seq data enter as TPM").
#
# Outcome (manuscript, METSIM data subsection): Matsuda Index; control >= 4, case < 4.
# The Matsuda values live in the GSE70353 sample characteristics only, and both series
# use the subject ID as the sample title, so the outcome is carried across by subject.
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
suppressPackageStartupMessages({library(AnnotationDbi); library(hgu219.db); library(org.Hs.eg.db)})
RAW <- "/work/users/x/u/xuejun1/Integration_paper_Sim/real_data_METSIM/raw"
OUT <- "/work/users/x/u/xuejun1/Integration_paper_Sim/real_data_METSIM"

collapse_by_median <- function(mat, sym) {           # highest row median wins, as in LUSC
  keep <- !is.na(sym) & sym != ""
  mat <- mat[keep, , drop = FALSE]; sym <- sym[keep]
  o <- order(sym, -matrixStats::rowMedians(mat))
  mat <- mat[o, , drop = FALSE]; sym <- sym[o]
  mat <- mat[!duplicated(sym), , drop = FALSE]
  rownames(mat) <- sym[!duplicated(sym)]
  mat[order(rownames(mat)), , drop = FALSE]
}

## ---- microarray -------------------------------------------------------------
ln  <- readLines(gzfile(file.path(RAW, "GSE70353_series_matrix.txt.gz")))
gsm <- strsplit(sub("^!Sample_geo_accession\t", "", grep("^!Sample_geo_accession", ln, value = TRUE)), "\t")[[1]]
gsm <- gsub('"', '', gsm)
ttl <- strsplit(sub("^!Sample_title\t", "", grep("^!Sample_title", ln, value = TRUE)), "\t")[[1]]
ttl <- gsub('"', '', ttl)
mat_line <- grep("matsuda", ln, value = TRUE)[1]
mts <- gsub('"', '', strsplit(sub("^!Sample_characteristics_ch1\t", "", mat_line), "\t")[[1]])
mts <- suppressWarnings(as.numeric(sub("^matsuda: ", "", mts)))
stopifnot(length(gsm) == length(ttl), length(ttl) == length(mts))
pheno <- data.frame(gsm = gsm, subject = ttl, matsuda = mts, stringsAsFactors = FALSE)

b <- grep("^!series_matrix_table_begin", ln); e <- grep("^!series_matrix_table_end", ln)
tab <- read.delim(text = paste(ln[(b + 1):(e - 1)], collapse = "\n"), row.names = 1, check.names = FALSE)
a <- as.matrix(tab); colnames(a) <- gsub('"', '', colnames(a))
cat("array raw:", dim(a), "\n")
sym_a <- AnnotationDbi::mapIds(hgu219.db, rownames(a), "SYMBOL", "PROBEID", multiVals = "first")
a <- collapse_by_median(a, unname(sym_a))
colnames(a) <- pheno$subject[match(colnames(a), pheno$gsm)]
cat("array after symbol collapse:", dim(a), "\n")

## ---- RNA-seq ----------------------------------------------------------------
s <- read.delim(gzfile(file.path(RAW, "GSE135134_TPM.txt.gz")), row.names = 1, check.names = FALSE)
s <- as.matrix(s)
cat("seq raw:", dim(s), "\n")
ens <- sub("\\..*$", "", rownames(s))
sym_s <- AnnotationDbi::mapIds(org.Hs.eg.db, ens, "SYMBOL", "ENSEMBL", multiVals = "first")
s <- collapse_by_median(s, unname(sym_s))
cat("seq after symbol collapse:", dim(s), "\n")

## ---- pair, define outcome ---------------------------------------------------
subj  <- intersect(colnames(a), colnames(s))
mt    <- pheno$matsuda[match(subj, pheno$subject)]
subj  <- subj[!is.na(mt)]; mt <- mt[!is.na(mt)]
genes <- intersect(rownames(a), rownames(s))
d_array <- a[genes, subj, drop = FALSE]
d_seq   <- s[genes, subj, drop = FALSE]
grp     <- ifelse(mt < 4, "case", "control")     # manuscript: control >= 4, case < 4
colnames(d_array) <- paste0(grp, "_", subj)
colnames(d_seq)   <- paste0(grp, "_", subj)
cat(sprintf("\nPAIRED: %d subjects, %d genes | cases %d controls %d\n",
            length(subj), length(genes), sum(grp == "case"), sum(grp == "control")))
cat("array range", round(range(d_array), 2), "| seq(TPM) range", round(range(d_seq), 2), "\n")
save(d_array, d_seq, grp, subj, mt, file = file.path(OUT, "METSIM_analysis.rda"))
cat("saved", file.path(OUT, "METSIM_analysis.rda"), "\n")
