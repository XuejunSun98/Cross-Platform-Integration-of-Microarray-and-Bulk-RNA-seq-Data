# Build TCGA-LUSC array / RNA-seq matrices + clinical for the cross-platform arm.
# Source: Broad Firehose stddata 2016_01_28 (GDC Legacy Archive, which held the
# Affymetrix data, has been retired; cBioPortal values are RSEM-processed, which
# would confound the preprocessing question in R2.5).
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
suppressPackageStartupMessages(library(data.table))
setwd("/work/users/x/u/xuejun1/Integration_paper_Sim/real_data_TCGA")
A <- Sys.glob("raw/gdac*transcriptome*/LUSC.transcriptome*data.data.txt")
R <- Sys.glob("raw/gdac*RSEM_genes*/LUSC.rnaseqv2*data.data.txt")

# ---- microarray: HT-HG-U133A, RMA ------------------------------------------
a  <- fread(A, skip = 2, header = FALSE)
ah <- strsplit(readLines(A, n = 1), "\t")[[1]][-1]
arr <- as.matrix(a[, -1, with = FALSE]); rownames(arr) <- a[[1]]; colnames(arr) <- ah
storage.mode(arr) <- "numeric"

# ---- RNA-seq: RSEM_genes, keep the raw_count column of each triplet ---------
r   <- fread(R, skip = 2, header = FALSE)
rh  <- strsplit(readLines(R, n = 1), "\t")[[1]][-1]
typ <- strsplit(readLines(R, n = 2)[2], "\t")[[1]][-1]
# typ/rh index the header WITHOUT its first field; r's column 1 is gene_id, so
# data column j of the header is r column j+1
keep <- which(typ == "raw_count")
seq <- as.matrix(r[, keep + 1L, with = FALSE]); rownames(seq) <- r[[1]]; colnames(seq) <- rh[keep]
storage.mode(seq) <- "numeric"

# RSEM gene ids are "SYMBOL|ENTREZ"; drop unmapped "?" symbols
sym <- sub("\\|.*$", "", rownames(seq))
ok  <- sym != "?" & !duplicated(sym)
seq <- seq[ok, ]; rownames(seq) <- sym[ok]

cat("array:", dim(arr), " range", round(range(arr, na.rm=TRUE), 2), "\n")
cat("seq  :", dim(seq), " range", round(range(seq, na.rm=TRUE), 1), "\n")

# ---- sample types: TCGA barcode position 14-15, 01 = primary tumour ---------
stype <- function(x) substr(x, 14, 15)
cat("\narray sample types:", paste(names(table(stype(colnames(arr)))),
    table(stype(colnames(arr))), sep="=", collapse="  "), "\n")
cat("seq   sample types:", paste(names(table(stype(colnames(seq)))),
    table(stype(colnames(seq))), sep="=", collapse="  "), "\n")

arr <- arr[, stype(colnames(arr)) == "01", drop = FALSE]
seq <- seq[, stype(colnames(seq)) == "01", drop = FALSE]
pid <- function(x) substr(x, 1, 12)
colnames(arr) <- pid(colnames(arr)); colnames(seq) <- pid(colnames(seq))
arr <- arr[, !duplicated(colnames(arr))]; seq <- seq[, !duplicated(colnames(seq))]

both <- intersect(colnames(arr), colnames(seq))
cat("\nprimary tumours -- array:", ncol(arr), " seq:", ncol(seq), " PAIRED:", length(both), "\n")

genes <- intersect(rownames(arr), rownames(seq))
cat("shared genes:", length(genes), "\n")

# ---- clinical ---------------------------------------------------------------
cl <- fread(Sys.glob("raw/gdac*Clinical_Pick_Tier1*/LUSC.clin.merged.picked.txt"), header = TRUE)
cl <- as.data.frame(cl); rownames(cl) <- cl[[1]]; cl <- cl[, -1, drop = FALSE]
colnames(cl) <- toupper(gsub("-", ".", colnames(cl)))
cat("\nclinical fields:", paste(head(rownames(cl), 20), collapse=", "), "\n")
save(arr, seq, cl, both, genes, file = "LUSC_raw_matrices.rda")
cat("\nsaved LUSC_raw_matrices.rda\n")

# ---- endpoints --------------------------------------------------------------
clt <- as.data.frame(t(cl), stringsAsFactors = FALSE)
rownames(clt) <- toupper(gsub("\\.", "-", rownames(clt)))
st <- tolower(clt$pathologic_stage)
clt$stage_grp <- ifelse(grepl("stage i[^ivx]*$|stage ia|stage ib", st), "I",
                 ifelse(grepl("stage ii|stage iii|stage iv", st), "II+", NA))
clt$os_time  <- suppressWarnings(as.numeric(ifelse(clt$vital_status == "1",
                   clt$days_to_death, clt$days_to_last_followup)))
clt$os_event <- suppressWarnings(as.numeric(clt$vital_status))

cat("\n---- pathologic stage, all LUSC patients ----\n"); print(table(st, useNA="ifany"))
cat("\n---- stage I vs II+ ----\n"); print(table(clt$stage_grp, useNA="ifany"))
for (nm in c("array","seq","paired")) {
  ids <- switch(nm, array = colnames(arr), seq = colnames(seq), paired = both)
  s <- clt[intersect(ids, rownames(clt)), ]
  cat(sprintf("\n%-7s n=%-4d  stage I/II+ = %d/%d   deaths %d/%d   median OS %.0f d\n",
      nm, length(ids), sum(s$stage_grp=="I",na.rm=TRUE), sum(s$stage_grp=="II+",na.rm=TRUE),
      sum(s$os_event==1,na.rm=TRUE), sum(!is.na(s$os_event)), median(s$os_time,na.rm=TRUE)))
}
save(arr, seq, clt, both, genes, file = "LUSC_raw_matrices.rda")
cat("\nre-saved with endpoints\n")
