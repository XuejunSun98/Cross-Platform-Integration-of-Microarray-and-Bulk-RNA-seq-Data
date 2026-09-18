# Shambhala2, Rank-In and COCONUT densities for the distributional-alignment figure, run
# as a separate array task per (dataset, method): Shambhala2 shells out to MATLAB and
# Rank-In to the vendor's Python runner, and both write fixed filenames, so each needs its
# own working directory. COCONUT rides along here because, like the other two, it is
# supervised and needs the class labels the fast script does not carry.
#
# Rank-In needs a class label per sample. SEQC uses A/B, CCLE breast vs large intestine,
# METSIM the Matsuda dichotomy and LUSC stage I vs II+.
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
source("/work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/code/sim_rerun_common.R")
source(file.path(ROOT, "revision_repo/code/longleaf_external_tools.R"))
suppressPackageStartupMessages({library(matrixStats); library(COCONUT)})
src <- grab(readLines(file.path(ROOT, "Shambhala2_upstream/Shambhala2.R"), warn = FALSE),
            "^Shambhala2 <- function")
src <- sub('system\\("matlab[^\n]*\\)', 'shambhala2_matlab()', src)
eval(parse(text = src), envir = globalenv())

RA <- new.env(); load(file.path(ROOT, "real_all2.rda"), envir = RA)
get_ds <- function(dn) {
  if (dn == "SEQC") { a <- as.matrix(RA$d_affy_SEQC); s <- as.matrix(RA$d_seq_SEQC)
    g <- c(ifelse(grepl("^A", colnames(a)), 0, 1), ifelse(grepl("^A", colnames(s)), 0, 1))
  } else if (dn == "CCLE") { a <- as.matrix(RA$d_affy_CCLE); s <- as.matrix(RA$d_seq_CCLE)
    g <- c(as.integer(grepl("^breast", colnames(a))), as.integer(grepl("^breast", colnames(s))))
  } else if (dn == "METSIM") { m <- suppressWarnings(as.numeric(RA$meta_BMI[["matsuda:ch1"]])); ok <- !is.na(m)
    a <- as.matrix(RA$d_affy_BMI)[, ok, drop=FALSE]; s <- as.matrix(RA$d_seq_BMI)[, ok, drop=FALSE]
    g <- rep(as.integer(m[ok] <= 4), 2)
  } else { e <- new.env(); load(file.path(ROOT, "real_data_TCGA/LUSC_analysis.rda"), envir = e)
    a <- as.matrix(e$d_array); s <- as.matrix(e$d_seq)
    g <- c(as.integer(!grepl("^control", colnames(a))), as.integer(!grepl("^control", colnames(s)))) }
  list(array = a, seq = s, grp = g)
}

COMBOS <- expand.grid(dataset = c("SEQC","CCLE","METSIM","TCGA-LUSC"),
                      method  = c("Shambhala2","Rank-In","COCONUT"), stringsAsFactors = FALSE)
i  <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
cb <- COMBOS[i, ]; dn <- cb$dataset; mn <- cb$method
cat("task", i, ":", dn, mn, "\n"); flush.console()

wd <- file.path(Sys.getenv("TMPDIR", unset = tempdir()), sprintf("hist_%d", i))
dir.create(wd, recursive = TRUE, showWarnings = FALSE); old <- setwd(wd); on.exit(setwd(old))
D <- get_ds(dn); a <- D$array; s <- D$seq

if (mn == "Shambhala2") {
  # RAW seq in, output used as is -- Shambhala2 maps onto Q0's (log2 array) scale.
  # Shambhala2 costs ~3.4 min per Input sample, so METSIM (331) and LUSC (497) would need
  # a day each. The Input is therefore capped at 50 randomly chosen RNA-seq samples; the
  # calibration sets P0 and Q0 stay whole, so the transformation itself is unchanged, and
  # a 50-sample subset is ample for estimating the density of ~10^6 corrected values.
  set.seed(101)
  CAP <- 50
  si  <- if (ncol(s) > CAP) sort(sample(ncol(s), CAP)) else seq_len(ncol(s))
  write.csv(data.frame(SYMBOL = rownames(a), a, check.names = FALSE), "Q0.csv", row.names = FALSE)
  write.csv(data.frame(SYMBOL = rownames(s), s, check.names = FALSE), "P0.csv", row.names = FALSE)
  write.csv(data.frame(SYMBOL = rownames(s), s[, si, drop = FALSE], check.names = FALSE),
            "Input.csv", row.names = FALSE)
  H  <- Shambhala2("Input.csv", "P0.csv", "Q0.csv", delete_buffer_files = TRUE, k = 5)
  cs <- as.matrix(apply(H[, -1, drop = FALSE], 2, as.numeric)); rownames(cs) <- as.character(H[, 1])
  cs <- cs[rownames(s), , drop = FALSE]; ca <- a
} else if (mn == "Rank-In") {
  # Input format follows Unbalanced_Type_I_RankIn.R exactly: a leading "gene" column, no
  # row names, no quoting, and RAW (unlogged) counts on the RNA-seq side.
  # LUSC reuses the same TCGA barcode on both platforms; Rank-In matches the class file to
  # the expression header by name, so the columns are made unique per platform first.
  aa <- a; ss <- s
  colnames(aa) <- paste0("array_", colnames(a)); colnames(ss) <- paste0("seq_", colnames(s))
  d_all <- cbind(gene = rownames(aa), as.data.frame(aa), as.data.frame(ss))
  write.table(d_all, "d_all_Rankin.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(data.frame(sample_name = colnames(d_all)[-1], Class = D$grp),
              "sample.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  RankIn(quiet = FALSE)
  r  <- read.table("Rank_In_r_result.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
  ca <- as.matrix(r[, seq_len(ncol(a)), drop = FALSE])
  cs <- as.matrix(r[, ncol(a) + seq_len(ncol(s)), drop = FALSE])
} else {
  # COCONUT, as in fig2_pca.R: controls coded 0, log(s+1) on the RNA-seq side, and the
  # control and corrected blocks pasted back together per platform.
  ga <- as.factor(D$grp[seq_len(ncol(a))])
  gs <- as.factor(D$grp[ncol(a) + seq_len(ncol(s))])
  la <- list(pheno = data.frame(group = ga, platform = "array", row.names = colnames(a)),
             genes = as.data.frame(a))
  ls_ <- list(pheno = data.frame(group = gs, platform = "seq", row.names = colnames(s)),
              genes = as.data.frame(log(s + 1)))
  o  <- COCONUT(GSEs = list(d_array = la, d_seq = ls_), control.0.col = "group", byPlatform = FALSE)
  ca <- as.matrix(cbind(o$controlList$GSEs$d_array$genes, o$COCONUTList$d_array$genes))
  cs <- as.matrix(cbind(o$controlList$GSEs$d_seq$genes,   o$COCONUTList$d_seq$genes))
}

SUB <- 2e5
dens <- function(x) { x <- x[is.finite(x)]; if (length(x) > SUB) x <- sample(x, SUB)
  d <- density(x, n = 512); data.frame(x = d$x, y = d$y) }
res <- rbind(cbind(dens(as.vector(ca)), platform = "Microarray", dataset = dn, method = mn),
             cbind(dens(as.vector(cs)), platform = "RNA-seq",    dataset = dn, method = mn))
out <- file.path(ROOT, "revision_repo/results_hist_sup")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
saveRDS(res, file.path(out, sprintf("hist_%s_%s.rds", gsub("[^A-Za-z0-9]", "", dn), gsub("[^A-Za-z0-9]", "", mn))))
cat("DONE", dn, mn, "\n")
