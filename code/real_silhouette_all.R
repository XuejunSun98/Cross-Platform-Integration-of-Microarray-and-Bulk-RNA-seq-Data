# R2.8, all methods: quantitative clustering metrics for the real-data evaluation,
# extended to the full method set (Shambhala2, Rank-In and COCONUT included) and to all
# four datasets. Supersedes real_silhouette.R, which covered the eleven unsupervised
# methods only.
#
# Two properties are scored on separate axes, because a method can succeed at one by
# destroying the other:
#   platform silhouette   -- near 0 means the two platforms are fully intermixed (good)
#   biological silhouette -- larger means samples of the same group stay together (good)
# plus the adjusted Rand index of hierarchical clustering against the biological labels.
# All on Euclidean distance over z-scored genes (see the note in score()).
#
# Every dataset is first reduced to at most CAP samples per platform, stratified by the
# biological label with a fixed seed. Shambhala2 costs ~3.4 min of MATLAB per RNA-seq
# sample, so the full METSIM (331) and LUSC (497) seq matrices are out of reach; capping
# keeps every method on *identical* data, which matters more here than sample count, and
# it also gives MMR the equal, paired column counts it requires.
#
# Each task also saves the dendrogram ordering behind its ARI, which is what the clustering
# tile figure (Supplementary Figure S6) draws, so the tiles and the metrics come from
# exactly the same corrected matrices.
#
# One array task per dataset; each writes one CSV plus one tile RDS.
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
source("/work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/code/sim_rerun_common.R")
source(file.path(ROOT, "revision_repo/code/longleaf_external_tools.R"))
suppressPackageStartupMessages({
  library(cluster); library(mclust); library(preprocessCore); library(limma)
  library(sva); library(MatchMixeR); library(TDM); library(batchelor); library(data.table)
  library(matrixStats); library(COCONUT)})
src <- grab(readLines(file.path(ROOT, "Shambhala2_upstream/Shambhala2.R"), warn = FALSE),
            "^Shambhala2 <- function")
src <- sub('system\\("matlab[^\n]*\\)', 'shambhala2_matlab()', src)
eval(parse(text = src), envir = globalenv())

CAP  <- 50
SEED <- 101

## ---- datasets: each gives array, seq (raw), and a biological label ----------
RA <- new.env(); load(file.path(ROOT, "real_all2.rda"), envir = RA)
load_seqc <- function() {
  a <- as.matrix(RA$d_affy_SEQC); s <- as.matrix(RA$d_seq_SEQC)
  list(array = a, seq = s, g_array = ifelse(grepl("^A", colnames(a)), "A", "B"),
       g_seq = ifelse(grepl("^A", colnames(s)), "A", "B")) }
load_ccle <- function() {
  a <- as.matrix(RA$d_affy_CCLE); s <- as.matrix(RA$d_seq_CCLE)
  list(array = a, seq = s, g_array = sub("_array_[0-9]+$", "", colnames(a)),
       g_seq = sub("_seq_[0-9]+$", "", colnames(s))) }
load_metsim <- function() {
  m  <- suppressWarnings(as.numeric(RA$meta_BMI[["matsuda:ch1"]])); ok <- !is.na(m)
  g  <- ifelse(m[ok] <= 4, "case", "control")
  list(array = as.matrix(RA$d_affy_BMI)[, ok, drop = FALSE],
       seq   = as.matrix(RA$d_seq_BMI)[, ok, drop = FALSE], g_array = g, g_seq = g) }
load_lusc <- function() {
  e <- new.env(); load(file.path(ROOT, "real_data_TCGA/LUSC_analysis.rda"), envir = e)
  a <- as.matrix(e$d_array); s <- as.matrix(e$d_seq)
  list(array = a, seq = s, g_array = ifelse(grepl("^control", colnames(a)), "I", "II+"),
       g_seq = ifelse(grepl("^control", colnames(s)), "I", "II+")) }
DATASETS <- list(SEQC = load_seqc, CCLE = load_ccle, METSIM = load_metsim, `TCGA-LUSC` = load_lusc)

# stratified cap; when the two platforms are the same subjects in the same order the same
# columns are kept on both sides, so MMR's pairing survives
strat <- function(g, cap) {
  if (length(g) <= cap) return(seq_along(g))
  by <- split(seq_along(g), g)
  # proportional allocation, at least 2 per group, adjusted to sum to exactly cap
  k <- pmax(2, round(cap * lengths(by) / length(g)))
  while (sum(k) > cap) { j <- which.max(k); k[j] <- k[j] - 1L }
  while (sum(k) < cap) { j <- which.max(lengths(by) - k); k[j] <- k[j] + 1L }
  sort(unlist(Map(function(ix, n) sample(ix, n), by, k), use.names = FALSE))
}
# The study design is paired, so where a sample can be identified on both platforms the two
# sides are restricted to the shared samples and kept in the same order. This is what MMR
# needs: it regresses one platform on the other sample by sample, so mismatched columns
# would produce a number that looks valid and means nothing. CCLE matches on cell line
# (107 array / 105 seq share 101) and TCGA-LUSC on barcode (129 of 132 array samples).
# METSIM has no shared column names but is the same 330 subjects in the same order, which
# the equal-size branch covers. SEQC pairs 4 array replicates against 8 sequencing lanes
# per condition and has no 1:1 correspondence, so MMR is not computable there and is
# reported as missing.
pkey <- function(x) sub("_(array|seq)_", "_", x)
subset_ds <- function(D) {
  set.seed(SEED)
  ka <- pkey(colnames(D$array)); ks <- pkey(colnames(D$seq))
  common <- intersect(ka, ks)
  if (length(common) >= 2) {
    ia_all <- match(common, ka); is_all <- match(common, ks)
    sel <- strat(D$g_array[ia_all], CAP)
    ia <- ia_all[sel]; is_ <- is_all[sel]; paired <- TRUE
  } else if (ncol(D$array) == ncol(D$seq) && identical(D$g_array, D$g_seq)) {
    ia <- strat(D$g_array, CAP); is_ <- ia; paired <- TRUE
  } else {
    # No shared sample identifiers and unequal counts: SEQC, where the two platforms measure
    # the same two MAQC reference RNA pools (A = UHRR, B = HBRR) with 4 microarray replicates
    # against 8 sequencing lanes per pool. Every replicate within a pool is the same RNA, so
    # the platforms are matched by truncating each biological group to the smaller of the two
    # counts; which replicate pairs with which is immaterial. This makes every dataset a
    # paired design and lets MMR run everywhere.
    ba <- split(seq_along(D$g_array), D$g_array); bs <- split(seq_along(D$g_seq), D$g_seq)
    gg <- intersect(names(ba), names(bs))
    n  <- pmin(lengths(ba[gg]), lengths(bs[gg]))
    ia <- unlist(Map(function(ix, k) sort(ix)[seq_len(k)], ba[gg], n), use.names = FALSE)
    is_ <- unlist(Map(function(ix, k) sort(ix)[seq_len(k)], bs[gg], n), use.names = FALSE)
    sel <- strat(D$g_array[ia], CAP)
    ia <- ia[sel]; is_ <- is_[sel]; paired <- TRUE
  }
  list(array = D$array[, ia, drop = FALSE], seq = D$seq[, is_, drop = FALSE],
       g_array = D$g_array[ia], g_seq = D$g_seq[is_], paired = paired)
}

## ---- methods ---------------------------------------------------------------
M <- list(
  No_correction = function(a,s,g) list(array=a, seq=log2(s+1)),
  QN = function(a,s,g) { tgt <- as.matrix(a)[,1]
    list(array=normalize.quantiles.use.target(as.matrix(a), target=tgt),
         seq  =normalize.quantiles.use.target(log2(as.matrix(s)+1), target=tgt)) },
  Angel = function(a,s,g) { d <- cbind(a, log2(s+1)); d <- apply(d, 2, function(x) rank(x)/length(x))
    list(array=d[,1:ncol(a),drop=FALSE], seq=d[,ncol(a)+seq_len(ncol(s)),drop=FALSE]) },
  TDM = function(a,s,g) { r <- tdm_transform(
      ref_data    = data.table(cbind(gene=rownames(a), as.data.frame(a))),
      target_data = data.table(cbind(gene=rownames(s), as.data.frame(round(s)))))
    m <- as.matrix(r[,-1]); rownames(m) <- r$gene; list(array=a, seq=m[rownames(a),,drop=FALSE]) },
  MMR        = function(a,s,g) METHODS$MMR(a,s),
  ComBat     = function(a,s,g) METHODS$ComBat(a,s),
  ComBat_seq = function(a,s,g) METHODS$ComBat_seq(a,s),
  RNABC      = function(a,s,g) METHODS$RNABC(a,s),
  limma      = function(a,s,g) METHODS$limma(a,s),
  XPN        = function(a,s,g) METHODS$XPN(a,s),
  MNN        = function(a,s,g) METHODS$MNN(a,s),
  Shambhala2 = function(a,s,g) {
    # RAW seq in, output used as is: Shambhala2 maps onto Q0's (log2 array) scale
    write.csv(data.frame(SYMBOL=rownames(a), a, check.names=FALSE), "Q0.csv", row.names=FALSE)
    write.csv(data.frame(SYMBOL=rownames(s), s, check.names=FALSE), "P0.csv", row.names=FALSE)
    file.copy("P0.csv", "Input.csv", overwrite = TRUE)
    H  <- Shambhala2("Input.csv", "P0.csv", "Q0.csv", delete_buffer_files = TRUE, k = 5)
    cs <- as.matrix(apply(H[, -1, drop=FALSE], 2, as.numeric)); rownames(cs) <- as.character(H[,1])
    list(array = a, seq = cs[rownames(s), , drop = FALSE]) },
  Rank_in = function(a,s,g) {
    # input format follows Unbalanced_Type_I_RankIn.R: leading "gene" column, no row
    # names, no quoting, RAW counts on the RNA-seq side
    # LUSC reuses the same TCGA barcode on both platforms; Rank-In matches the class file
    # to the expression header by name, so columns are made unique per platform first.
    aa <- a; ss <- s
    colnames(aa) <- paste0("array_", colnames(a)); colnames(ss) <- paste0("seq_", colnames(s))
    d_all <- cbind(gene = rownames(aa), as.data.frame(aa), as.data.frame(ss))
    write.table(d_all, "d_all_Rankin.txt", sep="\t", row.names=FALSE, quote=FALSE)
    write.table(data.frame(sample_name = colnames(d_all)[-1],
                           Class = as.integer(factor(g)) - 1L),
                "sample.txt", sep="\t", row.names=FALSE, quote=FALSE)
    RankIn(quiet = FALSE)
    r <- read.table("Rank_In_r_result.txt", header=TRUE, row.names=1, sep="\t", check.names=FALSE)
    list(array = as.matrix(r[, seq_len(ncol(a)), drop=FALSE]),
         seq   = as.matrix(r[, ncol(a) + seq_len(ncol(s)), drop=FALSE])) },
  COCONUT = function(a,s,g) {
    ga <- factor(as.integer(factor(g[seq_len(ncol(a))])) - 1L)
    gs <- factor(as.integer(factor(g[ncol(a) + seq_len(ncol(s))])) - 1L)
    la  <- list(pheno = data.frame(group=ga, platform="array", row.names=colnames(a)),
                genes = as.data.frame(a))
    ls_ <- list(pheno = data.frame(group=gs, platform="seq", row.names=colnames(s)),
                genes = as.data.frame(log(s + 1)))
    o <- COCONUT(GSEs=list(d_array=la, d_seq=ls_), control.0.col="group", byPlatform=FALSE)
    list(array = as.matrix(cbind(o$controlList$GSEs$d_array$genes, o$COCONUTList$d_array$genes)),
         seq   = as.matrix(cbind(o$controlList$GSEs$d_seq$genes,   o$COCONUTList$d_seq$genes))) })

## ---- metrics ---------------------------------------------------------------
# Per-gene two-sample Kolmogorov-Smirnov statistic between the two platforms: the
# quantitative counterpart of the overlaid histograms in Supplementary Figure S4, one value
# per gene, summarised by its median. Computed directly from the pooled ranks rather than
# through ks.test() so that 16,160 genes cost seconds rather than minutes.
ks_gene <- function(ca, cs) {
  na <- ncol(ca); ns <- ncol(cs)
  apply(cbind(ca, cs), 1, function(x) {
    if (anyNA(x)) return(NA_real_)
    o  <- order(x); xs <- x[o]
    cs <- cumsum(ifelse(o <= na, 1/na, -1/ns))
    # with ties the ECDF difference is only defined at the end of each tie block
    max(abs(cs[c(xs[-length(xs)] != xs[-1], TRUE)]))
  })
}

score <- function(ca, cs, ga, gs) {
  d <- cbind(ca, cs)
  ok <- complete.cases(d) & apply(d, 1, function(x) sd(x) > 0)
  d <- d[ok, , drop = FALSE]
  # Euclidean on standardised genes, NOT Spearman-correlation distance. Spearman is
  # invariant to monotone transformation, so QN, Angel and TDM -- which are exactly that --
  # cannot change it. Genes are z-scored so no single high-variance gene dominates, but the
  # metric stays sensitive to the values each method actually alters, and it matches the
  # Euclidean distance specified in the Methods.
  z <- t(scale(t(d))); z <- z[complete.cases(z), , drop = FALSE]
  D <- dist(t(z))
  plat <- factor(c(rep("array", ncol(ca)), rep("seq", ncol(cs))))
  bio  <- factor(c(ga, gs))
  sil  <- function(f) if (nlevels(f) < 2) NA_real_ else summary(silhouette(as.integer(f), D))$avg.width
  h    <- hclust(D, method = "average")
  hc   <- cutree(h, k = nlevels(bio))
  ks   <- ks_gene(ca, cs)
  # the sample ordering the dendrogram induces is what the clustering tile figure draws:
  # one column per sample, a biological bar above a platform bar
  list(metrics = c(sil_platform = sil(plat), sil_biological = sil(bio),
                   ARI = mclust::adjustedRandIndex(hc, bio),
                   KS_median = median(ks, na.rm = TRUE), KS_mean = mean(ks, na.rm = TRUE)),
       tile = data.frame(pos = seq_along(h$order), biological = as.character(bio)[h$order],
                         platform = as.character(plat)[h$order]))
}

## ---- run one dataset --------------------------------------------------------
i  <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
dn <- names(DATASETS)[i]
wd <- file.path(Sys.getenv("TMPDIR", unset = tempdir()), sprintf("sil_%d", i))
dir.create(wd, recursive = TRUE, showWarnings = FALSE); old <- setwd(wd); on.exit(setwd(old))

D <- subset_ds(DATASETS[[dn]]())
a <- D$array; s <- D$seq; g <- c(D$g_array, D$g_seq)
cat("===", dn, ":", nrow(a), "genes |", ncol(a), "array /", ncol(s), "seq | paired:",
    D$paired, "===\n"); flush.console()

out <- list(); tiles <- list()
for (mn in names(M)) {
  t0 <- Sys.time()
  z <- tryCatch({
    r  <- M[[mn]](a, s, g)
    ca <- as.matrix(r$array); cs <- as.matrix(r$seq)
    dimnames(ca) <- dimnames(a); dimnames(cs) <- dimnames(s)
    score(ca, cs, D$g_array, D$g_seq)
  }, error = function(e) { cat("   ", mn, "ERR:", substr(conditionMessage(e),1,70), "\n")
    list(metrics = c(sil_platform=NA_real_, sil_biological=NA_real_, ARI=NA_real_,
                     KS_median=NA_real_, KS_mean=NA_real_), tile = NULL) })
  v <- z$metrics
  if (!is.null(z$tile)) tiles[[length(tiles)+1]] <- cbind(z$tile, dataset = dn, method = mn)
  out[[length(out)+1]] <- data.frame(dataset = dn, method = mn, t(v),
    n_array = ncol(a), n_seq = ncol(s),
    secs = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
  cat(sprintf("  %-14s platform %7s  biological %7s  ARI %7s  KS %7s  (%ss)\n", mn,
      format(v[["sil_platform"]], digits=3), format(v[["sil_biological"]], digits=3),
      format(v[["ARI"]], digits=3), format(v[["KS_median"]], digits=3),
      round(as.numeric(difftime(Sys.time(), t0, units="secs"))))); flush.console()
}
res <- do.call(rbind, out)
o <- file.path(ROOT, "revision_repo/results_sil_all"); dir.create(o, showWarnings = FALSE)
write.csv(res, file.path(o, sprintf("sil_%s.csv", gsub("[^A-Za-z0-9]", "", dn))), row.names = FALSE)
saveRDS(do.call(rbind, tiles), file.path(o, sprintf("tile_%s.rds", gsub("[^A-Za-z0-9]", "", dn))))
cat("DONE", dn, "\n")
