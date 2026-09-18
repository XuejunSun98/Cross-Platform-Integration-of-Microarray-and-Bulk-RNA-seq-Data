# ============================================================
# Smoke test: does each of the compared methods actually RUN on Longleaf?
# One small simulated dataset (m=5 per group, 16160 genes, n_DE=1000),
# generated with the exact generator from data_simulation.R.
# Each method block is copied from the current comparison scripts.
# Reports: PASS/FAIL, wall time, and a sanity Type I error value.
# ============================================================
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
ROOT <- "/work/users/x/u/xuejun1/Integration_paper_Sim"
SIM  <- file.path(ROOT, "simulation_2024")
MMR_DIR <- file.path(ROOT, "Cross-Platform-Normalization-master/MatchMixeR/R")
setwd(SIM)

suppressPackageStartupMessages({
  library(preprocessCore); library(stringr); library(limma); library(sva)
})

# ---------- 1. simulate (generator verbatim from data_simulation.R) ----------
sim_array_seq <- function(m, n_DE, seed) {
  load(file.path(ROOT, "sim_BMI_control0.rda"))
  set.seed(seed)
  m0 <- sample(1:ncol(GSE70353_control),  2*m, replace = TRUE)
  m1 <- sample(1:ncol(GSE135134_control), 2*m, replace = TRUE)
  sim_array <- GSE70353_control[, m0]
  sim_seq   <- GSE135134_control[, m1]
  random_error_matrix <- matrix(rnorm(nrow(sim_array)*ncol(sim_array), 0, 0.3), nrow = nrow(sim_array))
  sim_array <- sim_array + random_error_matrix
  sim_seq   <- exp(log(sim_seq + 1) + random_error_matrix)
  colnames(sim_array) <- c(paste0("control",1:m,"_array"), paste0("case",1:m,"_array"))
  colnames(sim_seq)   <- c(paste0("control",1:m,"_seq"),   paste0("case",1:m,"_seq"))
  rownames(sim_array) <- paste0("gene_", 1:nrow(sim_array))
  rownames(sim_seq)   <- paste0("gene_", 1:nrow(sim_seq))
  lambda2 <- 2; sigma <- 0.5
  sim_array[1:n_DE,(m+1):(2*m)] <- sim_array[1:n_DE,(m+1):(2*m)] + abs(matrix(rnorm(n_DE*m, 0.5+lambda2*exp(-lambda2), sigma), ncol=m))*0.3
  sigma <- 1
  sim_seq[1:n_DE,(m+1):(2*m)]   <- sim_seq[1:n_DE,(m+1):(2*m)]   + exp(abs(matrix(rnorm(n_DE*m, 2.5+lambda2*exp(-lambda2), sigma), ncol=m)))*0.3
  sim_array[(n_DE+1):(2*n_DE),1:m]         <- sim_array[1:n_DE,(m+1):(2*m)]
  sim_seq[(n_DE+1):(2*n_DE),1:m]           <- sim_seq[1:n_DE,(m+1):(2*m)]
  sim_array[(n_DE+1):(2*n_DE),(m+1):(2*m)] <- sim_array[1:n_DE,1:m]
  sim_seq[(n_DE+1):(2*n_DE),(m+1):(2*m)]   <- sim_seq[1:n_DE,1:m]
  list(seed = seed, sim_array = sim_array, sim_seq = sim_seq)
}

ARG <- commandArgs(trailingOnly=TRUE)
M <- if (length(ARG)) as.integer(ARG[1]) else 5
N_DE <- 1000
cat("Simulating one dataset (m =", M, ", n_DE =", N_DE, ")...\n")
t0 <- Sys.time()
SIMDATA <- sim_array_seq(m = M, n_DE = N_DE, seed = 101)
cat("  done in", round(difftime(Sys.time(), t0, units="secs")), "s;",
    nrow(SIMDATA$sim_array), "genes x", ncol(SIMDATA$sim_array), "array +",
    ncol(SIMDATA$sim_seq), "seq samples\n\n")

# ---------- 2. scoring (verbatim from integration_comparison_FDR_Final.R) ----------
score <- function(d_array, d_seq, m = M, n_DE = N_DE) {
  d_array <- as.data.frame(d_array); d_seq <- as.data.frame(d_seq)
  n <- nrow(d_array)
  d_all <- cbind(d_array, d_seq)
  ctrl <- c(1:m, (2*m+1):(3*m)); case <- c((m+1):(2*m), (3*m+1):(4*m))
  X <- as.matrix(d_all)
  LFC <- log(rowMeans(X[, case, drop=FALSE]) / rowMeans(X[, ctrl, drop=FALSE]))
  p <- apply(X, 1, function(r) tryCatch(wilcox.test(r[ctrl], r[case])$p.value, error=function(e) NA))
  sig <- !is.na(p) & p <= 0.05 & !is.na(LFC)
  de_idx <- 1:(2*n_DE); null_idx <- setdiff(1:n, de_idx)
  list(TypeI = mean(sig[null_idx]), Power = mean(sig[de_idx]),
       naP = sum(is.na(p)), naLFC = sum(!is.finite(LFC)))
}

# ---------- 3. harness ----------
RESULTS <- list()
run <- function(label, f) {
  cat("---- ", label, " ----\n"); flush.console()
  t <- Sys.time()
  out <- tryCatch({
    r <- f(SIMDATA$sim_array, SIMDATA$sim_seq)
    s <- if (is.list(r) && !is.null(r$TypeI)) r else score(r$d_array, r$d_seq)
    list(status = "PASS", TypeI = s$TypeI, Power = s$Power,
         msg = if (!is.null(s$naP) && (s$naP > 0 || s$naLFC > 0))
                 sprintf("NA p=%d, non-finite LFC=%d", s$naP, s$naLFC) else "")
  }, error = function(e) list(status = "FAIL", TypeI = NA, Power = NA,
                              msg = paste(conditionMessage(e), collapse=" ")))
  if (is.null(out$msg)) out$msg <- ""
  out$secs <- round(as.numeric(difftime(Sys.time(), t, units = "secs")), 1)
  out$method <- label
  cat("     ", out$status, " ", out$secs, "s  TypeI=", format(out$TypeI, digits=3),
      " Power=", format(out$Power, digits=3), " ", substr(out$msg, 1, 160), "\n", sep="")
  flush.console()
  RESULTS[[label]] <<- out
}

# ---------- 4. methods ----------
run("None (no correction)", function(a, s) list(d_array = a, d_seq = s))

run("QN", function(a, s) {
  ds <- as.data.frame(normalize.quantiles.use.target(as.matrix(s), target = a$control1_array))
  da <- as.data.frame(normalize.quantiles.use.target(as.matrix(a), target = a$control1_array))
  dimnames(da) <- dimnames(a); dimnames(ds) <- dimnames(s)
  list(d_array = da, d_seq = ds)
})

run("TDM", function(a, s) {
  library(TDM); library(data.table)
  ds <- tdm_transform(ref_data    = data.table(cbind(gene = rownames(a), a)),
                      target_data = data.table(cbind(gene = rownames(s), s)))
  ds <- as.data.frame(ds); rownames(ds) <- ds$gene; ds <- ds[, -1]
  ds[] <- lapply(ds, as.numeric)
  list(d_array = a, d_seq = ds)
})

run("Angel (rank)", function(a, s) {
  d_all <- cbind(a, s)
  for (i in 1:ncol(d_all)) d_all[, i] <- rank(d_all[, i]) / length(d_all[, i])
  list(d_array = d_all[, 1:(2*M)], d_seq = d_all[, (2*M+1):(4*M)])
})

MMR_fit <- local({ cached <- NULL; function(a, s) {
  if (is.null(cached)) {
    owd <- getwd(); setwd(MMR_DIR)
    for (f in c("ckmeans.R","dwd.R","eb.R","functions.R","gq.R","MatchMixeR.R","xpn.R")) source(f)
    setwd(owd); set.seed(111)
    mm <- MM(as.matrix(a), as.matrix(log(s + 1)))
    cached <<- list(intercept = mm$betamat[,1] + mm$betahat[1],
                    coef      = mm$betamat[,2] + mm$betahat[2])
  }
  cached }})

run("MMR [as written in integration_comparison*: loop 1:ncol]", function(a, s) {
  b <- MMR_fit(a, s); ds <- s
  # loop bound is ncol(d_seq) but indexes rows -> only the first 4m genes get transformed
  for (i in 1:ncol(s)) ds[i, ] <- b$intercept[i] + b$coef[i] * as.numeric(s[i, ])
  list(d_array = a, d_seq = ds)
})

run("MMR [as in 2026 scripts: sweep, all genes, raw seq]", function(a, s) {
  b <- MMR_fit(a, s)
  ds <- sweep(sweep(s, 1, b$coef, "*"), 1, b$intercept, "+")
  list(d_array = a, d_seq = ds)
})

run("MMR [log-consistent: sweep on log(seq+1)]", function(a, s) {
  b <- MMR_fit(a, s)
  ds <- sweep(sweep(log(s + 1), 1, b$coef, "*"), 1, b$intercept, "+")
  list(d_array = a, d_seq = ds)
})

run("Combat [as written in integration_comparison*]", function(a, s) {
  d_all <- as.matrix(cbind(a, log(s + 1)))
  batch <- factor(c(rep(1, ncol(a)), rep(2, ncol(s))))
  d_Combat <- as.data.frame(sva::ComBat(dat = d_all, batch = batch, mod = model.matrix(~1, data.frame(batch))))
  list(d_array = d_all[, 1:(2*M)], d_seq = d_all[, (2*M+1):(4*M)])   # <- bug: d_all, not d_Combat
})

run("Combat [corrected split, as in 2026 scripts]", function(a, s) {
  d_all <- as.matrix(cbind(a, log(s + 1)))
  batch <- factor(c(rep(1, ncol(a)), rep(2, ncol(s))))
  d_Combat <- as.data.frame(sva::ComBat(dat = d_all, batch = batch, mod = model.matrix(~1, data.frame(batch))))
  list(d_array = d_Combat[, 1:(2*M)], d_seq = d_Combat[, (2*M+1):(4*M)])
})

run("Combat_seq [as written]", function(a, s) {
  d_all <- as.matrix(cbind(round(exp(a)), s))
  storage.mode(d_all) <- "integer"
  batch <- factor(c(rep(1, ncol(a)), rep(2, ncol(s))))
  d_Combat <- sva::ComBat_seq(counts = d_all, batch = batch)
  list(d_array = d_all[, 1:(2*M)], d_seq = d_all[, (2*M+1):(4*M)])   # <- same bug
})

run("Combat_seq [corrected split]", function(a, s) {
  d_all <- as.matrix(cbind(round(exp(a)), s))
  storage.mode(d_all) <- "integer"
  batch <- factor(c(rep(1, ncol(a)), rep(2, ncol(s))))
  d_Combat <- sva::ComBat_seq(counts = d_all, batch = batch)
  list(d_array = d_Combat[, 1:(2*M)], d_seq = d_Combat[, (2*M+1):(4*M)])
})

run("RNABC [as written: QN + ComBat_seq]", function(a, s) {
  s <- round(s)
  d_QN <- as.data.frame(normalize.quantiles.use.target(as.matrix(a), target = s[, 1]))
  dimnames(d_QN) <- dimnames(a)
  dmatrix <- as.matrix(cbind(d_QN, s)); storage.mode(dmatrix) <- "integer"
  batch <- factor(c(rep(1, ncol(a)), rep(2, ncol(s))))
  d <- sva::ComBat_seq(counts = dmatrix, batch = batch)
  list(d_array = d[, 1:(2*M)], d_seq = d[, (2*M+1):(4*M)])
})

run("RNABC [manuscript definition: QN + ComBat]", function(a, s) {
  d_QN <- as.data.frame(normalize.quantiles.use.target(as.matrix(log(s + 1)), target = a[[1]]))
  dimnames(d_QN) <- dimnames(s)
  d_all <- as.matrix(cbind(a, d_QN))
  batch <- factor(c(rep(1, ncol(a)), rep(2, ncol(s))))
  d <- as.data.frame(sva::ComBat(dat = d_all, batch = batch, mod = model.matrix(~1, data.frame(batch))))
  list(d_array = d[, 1:(2*M)], d_seq = d[, (2*M+1):(4*M)])
})

run("COCONUT", function(a, s) {
  library(COCONUT)
  group <- as.factor(c(rep(0, M), rep(1, M)))
  al <- list(pheno = as.data.frame(group), genes = a)
  sl <- list(pheno = as.data.frame(group), genes = log(s + 1))
  rownames(al$pheno) <- colnames(a); rownames(sl$pheno) <- colnames(s)
  al$pheno$platform <- "array"; sl$pheno$platform <- "seq"
  g <- COCONUT(GSEs = list(d_array = al, d_seq = sl), control.0.col = "group", byPlatform = FALSE)
  list(d_array = cbind(g$controlList$GSEs$d_array$genes, g$COCONUTList$d_array$genes),
       d_seq   = cbind(g$controlList$GSEs$d_seq$genes,   g$COCONUTList$d_seq$genes))
})

run("limma (removeBatchEffect)", function(a, s) {
  d_all <- as.matrix(cbind(a, log(s + 1)))
  batch <- factor(c(rep(1, ncol(a)), rep(2, ncol(s))))
  grp   <- factor(rep(c(rep("control", M), rep("case", M)), 2), levels = c("control","case"))
  d <- as.data.frame(limma::removeBatchEffect(d_all, batch = batch, design = model.matrix(~ grp)))
  list(d_array = d[, 1:(2*M)], d_seq = d[, (2*M+1):(4*M)])
})

run("Meta (MetaDE + ACAT)", function(a, s) {
  library(MetaDE); library(ACAT)
  ds <- round(s)
  meta <- data.frame(label = c(rep("control", M), rep("case", M)))
  mres <- MetaDE(data = list(d_array = a, d_seq = ds),
                 clin.data = list(meta_array = meta, meta_seq = meta),
                 data.type = c("continuous","discrete"), resp.type = "twoclass",
                 response = "label", meta.method = "Fisher", tail = "abs",
                 ind.method = c("limma","DESeq2"), paired = c(FALSE, FALSE),
                 parametric = TRUE, select.group = c("case","control"))
  # ACAT() errors on any NA; DESeq2 returns NA for Cook's-distance outliers.
  # Combine over the non-NA p-values per gene instead.
  P <- mres$ind.p
  p <- apply(P, 1, function(v) { v <- v[!is.na(v)]; if (!length(v)) NA_real_ else ACAT(v) })
  st <- mres$ind.stat
  ok   <- !is.na(p) & stats::complete.cases(st)
  up   <- which(ok & p < 0.05 & st[,1] < 0 & st[,2] < 0)
  down <- which(ok & p < 0.05 & st[,1] > 0 & st[,2] > 0)
  sig <- rep(FALSE, nrow(a)); sig[c(up, down)] <- TRUE
  de_idx <- 1:(2*N_DE); null_idx <- setdiff(1:nrow(a), de_idx)
  list(TypeI = mean(sig[null_idx]), Power = mean(sig[de_idx]),
       naP = sum(is.na(p)), naLFC = 0L)
})

# ---------- 5. external tools ----------
ext <- function(label, present, note) {
  RESULTS[[label]] <<- list(method = label, status = if (present) "PASS" else "BLOCKED",
                            TypeI = NA, Power = NA, secs = 0, msg = note)
  cat("---- ", label, " ----\n      ", RESULTS[[label]]$status, "  ", note, "\n", sep="")
}
rank_exe <- list.files(ROOT, pattern = "Rank-?In", recursive = TRUE, ignore.case = TRUE, full.names = TRUE)
rank_exe <- rank_exe[grepl("\\.exe$", rank_exe)]
ext("Rank_in", length(rank_exe) > 0,
    if (length(rank_exe)) paste("found", rank_exe) else "Rank-In.exe not present anywhere in the package (Windows binary; would not run on Linux regardless)")
mfile <- list.files(ROOT, pattern = "^Shambhala2\\.m$", recursive = TRUE, full.names = TRUE)
ext("Shambhala2", length(mfile) > 0,
    if (length(mfile)) paste("found", mfile) else "Shambhala2.m not present anywhere in the package; system() also points at a hardcoded local MATLAB path")

# ---------- 6. summary ----------
cat("\n\n================ SUMMARY ================\n")
tab <- do.call(rbind, lapply(RESULTS, function(r) data.frame(
  method = r$method, status = r$status, secs = r$secs,
  TypeI = round(r$TypeI, 4), Power = round(r$Power, 4),
  note = substr(r$msg, 1, 90), stringsAsFactors = FALSE)))
print(tab, row.names = FALSE)
saveRDS(tab, file.path(ROOT, sprintf("method_smoke_test_results_m%d.rds", M)))
write.csv(tab, file.path(ROOT, sprintf("method_smoke_test_results_m%d.csv", M)), row.names = FALSE)
