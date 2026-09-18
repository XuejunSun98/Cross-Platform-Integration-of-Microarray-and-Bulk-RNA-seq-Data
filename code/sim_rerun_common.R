# ============================================================
# Shared setup for the simulation re-run (RUNLIST items 1 and 4).
#
# Stage-1 .rda files (1-8 GB) were never shipped, so every grid is regenerated
# IN-SCRIPT from sim_BMI_control0.rda using the project's own generators, with
# the exact arguments the original drivers used. Verified faithful: the MMR
# pilot reproduced the published MMR power at all 12 (m, seed) points exactly,
# which is only possible if the regenerated data is identical.
#
# NOTE on grabbing the generators: data_simulation_unbalanced.R defines
# sim_array_seq_conditions TWICE (lines 474 and 838). source()-ing the file
# would silently take the SECOND, which is the prediction variant. We extract
# the first by brace balance.
# ============================================================
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
ROOT <- "/work/users/x/u/xuejun1/Integration_paper_Sim"
SIM  <- file.path(ROOT, "simulation_2024")

grab <- function(l, pat, which = 1) {
  s <- grep(pat, l)[which]
  if (is.na(s)) stop("not found: ", pat)
  bal <- 0; op <- FALSE
  for (k in s:length(l)) {
    t  <- gsub("#.*$", "", l[k])
    no <- lengths(regmatches(t, gregexpr("\\{", t)))
    nc <- lengths(regmatches(t, gregexpr("\\}", t)))
    if (no > 0) op <- TRUE
    bal <- bal + no - nc
    if (op && bal == 0) return(paste(l[s:k], collapse = "\n"))
  }
  stop("unterminated: ", pat)
}

local({
  gl <- readLines(file.path(SIM, "data_simulation_unbalanced.R"), warn = FALSE)
  # sim_array_seq_conditions is defined TWICE and the two are NOT interchangeable:
  #   [1] line 474 -- m=100, balanced 50/50 vs imbalanced 64/36  -> POWER by balance
  #   [2] line 838 -- m= 20, balanced 10/10 vs imbalanced 13/ 7  -> TYPE I (small sample)
  # source()-ing the file keeps only the second. Both are extracted, under
  # distinct names. Using [1] for the Type I grid fails the verification gate.
  for (spec in list(list("^sim_array_seq_balanced <- function",    1, NULL),
                    list("^sim_array_seq_effectratio <- function", 1, NULL),
                    list("^sim_array_seq_conditions <- function",  1, "sim_cond_power"),
                    list("^sim_array_seq_conditions <- function",  2, "sim_cond_typeI"))) {
    src <- grab(gl, spec[[1]], spec[[2]])
    src <- sub('load\\(".*sim_BMI_control0\\.rda"\\)',
               sprintf('load("%s")', file.path(ROOT, "sim_BMI_control0.rda")), src)
    if (!is.null(spec[[3]]))
      src <- sub("^sim_array_seq_conditions <- function", paste(spec[[3]], "<- function"), src)
    eval(parse(text = src), envir = globalenv())
  }
  # The scorer is copy-pasted and locally edited across the scripts, and the two
  # copies are NOT equivalent -- verified by diff:
  #   power  copy (94 lines): excludes gene_1..gene_2000 as true DE, so
  #                           Type I = FP / 14160, and guards LFC with 1e-8.
  #   Type I copy (52 lines): Type I = FP / 16160 over ALL genes, no DE concept,
  #                           no guard.
  # Using the power copy on the Type I grid fails the verification gate.
  cl <- readLines(file.path(SIM, "unbalanced_power_sample_size.R"), warn = FALSE)
  src <- sub("^Type_I_error_Power <- function", "score_power_fn <- function",
             grab(cl, "^Type_I_error_Power <- function"))
  eval(parse(text = src), envir = globalenv())
  tl <- readLines(file.path(SIM, "Unbalanced_Type_I_small.R"), warn = FALSE)
  src <- sub("^Type_I_error_Power <- function", "score_typeI_fn <- function",
             grab(tl, "^Type_I_error_Power <- function"))
  eval(parse(text = src), envir = globalenv())
})

# ---- the four grids, with the original drivers' arguments -------------------
GRIDS <- list(
  # Type I: SMALL sample (m = 20 per platform), no DE, no added noise --
  # driver at data_simulation_unbalanced.R:1016, saved as sim_unbalanced_type1_small.rda
  typeI_balance = list(levels = c("balanced","same_imbalance","opp_imbalance","half_imbalance"),
                       gen = function(lv, sd) sim_cond_typeI(
                         m = 20, n_case_bal = 10, n_case_imb = 13,
                         condition = lv, n_DE = 0, seed = sd)),
  power_size    = list(levels = c(5, 10, 30, 50),
                       gen = function(lv, sd) sim_array_seq_balanced(
                         n_per_group = as.numeric(lv), n_DE = 1000, seed = sd)),
  power_balance = list(levels = c("balanced","same_imbalance","opp_imbalance","half_imbalance"),
                       gen = function(lv, sd) sim_cond_power(
                         m = 100, n_case_bal = 50, n_case_imb = 64,
                         condition = lv, n_DE = 1000, seed = sd)),
  power_effect  = list(levels = c(0.3, 0.5, 0.7, 1),
                       gen = function(lv, sd) sim_array_seq_effectratio(
                         n_per_group = 50, n_DE = 1000,
                         effect_size = as.numeric(lv), seed = sd)))

# ---- methods ---------------------------------------------------------------
# Each returns list(array = <corrected array>, seq = <corrected seq>), both on a
# common scale, exactly as the existing method blocks hand them to the scorer.

m_MMR <- function(a, s) {
  # CORRECTED usage (RUNLIST item 4): Yhat is the corrected ARRAY, and betamat
  # already contains betahat -- the published arm added betahat twice and then
  # applied the result to the seq matrix.
  lg <- log(s + 1)
  mm <- MatchMixeR::MM(a, lg)
  list(array = mm$Yhat, seq = lg)
}

m_MNN <- function(a, s) {
  lg  <- log(s + 1)
  res <- batchelor::mnnCorrect(a, lg, k = 20, cos.norm.in = FALSE, cos.norm.out = FALSE)
  cm  <- SummarizedExperiment::assay(res, "corrected")
  list(array = cm[, seq_len(ncol(a)), drop = FALSE],
       seq   = cm[, ncol(a) + seq_len(ncol(s)), drop = FALSE])
}

m_ComBatseq <- function(a, s, protect = TRUE) {
  # ComBat-seq needs integer counts in BOTH batches, but the microarray is log2
  # RMA. We map it back with 2^x and round -- pseudo-counts, disclosed as a
  # limitation in the manuscript.
  #
  # The superseded integration_comparison* scripts used round(exp(a)), which is
  # the wrong inverse for log2 data: exp(log2(x)) = x^1.443, inflating the top
  # of the range from ~1.6e4 to ~1.1e6. Those results never reached a figure.
  # Depth-matching is REQUIRED, not cosmetic. 2^x gives the microarray a median
  # library size of 2.9e6 against the seq's 7.5e5 -- a 3.9x difference that
  # ComBat-seq reads as real sequencing depth. Measured at m=50 seed 101:
  #   2^x alone   : Type I 0.0004, power 0.503, 64.0% of genes still completely
  #                 platform-separated (uncorrected is 75.6%) -- i.e. barely corrected,
  #                 and the near-zero Type I is the MMR artefact, not calibration.
  #   depth-matched: Type I 0.0416, power 0.753, separation 17.8%.
  ca <- 2^as.matrix(a); cs <- round(as.matrix(s))
  ca <- round(ca * (median(colSums(cs)) / median(colSums(ca))))
  storage.mode(ca) <- "integer"; storage.mode(cs) <- "integer"
  d  <- cbind(ca, cs)
  b  <- factor(c(rep(1, ncol(ca)), rep(2, ncol(cs))))
  # ComBat_seq has no `mod`; the equivalent slot is `group`, which it turns into
  # model.matrix(~group) when full_mod = TRUE. Passing it keeps ComBat-seq on the
  # same footing as ComBat rather than giving the two different information.
  out <- if (protect) sva::ComBat_seq(counts = d, batch = b, group = grp_of(colnames(d)))
         else         sva::ComBat_seq(counts = d, batch = b)
  list(array = log2(out[, seq_len(ncol(ca)), drop = FALSE] + 1),
       seq   = log2(out[, ncol(ca) + seq_len(ncol(cs)), drop = FALSE] + 1))
}

# ---- preprocessing shared by the ComBat / limma arms -----------------------
# WITHIN-platform QN: array samples normalized among themselves, seq samples
# among themselves. This is NOT RNABC's step -- RNABC maps the seq onto the
# microarray's distribution (cross-platform); this only removes sample-to-sample
# drift inside each platform, leaving the platform offset for the method to fix.
# Log first, then QN, on the seq side: QN is rank-based so ranks are unchanged
# either way, but the reference distribution differs (mean(log x) != log(mean x)),
# and log-then-QN is the conventional order.
qn_within <- function(a, s) {
  qa <- preprocessCore::normalize.quantiles(as.matrix(a))
  qs <- preprocessCore::normalize.quantiles(as.matrix(log(s + 1)))
  dimnames(qa) <- dimnames(a); dimnames(qs) <- dimnames(s)
  list(array = qa, seq = qs)
}

grp_of <- function(cn) factor(ifelse(grepl("control", cn, ignore.case = TRUE), "control", "case"))

# ComBat and limma, with within-platform QN. `protect` selects the design term:
#   TRUE  -> ~group, the recommended usage: the condition is shielded from being
#            absorbed into the batch estimate. Only batch is ever subtracted, so
#            this is not supervision -- the outcome defines nothing.
#   FALSE -> ~1, which ComBat discards as an all-ones column (== mod = NULL).
#            Reported as a sensitivity analysis.
m_ComBat <- function(a, s, protect = TRUE) {
  q <- qn_within(a, s)
  d <- as.matrix(cbind(q$array, q$seq))
  b <- factor(c(rep("array", ncol(q$array)), rep("seq", ncol(q$seq))))
  md <- if (protect) model.matrix(~ grp_of(colnames(d))) else model.matrix(~ 1, data.frame(x = rep(1, ncol(d))))
  out <- sva::ComBat(dat = d, batch = b, mod = md)
  list(array = out[, seq_len(ncol(q$array)), drop = FALSE],
       seq   = out[, ncol(q$array) + seq_len(ncol(q$seq)), drop = FALSE])
}

m_limma <- function(a, s, protect = TRUE) {
  q <- qn_within(a, s)
  d <- as.matrix(cbind(q$array, q$seq))
  b <- factor(c(rep("array", ncol(q$array)), rep("seq", ncol(q$seq))))
  des <- if (protect) model.matrix(~ grp_of(colnames(d))) else NULL
  out <- if (is.null(des)) limma::removeBatchEffect(d, batch = b)
         else              limma::removeBatchEffect(d, batch = b, design = des)
  list(array = out[, seq_len(ncol(q$array)), drop = FALSE],
       seq   = out[, ncol(q$array) + seq_len(ncol(q$seq)), drop = FALSE])
}

# RNABC as PUBLISHED (Pedersen 2018; code/reference/RNABC_pipeline_upstream.R):
#   normalize.quantiles.use.target(seq, target = rowMeans(<microarray reference>))
#   ComBat(cbind(...), batch)            <- two positional args, so mod = NULL
# Our previous code instead used joint normalize.quantiles() and mod = ~group.
# The log(seq+1) is retained, a deliberate deviation: the paper never mixes
# platforms the way we do, and every other arm here receives log-scale seq.
m_RNABC <- function(a, s) {
  am <- as.matrix(a); lg <- as.matrix(log(s + 1))
  qs <- preprocessCore::normalize.quantiles.use.target(lg, target = rowMeans(am))
  d  <- cbind(am, qs)
  b  <- factor(c(rep("array", ncol(am)), rep("seq", ncol(qs))))
  out <- sva::ComBat(dat = d, batch = b)          # mod = NULL, as published
  list(array = out[, seq_len(ncol(am)), drop = FALSE],
       seq   = out[, ncol(am) + seq_len(ncol(qs)), drop = FALSE])
}

m_XPN <- function(a, s) {
  lg <- log(s + 1)
  z  <- MatchMixeR::xpn(as.matrix(a), as.matrix(lg))
  list(array = z$x, seq = z$y)      # XPN corrects BOTH platforms
}

METHODS <- list(
  MMR         = m_MMR,
  MNN         = m_MNN,
  # ComBat-seq: group = NULL is primary, matching ComBat/limma/RNABC. Passing the
  # outcome labels is actively harmful when the platforms have different group
  # proportions: on TCGA-LUSC (array 72/60 vs seq 244/253) it gave a perfect
  # in-sample fit (AUC 1.000) and an INVERTED test AUC of 0.048; with group=NULL
  # the same model gives 0.541. The group-specific NB adjustment goes one way on
  # each platform. This is also why the simulation's ComBat-seq Type I error was
  # 0.0995 under opp_imbalance -- the only design with opposite ratios -- against
  # ~0.03 elsewhere.
  ComBat_seq  = function(a, s) m_ComBatseq(a, s, protect = FALSE),
  XPN         = m_XPN,
  # ComBat: mod = NULL is the primary arm, matching RNABC's published usage and
  # keeping both ComBat variants unsupervised in the strict sense.
  ComBat      = function(a, s) m_ComBat(a, s, protect = FALSE),
  # limma: ~1 as primary, matching ComBat and RNABC -- no arm uses outcome labels,
  # so the unsupervised classification is unambiguous. Measured identical to
  # ~group under balanced and same_imbalance; can only differ under
  # opp_imbalance / half_imbalance.
  limma       = function(a, s) m_limma(a, s, protect = FALSE),
  RNABC       = m_RNABC,
  # sensitivity: the opposite design choice for each
  ComBat_group    = function(a, s) m_ComBat(a, s, protect = TRUE),
  limma_group     = function(a, s) m_limma(a, s, protect = TRUE),
  ComBat_seq_group = function(a, s) m_ComBatseq(a, s, protect = TRUE))

# `grid` selects the scorer: the Type I grid has no DE genes and is scored by the
# Type I script's own function; every other grid by the power script's.
score_pair <- function(res, dn_a, dn_s, grid = "power_size") {
  da <- as.data.frame(res$array); ds <- as.data.frame(res$seq)
  dimnames(da) <- dn_a; dimnames(ds) <- dn_s
  if (identical(grid, "typeI_balance")) {
    o <- score_typeI_fn(da, ds)
    c(typeI = o$Type_I_error, power = NA_real_)
  } else {
    o <- score_power_fn(da, ds)
    c(typeI = o$Type_I_error, power = if (is.null(o$Power_overall)) NA_real_ else o$Power_overall)
  }
}

run_seed <- function(seed, methods, out_csv) {
  rows <- list(); t0 <- Sys.time()
  for (gname in names(GRIDS)) {
    g <- GRIDS[[gname]]
    for (lv in g$levels) {
      r <- g$gen(lv, seed)
      a <- as.matrix(r$sim_array); s <- as.matrix(r$sim_seq)
      dn_a <- dimnames(r$sim_array); dn_s <- dimnames(r$sim_seq)
      for (mname in methods) {
        t1 <- Sys.time()
        v <- tryCatch(score_pair(METHODS[[mname]](a, s), dn_a, dn_s, gname),
                      error = function(e) c(typeI = NA_real_, power = NA_real_))
        rows[[length(rows) + 1]] <- data.frame(
          grid = gname, level = as.character(lv), seed = seed, method = mname,
          Type_I_error = v[["typeI"]], Power = v[["power"]],
          secs = round(as.numeric(difftime(Sys.time(), t1, units = "secs")), 1))
        write.csv(do.call(rbind, rows), out_csv, row.names = FALSE)
        cat(sprintf("[%5.1f min] %-14s %-15s seed=%d %-11s tI=%-8s pw=%-8s %5.1fs\n",
                    as.numeric(difftime(Sys.time(), t0, units = "mins")),
                    gname, lv, seed, mname,
                    format(v[["typeI"]], digits = 3), format(v[["power"]], digits = 3),
                    rows[[length(rows)]]$secs)); flush.console()
      }
    }
  }
  invisible(do.call(rbind, rows))
}
