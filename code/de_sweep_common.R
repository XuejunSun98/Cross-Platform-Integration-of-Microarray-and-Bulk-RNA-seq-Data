# R1.3: vary the TOTAL number of DE genes and the up:down ratio.
#
# The published generator injects n_DE up-regulated genes into rows 1:n_DE and
# n_DE down-regulated into rows (n_DE+1):(2*n_DE) -- always symmetric. Both the
# generator and the power scorer are generalised here by patching the EXTRACTED
# source, so the original scripts stay byte-identical and every published result
# remains reproducible from them.
#
# What this tests: QN, Angel, TDM and Rank-In map each sample onto a common
# distribution. Under symmetric DE that is harmless. As the ratio skews, the case
# samples' overall distribution genuinely shifts, and a distribution-matching
# method will remove that shift as though it were technical -- destroying real
# signal. ComBat and limma are gene-wise and should be unaffected.
source("/work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/code/sim_rerun_common.R")

# ---- generalised generator --------------------------------------------------
local({
  gl  <- readLines(file.path(SIM, "data_simulation_unbalanced.R"), warn = FALSE)
  src <- grab(gl, "^sim_array_seq_balanced <- function")
  src <- sub('load\\(".*sim_BMI_control0\\.rda"\\)',
             sprintf('load("%s")', file.path(ROOT, "sim_BMI_control0.rda")), src)
  src <- sub("^sim_array_seq_balanced <- function", "sim_array_seq_deratio <- function", src)
  src <- sub("n_DE = 100,", "n_up = 1000, n_down = 1000,", src, fixed = TRUE)
  # up block: rows 1:n_DE -> 1:n_up ; down block: (n_DE+1):(2*n_DE) -> (n_up+1):(n_up+n_down)
  # n_down = 0 (the all-up extreme) must not produce the DESCENDING sequence
  # (n_up+1):(n_up+0), which silently indexes two rows. Route the down block
  # through a guarded index that is integer(0) when there are no down genes.
  src <- gsub("(n_DE + 1):(2 * n_DE)", "DN_ROWS", src, fixed = TRUE)
  src <- gsub("1:n_DE", "1:n_up", src, fixed = TRUE)
  src <- gsub("n_DE * length(case_idx_array)", "n_this * length(case_idx_array)", src, fixed = TRUE)
  src <- gsub("n_DE * length(case_idx_seq)",   "n_this * length(case_idx_seq)",   src, fixed = TRUE)
  src <- sub("if (n_DE > 0) {",
             "if (n_up + n_down > 0) {\n    DN_ROWS <- if (n_down > 0) (n_up + 1):(n_up + n_down) else integer(0)",
             src, fixed = TRUE)
  # the two rnorm() draws differ in length between the up and down blocks
  src <- sub("n_this * length(case_idx_array)", "n_up * length(case_idx_array)", src, fixed = TRUE)
  src <- sub("n_this * length(case_idx_seq)",   "n_up * length(case_idx_seq)",   src, fixed = TRUE)
  src <- gsub("n_this * length(case_idx_array)", "n_down * length(case_idx_array)", src, fixed = TRUE)
  src <- gsub("n_this * length(case_idx_seq)",   "n_down * length(case_idx_seq)",   src, fixed = TRUE)
  eval(parse(text = src), envir = globalenv())
})

# ---- generalised scorer -----------------------------------------------------
local({
  cl  <- readLines(file.path(SIM, "unbalanced_power_sample_size.R"), warn = FALSE)
  src <- grab(cl, "^Type_I_error_Power <- function")
  src <- sub("^Type_I_error_Power <- function", "score_de_fn <- function", src)
  src <- sub("n_DE = 1000, alpha = 0.05", "n_up = 1000, n_down = 1000, alpha = 0.05", src, fixed = TRUE)
  src <- sub('true_up   <- paste0("gene_", 1:n_DE)',
             'true_up   <- paste0("gene_", 1:n_up)', src, fixed = TRUE)
  src <- sub('true_down <- paste0("gene_", (n_DE + 1):(2 * n_DE))',
             'true_down <- paste0("gene_", (n_up + 1):(n_up + n_down))', src, fixed = TRUE)
  eval(parse(text = src), envir = globalenv())
})

# ---- the sweep grid ---------------------------------------------------------
# total DE as a share of 16,160 genes: 500=3%, 1000=6%, 2000=12%, 4000=25%
# 2000 at 1:1 is the PUBLISHED scenario and serves as a cross-check.
SWEEP <- local({
  g <- expand.grid(total = c(500, 1000, 2000, 4000), ratio = c("1:1","2:1","4:1"),
                   stringsAsFactors = FALSE)
  g$r_up <- as.numeric(sub(":.*","",g$ratio)); g$r_dn <- as.numeric(sub(".*:","",g$ratio))
  g$n_up   <- round(g$total * g$r_up / (g$r_up + g$r_dn))
  g$n_down <- g$total - g$n_up
  g$level  <- paste0(g$total, "_", g$ratio)
  g[, c("level","total","ratio","n_up","n_down")]
})

GRIDS[["de_sweep"]] <- list(
  levels = SWEEP$level,
  gen = function(lv, sd) {
    k <- SWEEP[SWEEP$level == lv, ]
    sim_array_seq_deratio(n_per_group = 50, n_up = k$n_up, n_down = k$n_down, seed = sd)
  })

score_sweep <- function(res, dn_a, dn_s, lv) {
  k  <- SWEEP[SWEEP$level == lv, ]
  da <- as.data.frame(res$array); ds <- as.data.frame(res$seq)
  dimnames(da) <- dn_a; dimnames(ds) <- dn_s
  o <- score_de_fn(da, ds, n_up = k$n_up, n_down = k$n_down)
  c(typeI = o$Type_I_error, power = o$Power_overall,
    power_up = o$Power_up, power_down = o$Power_down)
}

# ---- methods for the sweep --------------------------------------------------
# QN, Angel and TDM are the ones the hypothesis is ABOUT -- they map samples onto
# a common distribution, so they are expected to degrade as the ratio skews.
# MMR is excluded: it requires matched samples, which holds here (m=50 both
# sides), so it IS included. Rank-In and COCONUT are supervised and excluded, as
# in the manuscript's own DE analysis.
suppressPackageStartupMessages({library(preprocessCore); library(data.table)})
SWEEP_METHODS <- c(
  list(
    QN = function(a, s) { tgt <- as.matrix(a)[,1]
      list(array = normalize.quantiles.use.target(as.matrix(a), target = tgt),
           seq   = normalize.quantiles.use.target(as.matrix(s), target = tgt)) },
    Angel = function(a, s) { d <- cbind(as.matrix(a), as.matrix(s))
      d <- apply(d, 2, function(x) rank(x)/length(x))
      list(array = d[, seq_len(ncol(a)), drop=FALSE],
           seq   = d[, ncol(a) + seq_len(ncol(s)), drop=FALSE]) },
    TDM = function(a, s) { r <- TDM::tdm_transform(
        ref_data    = data.table(cbind(gene = rownames(a), as.data.frame(a))),
        target_data = data.table(cbind(gene = rownames(s), as.data.frame(s))))
      m <- as.matrix(r[,-1]); rownames(m) <- r$gene
      list(array = as.matrix(a), seq = m[rownames(a), , drop=FALSE]) }),
  METHODS[c("MMR","ComBat","ComBat_seq","RNABC","limma","MNN","XPN")])
