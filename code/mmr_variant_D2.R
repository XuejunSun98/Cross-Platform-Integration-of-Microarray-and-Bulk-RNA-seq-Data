# ============================================================
# Phase 1b: a NON-degenerate direction-corrected MMR.
#
# Variant D (naive MM(log(seq+1), array)) returns all-NaN coefficients:
# ~3,176 of 16,160 genes are constant in the seq matrix, so OLS's
# beta1 = CovXY/VarX divides by zero, and MM()'s cov(betamat) across genes
# then propagates NaN to EVERY gene. Any earlier D numbers are meaningless.
#
# D2 fits on the non-degenerate genes and falls back to the global
# (fixed-effect) coefficients for the constant ones, so all genes are kept
# and the comparison against B stays on the same gene set.
# ============================================================
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
ROOT <- "/work/users/x/u/xuejun1/Integration_paper_Sim"
setwd(file.path(ROOT, "simulation_2024"))
suppressPackageStartupMessages(library(stringr))

gen <- readLines("data_simulation_unbalanced.R")[838:1006]
gen <- sub('load\\(".*sim_BMI_control0\\.rda"\\)',
           sprintf('load("%s")', file.path(ROOT, "sim_BMI_control0.rda")), gen)
eval(parse(text = paste(gen, collapse = "\n")))

cons <- readLines("Unbalanced_Type_I_small.R")
s <- grep("^Type_I_error_Power <- function", cons)[1]
bal <- 0; e <- NA
for (k in s:length(cons)) {
  txt <- gsub("#.*$", "", cons[k])
  bal <- bal + lengths(regmatches(txt, gregexpr("\\{", txt))) -
               lengths(regmatches(txt, gregexpr("\\}", txt)))
  if (k > s && bal == 0) { e <- k; break }
}
eval(parse(text = paste(cons[s:e], collapse = "\n")))

owd <- getwd(); setwd(file.path(ROOT, "Cross-Platform-Normalization-master/MatchMixeR/R"))
for (f in c("ckmeans.R","dwd.R","eb.R","functions.R","gq.R","MatchMixeR.R","xpn.R")) source(f)
setwd(owd)

B_published <- function(d_array, d_seq) {          # what the paper actually ran
  mm <- MM(as.matrix(d_array), as.matrix(log(d_seq + 1)))
  b0 <- mm$betamat[,1] + mm$betahat[1]; b1 <- mm$betamat[,2] + mm$betahat[2]
  sweep(sweep(d_seq, 1, b1, "*"), 1, b0, "+")
}

C_scalefix <- function(d_array, d_seq) {           # same direction, right scale
  mm <- MM(as.matrix(d_array), as.matrix(log(d_seq + 1)))
  b0 <- mm$betamat[,1] + mm$betahat[1]; b1 <- mm$betamat[,2] + mm$betahat[2]
  sweep(sweep(log(d_seq + 1), 1, b1, "*"), 1, b0, "+")
}

D2_direction <- function(d_array, d_seq) {         # seq -> array, degeneracy handled
  lg <- log(d_seq + 1)
  v  <- apply(lg, 1, var)
  ok <- is.finite(v) & v > 0
  mm <- MM(as.matrix(lg[ok, , drop = FALSE]), as.matrix(d_array[ok, , drop = FALSE]))
  b0 <- mm$betamat[,1] + mm$betahat[1]; b1 <- mm$betamat[,2] + mm$betahat[2]
  out <- lg
  out[ok, ]  <- sweep(sweep(lg[ok, , drop = FALSE], 1, b1, "*"), 1, b0, "+")
  # constant genes: global fixed-effect line only
  out[!ok, ] <- lg[!ok, , drop = FALSE] * mm$betahat[2] + mm$betahat[1]
  attr(out, "n_degenerate") <- sum(!ok)
  attr(out, "nan_coefs")    <- sum(!is.finite(b1))
  out
}

load("d_plot_Type_I_error_small.rda")
target <- d_plot[d_plot$method == "MMR", ]

CONDS <- c("balanced", "same_imbalance", "opp_imbalance", "half_imbalance")
SEEDS <- 101:105
rows <- list()
for (cond in CONDS) for (sd in SEEDS) {
  r <- sim_array_seq_conditions(m = 20, n_case_bal = 10, n_case_imb = 13,
                                condition = cond, n_DE = 0, seed = sd)
  a <- r$sim_array; sq <- r$sim_seq
  sc <- function(f) {
    ds <- as.data.frame(f(a, sq)); dimnames(ds) <- dimnames(sq)
    Type_I_error_Power(as.data.frame(a), ds)$Type_I_error
  }
  d2 <- D2_direction(a, sq)
  rows[[length(rows)+1]] <- data.frame(
    condition = cond, seed = sd,
    published = target$Type_I_error[target$balance == cond & target$seed == sd],
    B_published = sc(B_published), C_scalefix = sc(C_scalefix), D2_direction = sc(D2_direction),
    n_degenerate = attr(d2, "n_degenerate"), nan_coefs = attr(d2, "nan_coefs"))
  cat("."); flush.console()
}
tab <- do.call(rbind, rows)
cat("\n\n=============== MMR: published vs corrected variants ===============\n")
print(tab, row.names = FALSE, digits = 6)
cat("\n--- mean Type I error (nominal 0.05) ---\n")
print(round(colMeans(tab[, c("published","B_published","C_scalefix","D2_direction")]), 4))
cat("\n--- mean by condition ---\n")
print(aggregate(cbind(published, C_scalefix, D2_direction) ~ condition, tab,
                FUN = function(x) round(mean(x), 4)))
cat("\n--- degenerate genes per fit (constant in seq) ---\n")
print(summary(tab$n_degenerate))
cat("NaN coefficients after filtering:", sum(tab$nan_coefs), "\n")
write.csv(tab, file.path(ROOT, "mmr_variant_D2_results.csv"), row.names = FALSE)
