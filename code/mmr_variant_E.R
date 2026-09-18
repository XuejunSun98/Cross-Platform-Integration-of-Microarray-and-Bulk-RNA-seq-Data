# ============================================================
# Variant E: MM() used as its vignette documents.
#
#   "Xmat and Ymat are m x n dimensional matched samples gene expression
#    matrices. The MM function ... gives the Yhat, the transformed values for
#    the gene expression values in Xmat. Here, only the values in platform X
#    will be changed for data integration."
#
# The scripts call MM(d_array, log(d_seq+1)) -- so X = array, and the ARRAY is
# what should be transformed. But they then apply the coefficients to d_seq.
#
# Variants compared, all on the same datasets and the paper's own scorer:
#   B   published Fig 1 : b0 + b1 * d_seq          (raw seq)
#   C   published Fig 3 : b0 + b1 * log(d_seq+1)
#   E   vignette-correct: Yhat as corrected ARRAY, log-seq left unchanged
#   E2  same as E but computed as b0 + b1*d_array  (sanity: should equal E)
# ============================================================
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
ROOT <- "/work/users/x/u/xuejun1/Integration_paper_Sim"
setwd(file.path(ROOT, "simulation_2024"))
suppressPackageStartupMessages(library(stringr))

grab <- function(l, pat) {
  s <- grep(pat, l)[1]; bal <- 0; op <- FALSE
  for (k in s:length(l)) {
    t <- gsub("#.*$", "", l[k])
    no <- lengths(regmatches(t, gregexpr("\\{", t))); nc <- lengths(regmatches(t, gregexpr("\\}", t)))
    if (no > 0) op <- TRUE
    bal <- bal + no - nc
    if (op && bal == 0) return(paste(l[s:k], collapse = "\n"))
  }
}
gl  <- readLines("data_simulation_unbalanced.R", warn = FALSE)
src <- grab(gl, "^sim_array_seq_balanced <- function")
src <- sub('load\\(".*sim_BMI_control0\\.rda"\\)',
           sprintf('load("%s")', file.path(ROOT, "sim_BMI_control0.rda")), src)
eval(parse(text = src))
cl <- readLines("unbalanced_power_sample_size.R", warn = FALSE)
eval(parse(text = grab(cl, "^Type_I_error_Power <- function")))

owd <- getwd(); setwd(file.path(ROOT, "Cross-Platform-Normalization-master/MatchMixeR/R"))
for (f in c("ckmeans.R","dwd.R","eb.R","functions.R","gq.R","MatchMixeR.R","xpn.R")) source(f)
setwd(owd)

load("d_plot_Power_Sample_Size.rda"); pub <- d_plot

SS <- c(5, 10, 30, 50); SEEDS <- 101:103
rows <- list()
for (m in SS) for (sd in SEEDS) {
  r <- sim_array_seq_balanced(n_per_group = m, n_DE = 1000, seed = sd)
  a <- r$sim_array; s <- r$sim_seq
  lg <- log(s + 1)

  mm <- MM(as.matrix(a), as.matrix(lg))
  b0 <- mm$betamat[,1] + mm$betahat[1]
  b1 <- mm$betamat[,2] + mm$betahat[2]

  sc <- function(da, ds) {
    da <- as.data.frame(da); ds <- as.data.frame(ds)
    dimnames(da) <- dimnames(a); dimnames(ds) <- dimnames(s)
    o <- Type_I_error_Power(da, ds)
    c(o$Type_I_error, o$Power_overall)
  }
  B  <- sc(a, sweep(sweep(s,  1, b1, "*"), 1, b0, "+"))          # published Fig 1
  C  <- sc(a, sweep(sweep(lg, 1, b1, "*"), 1, b0, "+"))          # published Fig 3
  E  <- sc(mm$Yhat, lg)                                          # vignette: corrected ARRAY
  E2 <- sc(sweep(sweep(a, 1, b1, "*"), 1, b0, "+"), lg)          # sanity check on E

  pubv <- pub$Power[pub$method=="MMR" & pub$sample_size==m & pub$seed==sd]
  rows[[length(rows)+1]] <- data.frame(
    m = m, seed = sd, published = pubv,
    B_power = B[2], C_power = C[2], E_power = E[2], E2_power = E2[2],
    E_typeI = E[1],
    lambdahat = mm$lambdahat,
    b1_med = median(b1), b1_sd = sd(b1),
    Yhat_sd = median(apply(mm$Yhat, 1, sd)), array_sd = median(apply(a, 1, sd)))
  cat(sprintf("m=%-3d seed=%d pub=%.4f  B=%.4f C=%.4f  E=%.4f (E2=%.4f)  lambda=%.3g  b1 med=%.3f sd=%.3f\n",
              m, sd, pubv, B[2], C[2], E[2], E2[2], mm$lambdahat, median(b1), sd(b1)))
  flush.console()
}
tab <- do.call(rbind, rows)
write.csv(tab, file.path(ROOT, "mmr_variant_E_results.csv"), row.names = FALSE)

cat("\n=============== MMR power by variant ===============\n")
print(tab[, c("m","seed","published","B_power","C_power","E_power")], row.names = FALSE, digits = 4)
cat("\n--- mean power by sample size ---\n")
print(aggregate(cbind(published, B_power, C_power, E_power) ~ m, tab,
                FUN = function(x) round(mean(x), 4)))
cat("\n--- diagnostics: is the per-gene random effect degenerate? ---\n")
print(aggregate(cbind(lambdahat, b1_med, b1_sd, Yhat_sd, array_sd) ~ m, tab,
                FUN = function(x) signif(mean(x), 4)))
cat("\n--- E vs E2 (should be identical) : max |diff| =",
    max(abs(tab$E_power - tab$E2_power)), "\n")

cat("\n=============== ranking at each sample size (E included) ===============\n")
others <- subset(pub, method != "MMR" & seed %in% SEEDS & sample_size %in% SS)
om <- aggregate(Power ~ method + sample_size, others, mean)
for (m in SS) {
  o <- om[om$sample_size == m, c("method","Power")]
  mine <- tab[tab$m == m, ]
  for (nm in c("published","B_power","C_power","E_power"))
    o <- rbind(o, data.frame(method = paste0("MMR[", sub("_power","",nm), "]"),
                             Power = mean(mine[[nm]], na.rm = TRUE)))
  o <- o[order(-o$Power), ]; o$rank <- seq_len(nrow(o))
  cat("\n--- m =", m, "---\n"); print(o, row.names = FALSE, digits = 4)
}
