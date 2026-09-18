# ============================================================
# RUNLIST item 4 (R3.m1) -- MMR re-run with the CORRECT usage,
# full seed set (101:200), all sample sizes.
#
# The published MMR arm is wrong in two ways:
#   (i)  betamat ALREADY contains betahat
#          flmer: betamat <- t(rbind(gamma0hat, gamma1hat) + betahat)
#        so the scripts' `betamat[,1] + betahat[1]` double-counts the fixed effect;
#   (ii) the coefficients are applied to the SEQ matrix, but MM(array, log(seq+1))
#        returns Yhat = the corrected ARRAY. The vignette: "only the values in
#        platform X will be changed".
#
# Variants scored here, all on the paper's own Type_I_error_Power():
#   pub  reference values read from the shipped d_plot (no recomputation)
#   E    vignette-correct, as fitted : MM(array, log-seq) -> Yhat is corrected ARRAY
#   F    direction-consistent        : MM(log-seq, array) -> Yhat is corrected SEQ,
#        which is the direction every other method uses. Needs a zero-variance
#        filter -- ~3,176 genes are constant in seq.
#
# Writes after every (m, seed) so a timeout still leaves usable results.
# ============================================================
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
ROOT <- "/work/users/x/u/xuejun1/Integration_paper_Sim"
OUT  <- file.path(ROOT, "revision_repo/code/mmr_rerun_full_results.csv")
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
  stop("unterminated: ", pat)
}
gl  <- readLines("data_simulation_unbalanced.R", warn = FALSE)
src <- grab(gl, "^sim_array_seq_balanced <- function")
src <- sub('load\\(".*sim_BMI_control0\\.rda"\\)',
           sprintf('load("%s")', file.path(ROOT, "sim_BMI_control0.rda")), src)
eval(parse(text = src))
eval(parse(text = grab(readLines("unbalanced_power_sample_size.R", warn = FALSE),
                       "^Type_I_error_Power <- function")))

owd <- getwd(); setwd(file.path(ROOT, "Cross-Platform-Normalization-master/MatchMixeR/R"))
for (f in c("ckmeans.R","dwd.R","eb.R","functions.R","gq.R","MatchMixeR.R","xpn.R")) source(f)
setwd(owd)

load("d_plot_Power_Sample_Size.rda"); pub <- d_plot

score <- function(da, ds, dn_a, dn_s) {
  da <- as.data.frame(da); ds <- as.data.frame(ds)
  dimnames(da) <- dn_a; dimnames(ds) <- dn_s
  o <- Type_I_error_Power(da, ds)
  c(typeI = o$Type_I_error, power = o$Power_overall)
}

SS <- c(5, 10, 30, 50); SEEDS <- 101:200
rows <- list(); t_start <- Sys.time()

for (m in SS) for (sd_i in SEEDS) {
  r  <- sim_array_seq_balanced(n_per_group = m, n_DE = 1000, seed = sd_i)
  a  <- as.matrix(r$sim_array); s <- as.matrix(r$sim_seq)
  lg <- log(s + 1)
  dn_a <- dimnames(r$sim_array); dn_s <- dimnames(r$sim_seq)

  # ---- E: as fitted; corrected ARRAY vs log-seq ----------------------------
  E <- tryCatch({
    mm <- MM(a, lg)
    c(score(mm$Yhat, lg, dn_a, dn_s),
      sd_ratio = median(apply(mm$Yhat, 1, sd)) / median(apply(lg, 1, sd)),
      lambdahat = mm$lambdahat)
  }, error = function(e) c(typeI = NA, power = NA, sd_ratio = NA, lambdahat = NA))

  # ---- F: direction-consistent; corrected SEQ vs array ---------------------
  F_ <- tryCatch({
    v  <- apply(lg, 1, var); ok <- is.finite(v) & v > 0
    mm2 <- MM(lg[ok, , drop = FALSE], a[ok, , drop = FALSE])
    ds  <- lg                      # genes constant in seq are left as they are
    ds[ok, ] <- mm2$Yhat
    c(score(a, ds, dn_a, dn_s),
      sd_ratio = median(apply(mm2$Yhat, 1, sd)) / median(apply(a, 1, sd)),
      lambdahat = mm2$lambdahat, n_const = sum(!ok))
  }, error = function(e) c(typeI = NA, power = NA, sd_ratio = NA, lambdahat = NA, n_const = NA))

  pv <- pub$Power[pub$method == "MMR" & pub$sample_size == m & pub$seed == sd_i]
  rows[[length(rows) + 1]] <- data.frame(
    m = m, seed = sd_i,
    pub_power   = if (length(pv)) pv else NA_real_,
    E_typeI = E[["typeI"]], E_power = E[["power"]],
    E_sd_ratio = E[["sd_ratio"]], E_lambda = E[["lambdahat"]],
    F_typeI = F_[["typeI"]], F_power = F_[["power"]],
    F_sd_ratio = F_[["sd_ratio"]], F_lambda = F_[["lambdahat"]],
    F_n_const = if ("n_const" %in% names(F_)) F_[["n_const"]] else NA)

  write.csv(do.call(rbind, rows), OUT, row.names = FALSE)   # checkpoint every iteration
  cat(sprintf("m=%-3d seed=%d  pub=%.4f | E: tI=%.4f pw=%.4f | F: tI=%.4f pw=%.4f | %.1f min\n",
              m, sd_i, if (length(pv)) pv else NA_real_,
              E[["typeI"]], E[["power"]], F_[["typeI"]], F_[["power"]],
              as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
  flush.console()
}

tab <- do.call(rbind, rows)
cat("\n=============== MMR, 100 seeds, corrected usage ===============\n")
print(aggregate(cbind(pub_power, E_typeI, E_power, F_typeI, F_power) ~ m, tab,
                FUN = function(x) round(mean(x, na.rm = TRUE), 4)))
cat("\n--- Type I variability across seeds (SD), which 3 seeds could not show ---\n")
print(aggregate(cbind(E_typeI, F_typeI) ~ m, tab,
                FUN = function(x) round(sd(x, na.rm = TRUE), 4)))
cat("\n--- variance shrinkage: Yhat is a fitted value, so check it has not collapsed ---\n")
print(aggregate(cbind(E_sd_ratio, F_sd_ratio) ~ m, tab,
                FUN = function(x) round(mean(x, na.rm = TRUE), 4)))

cat("\n=============== ranking against the published methods ===============\n")
for (m in SS) {
  o <- aggregate(Power ~ method, subset(pub, sample_size == m & seed %in% SEEDS), mean)
  t_m <- tab[tab$m == m, ]
  o <- rbind(o, data.frame(method = "MMR[E corrected]", Power = mean(t_m$E_power, na.rm = TRUE)),
                data.frame(method = "MMR[F direction]", Power = mean(t_m$F_power, na.rm = TRUE)))
  o <- o[order(-o$Power), ]; o$rank <- seq_len(nrow(o))
  cat("\n--- m =", m, "---\n"); print(o, row.names = FALSE, digits = 4)
}
cat("\nNOTE: power is only comparable where Type I is controlled -- read the Type I table first.\n")
