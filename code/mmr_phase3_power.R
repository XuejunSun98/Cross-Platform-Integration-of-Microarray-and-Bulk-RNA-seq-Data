# ============================================================
# Phase 3: does the MMR implementation question change POWER,
# and does it move MMR's rank against the other methods?
#
# Figure 3 path: data_simulation_unbalanced.R::sim_array_seq_balanced
#   (n_per_group in {5,10,30,50}, n_DE = 1000, seeds 101:200)
#   -> sim_balanced_power_by_sample_size.rda
#   -> unbalanced_power_sample_size.R -> d_plot_Power_Sample_Size.rda
#
# Other methods' published values are read straight from the shipped
# d_plot, so only MMR is recomputed.
# ============================================================
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
ROOT <- "/work/users/x/u/xuejun1/Integration_paper_Sim"
setwd(file.path(ROOT, "simulation_2024"))
suppressPackageStartupMessages(library(stringr))

# ---- extract a function verbatim, by brace balance ----------------------
grab <- function(lines, pattern) {
  s <- grep(pattern, lines)[1]
  bal <- 0; opened <- FALSE
  for (k in s:length(lines)) {
    txt <- gsub("#.*$", "", lines[k])
    no <- lengths(regmatches(txt, gregexpr("\\{", txt)))
    nc <- lengths(regmatches(txt, gregexpr("\\}", txt)))
    if (no > 0) opened <- TRUE
    bal <- bal + no - nc
    # only stop once the body's opening brace has actually been seen --
    # a multi-line function signature keeps bal at 0 for several lines
    if (opened && bal == 0) return(paste(lines[s:k], collapse = "\n"))
  }
  stop("unterminated: ", pattern)
}

gl <- readLines("data_simulation_unbalanced.R", warn = FALSE)
src <- grab(gl, "^sim_array_seq_balanced <- function")
src <- sub('load\\(".*sim_BMI_control0\\.rda"\\)',
           sprintf('load("%s")', file.path(ROOT, "sim_BMI_control0.rda")), src)
eval(parse(text = src))

cl <- readLines("unbalanced_power_sample_size.R", warn = FALSE)
eval(parse(text = grab(cl, "^Type_I_error_Power <- function")))

owd <- getwd(); setwd(file.path(ROOT, "Cross-Platform-Normalization-master/MatchMixeR/R"))
for (f in c("ckmeans.R","dwd.R","eb.R","functions.R","gq.R","MatchMixeR.R","xpn.R")) source(f)
setwd(owd)

# ---- MMR variants -------------------------------------------------------
B_published <- function(a, s) {
  mm <- MM(as.matrix(a), as.matrix(log(s + 1)))
  b0 <- mm$betamat[,1] + mm$betahat[1]; b1 <- mm$betamat[,2] + mm$betahat[2]
  sweep(sweep(s, 1, b1, "*"), 1, b0, "+")
}
C_scalefix <- function(a, s) {
  mm <- MM(as.matrix(a), as.matrix(log(s + 1)))
  b0 <- mm$betamat[,1] + mm$betahat[1]; b1 <- mm$betamat[,2] + mm$betahat[2]
  sweep(sweep(log(s + 1), 1, b1, "*"), 1, b0, "+")
}
D2_direction <- function(a, s) {
  lg <- log(s + 1); v <- apply(lg, 1, var); ok <- is.finite(v) & v > 0
  mm <- MM(as.matrix(lg[ok, , drop = FALSE]), as.matrix(a[ok, , drop = FALSE]))
  b0 <- mm$betamat[,1] + mm$betahat[1]; b1 <- mm$betamat[,2] + mm$betahat[2]
  out <- lg
  out[ok, ]  <- sweep(sweep(lg[ok, , drop = FALSE], 1, b1, "*"), 1, b0, "+")
  out[!ok, ] <- lg[!ok, , drop = FALSE] * mm$betahat[2] + mm$betahat[1]
  out
}
VAR <- list(B_published = B_published, C_scalefix = C_scalefix, D2_direction = D2_direction)

load("d_plot_Power_Sample_Size.rda")
pub <- d_plot

SS <- c(5, 10, 30, 50); SEEDS <- 101:103
rows <- list()
for (n in SS) for (sd in SEEDS) {
  r <- sim_array_seq_balanced(n_per_group = n, n_DE = 1000, seed = sd)
  a <- r$sim_array; s <- r$sim_seq
  row <- data.frame(sample_size = n, seed = sd,
                    published = pub$Power[pub$method == "MMR" & pub$sample_size == n & pub$seed == sd])
  for (nm in names(VAR)) {
    row[[nm]] <- tryCatch({
      ds <- as.data.frame(VAR[[nm]](a, s)); dimnames(ds) <- dimnames(s)
      out <- Type_I_error_Power(as.data.frame(a), ds)
      if (!is.null(out$Power_overall)) out$Power_overall else NA_real_
    }, error = function(e) NA_real_)
  }
  rows[[length(rows)+1]] <- row
  cat(sprintf("n=%-3d seed=%d published=%.4f B=%.4f C=%.4f D2=%.4f\n",
              n, sd, row$published, row$B_published, row$C_scalefix, row$D2_direction))
  flush.console()
}
tab <- do.call(rbind, rows)
write.csv(tab, file.path(ROOT, "mmr_phase3_power.csv"), row.names = FALSE)

cat("\n=============== MMR power: published vs variants ===============\n")
print(tab, row.names = FALSE, digits = 5)
cat("\n--- does B reproduce the published power? ---\n")
cat("max |B - published| =", format(max(abs(tab$B_published - tab$published), na.rm = TRUE)), "\n")

cat("\n--- mean power by sample size ---\n")
print(aggregate(cbind(published, B_published, C_scalefix, D2_direction) ~ sample_size, tab,
                FUN = function(x) round(mean(x), 4)))

# ---- ranking: MMR variants against the other methods (published) --------
cat("\n=============== ranking at each sample size ===============\n")
others <- subset(pub, method != "MMR" & seed %in% SEEDS & sample_size %in% SS)
om <- aggregate(Power ~ method + sample_size, others, mean)
for (n in SS) {
  o <- om[om$sample_size == n, c("method","Power")]
  mine <- tab[tab$sample_size == n, ]
  for (nm in c("published", names(VAR)))
    o <- rbind(o, data.frame(method = paste0("MMR[", nm, "]"), Power = mean(mine[[nm]], na.rm = TRUE)))
  o <- o[order(-o$Power), ]
  o$rank <- seq_len(nrow(o))
  cat("\n--- n =", n, "---\n"); print(o, row.names = FALSE, digits = 4)
}
