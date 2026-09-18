# ============================================================
# Is MNN usable on bulk cross-platform data?
#
# Three questions, not just "does it run":
#   1. Does it run at all at the sample sizes used here (m = 5, 10, 30, 50)?
#      mnnCorrect's default k = 20 exceeds the batch size when m = 5.
#   2. How many mutual-nearest-neighbour pairs are actually found? MNN assumes
#      shared discrete populations; a bulk case/control design has one
#      homogeneous cloud, so the structure it keys on may not exist.
#   3. Does it preserve the condition signal, or remove it? MNN assumes the
#      batch effect is SMALL relative to biological variation -- here the
#      platform effect dominates, which is the regime the method warns about.
#
# Scored with the paper's own Type I error / power function so the numbers are
# comparable with the published arms.
# ============================================================
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
ROOT <- "/work/users/x/u/xuejun1/Integration_paper_Sim"
setwd(file.path(ROOT, "simulation_2024"))
suppressPackageStartupMessages({ library(batchelor); library(stringr) })

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

# published values for context
load("d_plot_Power_Sample_Size.rda"); pub <- d_plot

SS <- c(5, 10, 30, 50); SEED <- 101
rows <- list()

for (m in SS) {
  r <- sim_array_seq_balanced(n_per_group = m, n_DE = 1000, seed = SEED)
  a  <- as.matrix(r$sim_array)
  s  <- as.matrix(log(r$sim_seq + 1))      # common scale, as ComBat/limma get
  nb <- ncol(a)                            # samples per batch

  for (k in c(20, max(2, floor(nb / 4)))) {   # default k, and k scaled to batch
    lab <- sprintf("m=%d k=%d (batch n=%d)", m, k, nb)
    t0 <- Sys.time()
    out <- tryCatch({
      res <- batchelor::mnnCorrect(a, s, k = k, cos.norm.in = FALSE, cos.norm.out = FALSE)
      cor_mat <- SummarizedExperiment::assay(res, "corrected")
      pairs   <- S4Vectors::metadata(res)$merge.info$pairs
      npair   <- if (is.null(pairs)) NA_integer_ else sum(sapply(pairs, nrow))
      da <- cor_mat[, 1:nb, drop = FALSE]
      ds <- cor_mat[, (nb + 1):(2 * nb), drop = FALSE]
      sc <- Type_I_error_Power(as.data.frame(da), as.data.frame(ds))
      list(status = "PASS", npair = npair,
           typeI = sc$Type_I_error, power = sc$Power_overall,
           sd_ratio = median(apply(ds, 1, sd)) / median(apply(a, 1, sd)))
    }, error = function(e) list(status = paste("FAIL:", conditionMessage(e)),
                                npair = NA, typeI = NA, power = NA, sd_ratio = NA))
    out$secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    rows[[length(rows) + 1]] <- data.frame(
      m = m, batch_n = nb, k = k, status = substr(out$status, 1, 60),
      mnn_pairs = out$npair, TypeI = out$typeI, Power = out$power,
      sd_ratio = out$sd_ratio, secs = out$secs)
    cat(sprintf("%-28s %-10s pairs=%-6s TypeI=%-8s Power=%-8s %ss\n", lab,
                substr(out$status, 1, 10), out$npair,
                format(out$typeI, digits = 3), format(out$power, digits = 3), out$secs))
    flush.console()
  }
}

tab <- do.call(rbind, rows)
cat("\n================ MNN on bulk cross-platform data ================\n")
print(tab, row.names = FALSE, digits = 4)

cat("\n--- published power for context (seed 101) ---\n")
ctx <- subset(pub, seed == SEED & sample_size %in% SS &
                method %in% c("limma", "Combat", "QN", "MMR"))
print(reshape(ctx[, c("sample_size", "method", "Power")], idvar = "sample_size",
              timevar = "method", direction = "wide"), row.names = FALSE)

write.csv(tab, file.path(ROOT, "mnn_bulk_test_results.csv"), row.names = FALSE)
cat("\nNote: sd_ratio = median per-gene SD of the corrected seq data divided by that of\n",
    "the array data. Values far below 1 indicate the correction has collapsed variation.\n")
