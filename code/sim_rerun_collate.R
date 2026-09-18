# Bind the per-seed CSVs, check completeness, and report against the shipped results.
#   Rscript sim_rerun_collate.R
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
ROOT <- "/work/users/x/u/xuejun1/Integration_paper_Sim"
RES  <- file.path(ROOT, "revision_repo/results")
SIM  <- file.path(ROOT, "simulation_2024")

f <- list.files(RES, pattern = "_seed[0-9]+\\.csv$", full.names = TRUE)
if (!length(f)) stop("no result files in ", RES)
d <- do.call(rbind, lapply(f, read.csv, stringsAsFactors = FALSE))

cat("files:", length(f), "  rows:", nrow(d), "  methods:", length(unique(d$method)), "\n")

# ---- completeness -----------------------------------------------------------
cat("\n=============== completeness ===============\n")
exp_lv <- c(typeI_balance = 4, power_size = 4, power_balance = 4, power_effect = 4)
for (m in sort(unique(d$method))) {
  sub <- d[d$method == m, ]
  got <- length(unique(sub$seed)); need <- sum(exp_lv) * got
  cat(sprintf("  %-16s seeds=%-4d rows=%-5d %s\n", m, got, nrow(sub),
              if (nrow(sub) == need) "complete" else sprintf("*** expected %d", need)))
}
miss <- setdiff(101:200, unique(d$seed))
if (length(miss)) cat("\n*** MISSING SEEDS:", paste(miss, collapse = ","),
                      "\n    re-run: sbatch --array=", paste(miss, collapse = ","), " <script>\n")

# ---- results ----------------------------------------------------------------
for (g in c("typeI_balance", "power_size", "power_balance", "power_effect")) {
  sub <- d[d$grid == g, ]
  if (!nrow(sub)) next
  col <- if (g == "typeI_balance") "Type_I_error" else "Power"
  cat("\n===============", g, "--", col, "( mean over", length(unique(sub$seed)), "seeds ) ===============\n")
  a <- aggregate(sub[[col]], list(method = sub$method, level = sub$level),
                 function(x) mean(x, na.rm = TRUE))
  w <- reshape(a, idvar = "method", timevar = "level", direction = "wide")
  names(w) <- sub("^x\\.", "", names(w))
  print(w[order(-rowMeans(w[, -1, drop = FALSE], na.rm = TRUE)), ], row.names = FALSE, digits = 4)
  # seed-to-seed spread: Shambhala2 taught us to look at this, not just the mean
  sd_a <- aggregate(sub[[col]], list(method = sub$method), function(x) sd(x, na.rm = TRUE))
  cat("  across-seed SD:", paste(sprintf("%s=%.3f", sd_a$method, sd_a$x), collapse = "  "), "\n")
}

# ---- ranking against the shipped published values ---------------------------
cat("\n=============== power by sample size, new arms vs shipped ===============\n")
e <- new.env(); load(file.path(SIM, "d_plot_Power_Sample_Size.rda"), envir = e)
pub <- e$d_plot
ps <- d[d$grid == "power_size", ]
for (lv in c("5", "10", "30", "50")) {
  o <- aggregate(Power ~ method, subset(pub, sample_size == as.numeric(lv)), mean)
  o$source <- "published"
  n <- aggregate(Power ~ method, subset(ps, level == lv), function(x) mean(x, na.rm = TRUE))
  if (nrow(n)) { n$source <- "rerun"; o <- rbind(o, n) }
  o <- o[order(-o$Power), ]; o$rank <- seq_len(nrow(o))
  cat("\n--- m =", lv, "---\n"); print(o, row.names = FALSE, digits = 4)
}
cat("\nNOTE: Type I in the three power grids uses the POWER scorer (FP / 14160, excluding\n",
    "gene_1..gene_2000). Only typeI_balance uses the Type I scorer (FP / 16160 on null\n",
    "data). The two columns are not comparable across grids.\n")
