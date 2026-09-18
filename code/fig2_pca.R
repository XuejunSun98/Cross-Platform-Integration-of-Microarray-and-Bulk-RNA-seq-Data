# Figure 2 (R3.m3): PCA of samples under the NULL setting, before and after
# correction, for BOTH supervised methods.
#
# The published figure has two rows -- uncorrected and Rank-In -- so the
# conclusion about "supervised methods" rested on one method. This adds a third
# row for COCONUT, and corrects the caption: the panels are a principal
# component analysis, not a clustering result.
#
# Data: the four balance designs with n_DE = 0 (no truly DE genes), m = 100 per
# platform, matching the 100 array + 100 seq samples of the published figure.
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
source("/work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/code/sim_rerun_common.R")
source(file.path(ROOT, "longleaf_external_tools.R"))
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(COCONUT)})
OUT <- file.path(ROOT, "revision_repo/figures"); dir.create(OUT, showWarnings = FALSE)
WD  <- file.path(ROOT, "revision_repo/scratch_fig2"); dir.create(WD, showWarnings = FALSE)

DESIGNS <- c(balanced = "Balanced", same_imbalance = "Matched Imbalance",
             opp_imbalance = "Reversed Imbalance", half_imbalance = "Mixed Balance-Imbalance")
SEED <- 101

grp_of <- function(cn) ifelse(grepl("control", cn, ignore.case = TRUE), "Control", "Case")

pca_df <- function(a, s, design, row) {
  d <- cbind(as.matrix(a), as.matrix(s))
  d <- d[apply(d, 1, function(x) var(x) > 0 & all(is.finite(x))), , drop = FALSE]
  pc <- prcomp(t(d), center = TRUE, scale. = FALSE)
  ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  data.frame(PC1 = pc$x[,1], PC2 = pc$x[,2],
             group = grp_of(colnames(d)),
             platform = c(rep("Array", ncol(a)), rep("Seq", ncol(s))),
             design = DESIGNS[[design]], row = row,
             lab = sprintf("PC1 (%.1f%%)", ve[1]), lab2 = sprintf("PC2 (%.1f%%)", ve[2]))
}

# ---- COCONUT: controls coded 0, applied per platform ------------------------
coconut_correct <- function(a, s) {
  g <- as.factor(ifelse(grepl("control", colnames(a), ignore.case = TRUE), 0, 1))
  gs <- as.factor(ifelse(grepl("control", colnames(s), ignore.case = TRUE), 0, 1))
  la <- list(pheno = data.frame(group = g,  platform = "array", row.names = colnames(a)),
             genes = as.data.frame(a))
  ls_ <- list(pheno = data.frame(group = gs, platform = "seq",  row.names = colnames(s)),
              genes = as.data.frame(log(s + 1)))
  out <- COCONUT(GSEs = list(d_array = la, d_seq = ls_),
                 control.0.col = "group", byPlatform = FALSE)
  list(array = cbind(out$controlList$GSEs$d_array$genes, out$COCONUTList$d_array$genes),
       seq   = cbind(out$controlList$GSEs$d_seq$genes,   out$COCONUTList$d_seq$genes))
}

# ---- Rank-In ----------------------------------------------------------------
rankin_correct <- function(a, s) {
  old <- setwd(WD); on.exit(setwd(old))
  # exact input format used by Unbalanced_Type_I_RankIn.R: a leading "gene"
  # column, no row names, and RAW seq counts (not logged)
  d_all  <- cbind(gene = rownames(a), as.data.frame(a), as.data.frame(s))
  sample <- data.frame(sample_name = colnames(d_all)[-1])
  sample$Class <- ifelse(grepl("^case", sample$sample_name), 1, 0)
  write.table(d_all,  "d_all_Rankin.txt", row.names = FALSE, sep = "\t", quote = FALSE)
  write.table(sample, "sample.txt",       row.names = FALSE, sep = "\t", quote = FALSE)
  RankIn(quiet = FALSE)
  r <- read.table("Rank_In_r_result.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
  list(array = as.matrix(r[, 1:ncol(a), drop = FALSE]),
       seq   = as.matrix(r[, ncol(a) + seq_len(ncol(s)), drop = FALSE]))
}

rows <- list()
for (dn in names(DESIGNS)) {
  r <- sim_cond_power(m = 100, n_case_bal = 50, n_case_imb = 64,
                      condition = dn, n_DE = 0, seed = SEED)
  a <- as.matrix(r$sim_array); s <- as.matrix(r$sim_seq)
  cat(dn, ": array", ncol(a), " seq", ncol(s), "\n"); flush.console()
  rows[[length(rows)+1]] <- pca_df(a, log(s + 1), dn, "Before correction")
  ck <- coconut_correct(a, s)
  rows[[length(rows)+1]] <- pca_df(ck$array, ck$seq, dn, "After COCONUT")
  rk <- tryCatch(rankin_correct(a, s), error = function(e) { cat("  Rank-In ERR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(rk)) rows[[length(rows)+1]] <- pca_df(rk$array, rk$seq, dn, "After Rank-In")
  cat("  done\n"); flush.console()
}
d <- bind_rows(rows)
d$design <- factor(d$design, levels = DESIGNS)
d$row    <- factor(d$row, levels = c("Before correction", "After Rank-In", "After COCONUT"))

p <- ggplot(d, aes(PC1, PC2, colour = group, shape = platform)) +
  geom_point(size = 2.1, alpha = .85) +
  facet_grid(row ~ design, scales = "free") +
  scale_colour_manual(values = c(Case = "#eb6834", Control = "#2a78d6"), name = "group") +
  scale_shape_manual(values = c(Array = 16, Seq = 17), name = "platform") +
  labs(title = "Principal component analysis of samples under the null setting",
       subtitle = "No truly differentially expressed genes; 100 samples per platform") +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 16),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 17),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        plot.title = element_text(face = "bold", size = 22),
        plot.subtitle = element_text(size = 15))
ggsave(file.path(OUT, "Fig2_PCA_supervised.png"), p, width = 17, height = 12, dpi = 300)
cat("\nwritten:", file.path(OUT, "Fig2_PCA_supervised.png"), "\n")

cat("\n=== group separation on PC1 (|standardised mean difference|) ===\n")
print(d %>% group_by(row, design) %>%
  summarise(sep = abs(diff(tapply(PC1, group, mean))) / sd(PC1), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = design, values_from = sep) %>% as.data.frame(),
  row.names = FALSE, digits = 3)
