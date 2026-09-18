# Supplementary Figures S2 and S3.
#   S2 (R2.5): does within-platform quantile normalization, rather than the batch-correction
#              algorithm, drive the Type I error and power of ComBat and limma?
#   S3 (R2.7): the joint-regression baseline -- platform as a covariate in one model
#              (limma lmFit + eBayes, no corrected matrix) against correct-then-test
#              (removeBatchEffect + Wilcoxon).
#
# Uses the manuscript's own grids and scorers, so Type I error keeps its definition (the
# null grid has n_DE = 0) and power is the usual power. Only the QN step differs between
# arms. Source: results_noqn_grids/ (job 1237025, 100 seeds).
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(patchwork)})
ROOT <- "/work/users/x/u/xuejun1/Integration_paper_Sim"
OUT  <- file.path(ROOT, "revision_repo/figures")
d <- bind_rows(lapply(list.files(file.path(ROOT, "revision_repo/results_noqn_grids"),
                                 "\\.csv$", full.names = TRUE), read.csv))
COND <- c(balanced = "Balanced", same_imbalance = "Matched Imbalance",
          opp_imbalance = "Reversed Imbalance", half_imbalance = "Mixed Balance-Imbalance")
lab <- c(ComBat = "ComBat", ComBat_noQN = "ComBat, no QN",
         limma = "limma",  limma_noQN  = "limma, no QN")
PAL <- c("ComBat" = "#1baf7a", "ComBat, no QN" = "#a8e6cf",
         "limma"  = "#eb6834", "limma, no QN"  = "#f7c1a8")
d$arm <- factor(lab[d$method], levels = lab)

base <- 19
mk <- function(dd, xlab_, ylab_, yvar, ttl) {
  ggplot(dd, aes(level, .data[[yvar]], fill = arm)) +
    { if (yvar == "Type_I_error")
        geom_hline(yintercept = 0.05, linetype = "dashed", colour = "red", linewidth = .5) } +
    geom_boxplot(outlier.size = .4, linewidth = .3, width = .7,
                 position = position_dodge(width = .8)) +
    scale_fill_manual(values = PAL, name = NULL) +
    labs(title = ttl, x = xlab_, y = ylab_) +
    theme_bw(base_size = base) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 20, hjust = 1, size = base),
          axis.title = element_text(size = base * 1.1),
          legend.position = "top", legend.text = element_text(size = base),
          strip.text = element_text(face = "bold", size = base * 0.95),
          plot.title = element_text(face = "bold", size = base * 1.2, hjust = 0.5))
}
ti <- d %>% filter(grid == "typeI_balance") %>% mutate(level = factor(COND[level], levels = COND))
p1 <- mk(ti, "Balance design", "Type I error", "Type_I_error",
         "Type I error under the null hypothesis")
LVMAP <- c("5"="20","10"="40","30"="120","50"="200")
pw <- d %>% filter(grid != "typeI_balance") %>%
  mutate(panel = recode(grid, power_size = "Total sample size",
                        power_effect = "Effect size", power_balance = "Balance design"),
         panel = factor(panel, levels = c("Total sample size","Effect size","Balance design")),
         level = ifelse(grid == "power_size", LVMAP[level],
                 ifelse(grid == "power_balance", COND[level], level)),
         level = factor(level, levels = c("20","40","120","200","0.3","0.5","0.7","1", COND)))
p2 <- mk(pw, NULL, "Power", "Power", "Power under the alternative hypothesis") +
  facet_wrap(~ panel, nrow = 1, scales = "free_x")
ggsave(file.path(OUT, "Supp2_QN_effect.png"), (p1 / p2) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "top"), width = 17, height = 12, dpi = 300)
cat("written Supp2_QN_effect.png\n")

# ---- Supplementary Figure S3: one-step vs two-step limma (R2.7) --------------
d3 <- bind_rows(lapply(list.files(file.path(ROOT, "revision_repo/results_limma_1step"),
                                  "\\.csv$", full.names = TRUE), read.csv))
lab3 <- c(limma_2step = "Correct, then test (removeBatchEffect + Wilcoxon)",
          limma_1step = "Platform as covariate (lmFit + eBayes)")
PAL3 <- setNames(c("#eb6834", "#2a78d6"), lab3)
d3$arm <- factor(lab3[d3$method], levels = lab3)
mk3 <- function(dd, xlab_, ylab_, yvar, ttl) {
  ggplot(dd, aes(level, .data[[yvar]], fill = arm)) +
    { if (yvar == "Type_I_error")
        geom_hline(yintercept = 0.05, linetype = "dashed", colour = "red", linewidth = .5) } +
    geom_boxplot(outlier.size = .4, linewidth = .3, width = .62,
                 position = position_dodge(width = .72)) +
    scale_fill_manual(values = PAL3, name = NULL) +
    labs(title = ttl, x = xlab_, y = ylab_) +
    theme_bw(base_size = base) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 20, hjust = 1, size = base),
          axis.title = element_text(size = base * 1.1),
          legend.position = "top", legend.text = element_text(size = base * 0.9),
          strip.text = element_text(face = "bold", size = base * 0.95),
          plot.title = element_text(face = "bold", size = base * 1.2, hjust = 0.5))
}
ti3 <- d3 %>% filter(grid == "typeI_balance") %>% mutate(level = factor(COND[level], levels = COND))
q1 <- mk3(ti3, "Balance design", "Type I error", "Type_I_error",
          "Type I error under the null hypothesis")
pw3 <- d3 %>% filter(grid != "typeI_balance") %>%
  mutate(panel = recode(grid, power_size = "Total sample size",
                        power_effect = "Effect size", power_balance = "Balance design"),
         panel = factor(panel, levels = c("Total sample size","Effect size","Balance design")),
         level = ifelse(grid == "power_size", LVMAP[level],
                 ifelse(grid == "power_balance", COND[level], level)),
         level = factor(level, levels = c("20","40","120","200","0.3","0.5","0.7","1", COND)))
q2 <- mk3(pw3, NULL, "Power", "Power", "Power under the alternative hypothesis") +
  facet_wrap(~ panel, nrow = 1, scales = "free_x")
ggsave(file.path(OUT, "Supp3_joint_model.png"), (q1 / q2) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "top"), width = 17, height = 12, dpi = 300)
cat("written Supp3_joint_model.png\n")
