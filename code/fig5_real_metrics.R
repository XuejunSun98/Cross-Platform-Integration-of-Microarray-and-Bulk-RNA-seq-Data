# Main Figure 5 (R2.8): quantitative real-data metrics, all thirteen methods.
#
# Form: four magnitudes across method identity, so a dot plot faceted by dataset. Platform
# silhouette, biological silhouette and ARI are DIFFERENT measures and each gets its own
# row -- never two y-scales on one panel. The uncorrected value is drawn as a reference
# line per panel, so "does correction help?" is read directly rather than inferred from a
# legend entry. Methods are ordered and separated by the taxonomy of Table 1.
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
ROOT <- "/work/users/x/u/xuejun1/Integration_paper_Sim"
OUT  <- file.path(ROOT, "revision_repo/figures")

f <- list.files(file.path(ROOT, "revision_repo/results_sil_all"), "\\.csv$", full.names = TRUE)
stopifnot(length(f) > 0)
d <- bind_rows(lapply(f, read.csv))
write.csv(d, file.path(ROOT, "revision_repo/results/real_silhouette_all.csv"), row.names = FALSE)

GROUPS <- list("Uncorrected"  = "No_correction",
               "Subject-wise" = c("QN","Angel","TDM"),
               "Gene-wise"    = c("MMR","ComBat","ComBat_seq","RNABC","Shambhala2","limma","XPN","MNN"),
               "Supervised"   = c("COCONUT","Rank_in"))
LAB <- c(No_correction="No correction", QN="QN", Angel="Angel", TDM="TDM", MMR="MMR", ComBat="ComBat",
         ComBat_seq="ComBat-seq", RNABC="RNABC", Shambhala2="Shambhala2", limma="limma",
         XPN="XPN", MNN="MNN", COCONUT="COCONUT", Rank_in="Rank-In")
LEV  <- unlist(GROUPS, use.names = FALSE)
ends <- cumsum(lengths(GROUPS)); seps <- head(ends, -1) + 0.5
DS   <- c("SEQC","CCLE","METSIM","TCGA-LUSC")
MET  <- c("sil_platform","sil_biological","ARI")
MLAB <- c("Platform silhouette\n(0 = platforms intermixed)",
          "Biological silhouette\n(larger = groups preserved)",
          "Adjusted Rand index\n(larger = clusters match groups)")

long <- function(x) x %>% select(dataset, method, all_of(MET)) %>%
  pivot_longer(all_of(MET), names_to = "metric", values_to = "value") %>%
  mutate(dataset = factor(dataset, levels = DS),
         metric  = factor(metric, levels = MET, labels = MLAB))

pl <- d %>% long() %>%
  mutate(method = factor(LAB[as.character(method)], levels = unname(LAB[LEV]))) %>%
  filter(!is.na(method))
na_pts <- pl %>% filter(is.na(value))
if (nrow(na_pts)) cat("not computable:", nrow(na_pts), "cells\n")

b <- 17
p <- ggplot(pl, aes(method, value)) +
  geom_vline(xintercept = seps, linetype = "dotted", colour = "grey60", linewidth = .35) +
  geom_point(aes(colour = method == "No correction"), size = 3.1, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "#d6336c", `FALSE` = "#1baf7a")) +
  facet_grid(metric ~ dataset, scales = "free_y", switch = "y") +
  labs(title = "Quantitative evaluation of the real-data comparison",
       subtitle = paste("Leftmost point in each panel is the uncorrected data.",
                        "Methods are grouped as in Table 1: subject-wise | gene-wise | supervised."),
       x = "Method", y = NULL) +
  theme_bw(base_size = b) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = b * 0.8),
        strip.text.x = element_text(face = "bold", size = b),
        strip.text.y.left = element_text(size = b * 0.75, angle = 90),
        strip.placement = "outside", strip.background.y = element_blank(),
        plot.title = element_text(face = "bold", size = b * 1.35, hjust = 0.5),
        plot.subtitle = element_text(size = b * 0.78, hjust = 0.5))
ggsave(file.path(OUT, "Fig5_real_metrics.png"), p, width = 18, height = 12.5, dpi = 300)
cat("written Fig5_real_metrics.png\n")
