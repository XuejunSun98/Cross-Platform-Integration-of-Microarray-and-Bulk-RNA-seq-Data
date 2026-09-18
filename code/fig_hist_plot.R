# Distributional alignment across all four real datasets and all thirteen correction
# methods (plus the uncorrected baseline). Replaces hist_all_real_3datasets_sorted.png,
# which covered three datasets and the original method set only.
#
# Densities are computed upstream: the eleven unsupervised methods by real_hist.R (one
# job), and Shambhala2, Rank-In and COCONUT by real_hist_supervised.R (one array task per
# dataset x method, because they shell out to MATLAB/Python and need class labels).
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(ggh4x)})
ROOT <- "/work/users/x/u/xuejun1/Integration_paper_Sim"
OUT  <- file.path(ROOT, "revision_repo/figures")

d <- readRDS(file.path(ROOT, "revision_repo/results/real_hist_density.rds"))
sup <- list.files(file.path(ROOT, "revision_repo/results_hist_sup"), "\\.rds$", full.names = TRUE)
if (length(sup)) d <- bind_rows(d, do.call(rbind, lapply(sup, readRDS)))

# same order as Table 1 and Supplementary Figures S5/S6, and the class each method belongs
# to is drawn as an outer row strip so the taxonomy is visible in the figure itself
CLASS <- list(
  `Uncorrected`                  = "No correction",
  `Unsupervised\nsubject-wise`   = c("QN", "Angel", "TDM"),
  `Unsupervised\ngene-wise`      = c("MMR", "ComBat", "ComBat-seq", "RNABC", "Shambhala2",
                                     "limma", "XPN", "MNN"),
  `Supervised`                   = c("COCONUT", "Rank-In"))
MORDER <- unlist(CLASS, use.names = FALSE)
class_of <- setNames(rep(names(CLASS), lengths(CLASS)), MORDER)
DORDER <- c("SEQC", "CCLE", "METSIM", "TCGA-LUSC")
miss <- setdiff(MORDER, unique(d$method))
if (length(miss)) cat("NOTE: not yet available:", paste(miss, collapse = ", "), "\n")
d <- d %>% filter(method %in% MORDER, dataset %in% DORDER) %>%
  mutate(method   = factor(method,  levels = MORDER),
         class    = factor(class_of[as.character(method)], levels = names(CLASS)),
         dataset  = factor(dataset, levels = DORDER),
         platform = factor(platform, levels = c("Microarray", "RNA-seq")))

# Rank-In and the uncorrected RNA-seq can carry a long right tail that squashes everything
# else; each panel is trimmed to the 0.1-99.9% mass of its own curves.
d <- d %>% group_by(method, dataset) %>%
  mutate(keep = { cw <- ave(y, platform, FUN = function(v) cumsum(v) / sum(v))
                  cw > 0.0005 & cw < 0.9995 }) %>% ungroup() %>% filter(keep)

PAL <- c("Microarray" = "#2a78d6", "RNA-seq" = "#eb6834")
p <- ggplot(d, aes(x, y, colour = platform, fill = platform)) +
  geom_area(position = "identity", alpha = 0.35, linewidth = 0) +
  geom_line(linewidth = 0.6) +
  facet_nested(class + method ~ dataset, scales = "free", independent = "all", switch = "y",
               strip = strip_nested(background_y = list(
                 element_rect(fill = "grey80", colour = NA),
                 element_rect(fill = "grey94", colour = NA)),
                 by_layer_y = TRUE)) +
  scale_colour_manual(values = PAL, name = NULL) +
  scale_fill_manual(values = PAL, name = NULL) +
  scale_y_continuous(breaks = NULL) +
  labs(x = "Expression value", y = NULL) +
  theme_bw(base_size = 15) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.background = element_rect(fill = "grey93", colour = NA),
        strip.text.y.left = element_text(angle = 0, hjust = 0.5, face = "bold"),
        strip.placement = "outside",
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 10))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(OUT, "Supp4_hist_alignment.png"), p, width = 14, height = 20, dpi = 300, limitsize = FALSE)
cat("wrote figures/Supp4_hist_alignment.png:",
    length(unique(d$method)), "methods x", length(unique(d$dataset)), "datasets\n")
