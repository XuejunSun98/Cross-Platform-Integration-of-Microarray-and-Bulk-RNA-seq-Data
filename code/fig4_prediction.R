# Figure 4: cross-platform prediction AUC, three penalties per method.
#
# Layout per the author's design: method on the x axis, THREE boxes per method
# (Lasso / elastic net / ridge) dodged side by side, colour encoding the penalty.
# Methods keep the grouping and dashed separators used in Figures 1 and 3, and
# the uncorrected baseline is shown first as a reference.
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
ROOT <- "/work/users/x/u/xuejun1/Integration_paper_Sim"
RES  <- file.path(ROOT, "revision_repo/results_pred_sweep")
OUT  <- file.path(ROOT, "revision_repo/figures"); dir.create(OUT, showWarnings = FALSE)

d <- bind_rows(lapply(list.files(RES, "\\.csv$", full.names = TRUE), read.csv))
stopifnot(nrow(d) > 0)
cat("rows:", nrow(d), " seeds:", length(unique(d$seed)), " methods:", length(unique(d$method)), "\n")

GROUPS <- list("None"        = "No_correction",
               "Subject-wise"= c("QN","Angel","TDM"),
               "Gene-wise"   = c("MMR","ComBat","ComBat_seq","RNABC","limma","XPN","MNN"))
PEN <- c(lasso = "Lasso", elastic_net = "Elastic net", ridge = "Ridge")
# Penalty gets its OWN palette. The previous one reused the condition colours from Figures 1
# and 3, so green meant "Balanced design" there and "Ridge" here. Purple/crimson/gold is
# disjoint from the condition set and passes all six checks (CVD dE 16.9, normal-vision 27.1).
# The two cool hues sit on Lasso and elastic net, with amber marking ridge.
PAL <- c("Lasso" = "#6f7bd6", "Elastic net" = "#9c5fa8", "Ridge" = "#d99a2b")

present <- unique(d$method)
g   <- lapply(GROUPS, function(m) m[m %in% present]); g <- g[lengths(g) > 0]
lev <- unlist(g, use.names = FALSE)
ends <- cumsum(lengths(g)); n <- lengths(g)
seps <- head(ends, -1) + 0.5
# a label centred on a one-method group at the panel edge spills outside it, so
# the first group is left-aligned and the last right-aligned
labs <- data.frame(grp = names(g), x = ends - n/2 + 0.5, tier = 0,
                   hj = c(0, rep(0.5, length(g) - 1))[seq_along(g)])
labs$x[1] <- 0.6   # only the one-method group at the left edge needs shifting

d$method  <- factor(d$method, levels = lev)
d$penalty <- factor(PEN[d$penalty], levels = PEN)
d$panel   <- factor(paste0(ifelse(d$condition == "balanced", "Balanced", "Matched Imbalance"),
                           "\n", ifelse(d$direction == "array_to_seq",
                                        "train microarray → test RNA-seq",
                                        "train RNA-seq → test microarray")))
d <- d[!is.na(d$method), ]

top <- 1.10; drop <- 0.05
p <- ggplot(d, aes(method, AUC, fill = penalty)) +
  geom_vline(xintercept = seps, linetype = "dashed", colour = "grey45", linewidth = .4) +
  geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey40", linewidth = .4) +
  geom_boxplot(outlier.size = .35, linewidth = .25, width = .72,
               position = position_dodge(width = .78)) +
  geom_label(data = merge(labs, unique(d["panel"]), by = NULL),
             aes(x = x, y = top - tier*drop, label = grp, hjust = hj), inherit.aes = FALSE,
             fontface = "bold", size = 5, vjust = 1,
             fill = "white", label.size = 0, label.padding = unit(0.15, "lines")) +
  facet_wrap(~ panel, nrow = 2) +
  scale_fill_manual(values = PAL, name = "Penalty") +
  coord_cartesian(ylim = c(0.35, top)) +
  labs(title = "Cross-platform outcome prediction", x = "Method", y = "AUC") +
  theme_bw(base_size = 19) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 19),
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 17),
        strip.text = element_text(face = "bold", size = 16),
        legend.position = "top",
        plot.title = element_text(face = "bold", size = 24, hjust = 0.5),
)
ggsave(file.path(OUT, "Fig4_Prediction.png"), p, width = 16, height = 12, dpi = 300)
cat("written:", file.path(OUT, "Fig4_Prediction.png"), "\n")

cat("\n=== mean AUC by method and penalty ===\n")
print(as.data.frame(d %>% group_by(panel, method, penalty) %>%
  summarise(AUC = mean(AUC, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = penalty, values_from = AUC)), row.names = FALSE, digits = 3)
