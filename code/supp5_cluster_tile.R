# Supplementary Figure S5: sample-level clustering after correction, rebuilt for all four
# real datasets and all thirteen methods. Replaces cluster_all_real_tile_sorted.png, which
# covered three datasets and the original method set and had no generating script in the
# package.
#
# Each panel is one dataset x method. Samples are columns, ordered by the average-linkage
# dendrogram on Euclidean distance over standardised genes -- the same clustering whose ARI
# is reported in Supplementary Figure S5, computed on the same corrected matrices, so the
# picture and the number cannot disagree. The upper bar is the biological label, the lower
# bar the platform: correction has worked when the lower bar is interleaved (platform no
# longer drives the ordering) while the upper bar stays in blocks.
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(ggh4x)})
ROOT <- "/work/users/x/u/xuejun1/Integration_paper_Sim"
OUT  <- file.path(ROOT, "revision_repo/figures")

f <- list.files(file.path(ROOT, "revision_repo/results_sil_all"), "^tile_.*\\.rds$", full.names = TRUE)
stopifnot(length(f) > 0)
d <- bind_rows(lapply(f, readRDS))

# the class each method belongs to is drawn as an outer row strip, so the taxonomy of
# Table 1 is visible in the figure itself rather than only in the caption
CLASS <- list(
  `Uncorrected`                 = "No_correction",
  `Unsupervised\nsubject-wise`  = c("QN","Angel","TDM"),
  `Unsupervised\ngene-wise`     = c("MMR","ComBat","ComBat_seq","RNABC","Shambhala2",
                                    "limma","XPN","MNN"),
  `Supervised`                  = c("COCONUT","Rank_in"))
MORDER <- unlist(CLASS, use.names = FALSE)
class_of <- setNames(rep(names(CLASS), lengths(CLASS)), MORDER)
MLAB <- c(No_correction="No correction", QN="QN", Angel="Angel", TDM="TDM", MMR="MMR",
          ComBat="ComBat", ComBat_seq="ComBat-seq", RNABC="RNABC", Shambhala2="Shambhala2",
          limma="limma", XPN="XPN", MNN="MNN", COCONUT="COCONUT", Rank_in="Rank-In")
DS <- c("SEQC","CCLE","METSIM","TCGA-LUSC")
miss <- setdiff(MORDER, unique(d$method))
if (length(miss)) cat("NOTE: not available:", paste(miss, collapse = ", "), "\n")

# The biological labels differ by dataset (A/B, breast/large intestine, case/control,
# stage I/II+), so they are recoded to a common two-level scale and named in the caption.
# Both bars then share one fill scale and one legend.
d <- d %>% group_by(dataset) %>%
  mutate(bio2 = ifelse(biological == sort(unique(biological))[1], "Group 1", "Group 2")) %>%
  ungroup()
pl <- bind_rows(
  d %>% transmute(dataset, method, pos, row = "Group",    lab = bio2),
  d %>% transmute(dataset, method, pos, row = "Platform",
                  lab = ifelse(platform == "array", "Microarray", "RNA-seq"))) %>%
  mutate(class   = factor(class_of[as.character(method)], levels = names(CLASS)),
         method  = factor(MLAB[as.character(method)], levels = unname(MLAB[MORDER])),
         dataset = factor(dataset, levels = DS),
         row     = factor(row, levels = c("Platform", "Group")),
         lab     = factor(lab, levels = c("Group 1","Group 2","Microarray","RNA-seq"))) %>%
  filter(!is.na(method))

# Green/magenta for the group bar rather than green/purple: at tile width the old pair was
# hard to separate. Within each bar -- which is the comparison the figure actually asks the
# reader to make -- both pairs are well clear of the floor (group normal dE 38.2, deutan 21.9;
# platform normal 33.6, deutan 24.7). The two bars sit in separate, y-axis-labelled rows, so a
# cross-row pair is never compared directly.
PAL <- c("Group 1" = "#2f9e44", "Group 2" = "#ae3ec9",
         "Microarray" = "#2a78d6", "RNA-seq" = "#eb6834")
b <- 15
p <- ggplot(pl, aes(pos, row, fill = lab)) +
  geom_tile(height = 0.72) +
  # equal-width columns: space = "free_x" would shrink SEQC to a fraction of the others and
  # make its bars unreadable
  facet_nested(class + method ~ dataset, scales = "free_x", switch = "y",
               strip = strip_nested(background_y = list(
                 element_rect(fill = "grey80", colour = NA),
                 element_rect(fill = "grey94", colour = NA)),
                 by_layer_y = TRUE)) +
  scale_fill_manual(values = PAL, name = NULL) +
  scale_x_continuous(expand = c(0, 0)) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 4))) +
  labs(x = "Samples, ordered by hierarchical clustering of the corrected data", y = NULL) +
  theme_bw(base_size = b) +
  theme(legend.position = "top",
        panel.grid = element_blank(), panel.spacing.y = unit(0.28, "lines"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = b * 0.7),
        strip.background = element_rect(fill = "grey93", colour = NA),
        strip.text.x = element_text(face = "bold", size = b),
        strip.text.y.left = element_text(angle = 0, hjust = 0.5, face = "bold", size = b * 0.9),
        strip.placement = "outside")
ggsave(file.path(OUT, "Supp5_cluster_tile.png"), p, width = 18, height = 17, dpi = 300, limitsize = FALSE)
cat("wrote Supp5_cluster_tile.png:", nlevels(droplevels(pl$method)), "methods x",
    nlevels(droplevels(pl$dataset)), "datasets\n")
