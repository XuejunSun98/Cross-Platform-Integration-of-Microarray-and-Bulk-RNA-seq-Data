# Type I error and power figures from the 100-seed re-run, in the manuscript's
# own layout: methods ordered and GROUPED by category, dashed separators between
# groups, group labels inside each panel, fill by condition.
#
# Categories follow Figure_Power_final_typeI_colors.png. Note MatchMixeR is
# GENE-wise (it fits a per-gene regression), not subject-wise.
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(patchwork); library(grid)})
ROOT <- "/work/users/x/u/xuejun1/Integration_paper_Sim"
RES  <- file.path(ROOT, "revision_repo/results"); SIM <- file.path(ROOT, "simulation_2024")
OUT  <- file.path(ROOT, "revision_repo/figures"); dir.create(OUT, showWarnings = FALSE)

new <- bind_rows(lapply(list.files(RES, "\\.csv$", full.names = TRUE), read.csv))

# NOTE: results_sbl is NOT read here. Those runs double-logged the RNA-seq (log on input
# AND on Shambhala2's output, which is already on the array scale), so they are invalid.
# Re-enable once the corrected runs land.
shipped <- function(f) { e <- new.env(); load(file.path(SIM, f), envir = e); e$d_plot }
RERUN <- c("MMR","MNN","ComBat_seq","ComBat","limma","RNABC","XPN")
carry <- function(d, val) d %>% mutate(method = recode(as.character(method), Combat = "ComBat")) %>%
  filter(!method %in% RERUN) %>% rename(value = !!val)

GROUPS <- list(
  "Subject-wise" = c("QN","Angel","TDM"),
  "Gene-wise"    = c("MMR","ComBat","ComBat_seq","RNABC","Shambhala2","limma","XPN","MNN"),
  "Supervised"   = c("COCONUT","Rank_in"),
  "Meta"         = c("Meta"))
COND <- c(balanced = "Balanced", same_imbalance = "Matched Imbalance",
          opp_imbalance = "Reversed Imbalance", half_imbalance = "Mixed Balance-Imbalance")
PAL  <- c("Balanced"="#1baf7a","Matched Imbalance"="#eb6834",
          "Reversed Imbalance"="#2a78d6","Mixed Balance-Imbalance"="#e87ba4")  # validated

# order methods by group, keeping only those present; return levels + separators + labels
layout_of <- function(present) {
  g <- lapply(GROUPS, function(m) m[m %in% present])
  g <- g[lengths(g) > 0]
  lev <- unlist(g, use.names = FALSE)
  n   <- lengths(g); ends <- cumsum(n)
  # stagger alternate labels vertically: "Supervised" (2 slots) and "Meta" (1 slot)
  # are wider than their groups and collide if placed at the same height
  list(levels = lev,
       seps   = head(ends, -1) + 0.5,
       labs   = data.frame(grp = names(g), x = ends - n/2 + 0.5,
                           tier = seq_along(g) %% 2))
}

grouped_panel <- function(d, facet, nrow, title, sub, ylab, hline = NA, base = 13,
                          ytop = NULL, ybot = NULL) {
  L <- layout_of(unique(d$method))
  d$method <- factor(d$method, levels = L$levels)
  d <- d[!is.na(d$method), ]
  top  <- if (is.null(ytop)) max(d$value, na.rm = TRUE) * 1.18 else ytop
  drop <- 0.055 * diff(range(d$value, na.rm = TRUE))
  ggplot(d, aes(method, value, fill = cond)) +
    geom_vline(xintercept = L$seps, linetype = "dashed", colour = "grey45", linewidth = .4) +
    { if (!is.na(hline)) geom_hline(yintercept = hline, linetype = "dashed",
                                    colour = "red", linewidth = .45) } +
    geom_boxplot(outlier.size = .5, linewidth = .3, width = .68) +
    geom_text(data = merge(L$labs, unique(d[facet]), by = NULL),
              aes(x = x, y = top - tier * drop, label = grp), inherit.aes = FALSE,
              fontface = "bold", size = base * 0.25, vjust = 1) +
    facet_wrap(as.formula(paste("~", facet)), nrow = nrow) +
    scale_fill_manual(values = PAL, name = "Condition") +
    expand_limits(y = top) +
    # a restricted lower bound is applied with coord_cartesian so nothing is dropped --
    # ybot is always set below the data minimum, so no box or outlier is clipped
    { if (!is.null(ybot)) coord_cartesian(ylim = c(ybot, top)) } +
    labs(title = title, subtitle = sub, x = "Method", y = ylab) +
    theme_bw(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = base * 1.0),
          axis.text.y = element_text(size = base * 1.0),
          axis.title  = element_text(size = base * 1.15),
          strip.text = element_text(face = "bold", size = base),
          # a single condition needs no legend box -- the subtitle names it
          legend.position = if (nlevels(droplevels(factor(d$cond))) > 1) "top" else "none",
          plot.title = element_text(face = "bold", size = base * 1.3, hjust = 0.5))
}

# ================= Type I : 2 x 2 =================
ti <- bind_rows(
  new %>% filter(grid == "typeI_balance") %>% select(method, level, value = Type_I_error),
  carry(shipped("d_plot_Type_I_error_small.rda"), "Type_I_error") %>%
    select(method, level = balance, value) %>% mutate(level = as.character(level)))
ti$cond  <- factor(COND[ti$level], levels = COND)
ti$level <- ti$cond
# R3.m2 asked for a grouped layout comparing all methods in the same panel. One panel,
# methods on x, the four designs as dodged fills -- the same structure as Figure 4.
L1 <- layout_of(unique(ti$method))
ti$method <- factor(ti$method, levels = L1$levels)
ti <- ti[!is.na(ti$method), ]
ti$cond <- factor(ti$cond, levels = COND)
# The axis spans the full range so Rank-In's Reversed (median 0.688) and Mixed (0.547)
# boxes are shown in full rather than clipped.
base1 <- 24
top1  <- max(ti$value, na.rm = TRUE) * 1.06
p1fig <- ggplot(ti, aes(method, value, fill = cond)) +
  geom_vline(xintercept = L1$seps, linetype = "dashed", colour = "grey45", linewidth = .4) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "red", linewidth = .5) +
  geom_boxplot(outlier.size = .4, linewidth = .3, width = .74,
               position = position_dodge(width = .8)) +
  geom_text(data = L1$labs, aes(x = x, y = top1 * 0.995, label = grp), inherit.aes = FALSE,
            fontface = "bold", size = base1 * 0.26, vjust = 1) +
  scale_fill_manual(values = PAL, name = "Design") +
  expand_limits(y = top1) +
  labs(title = "Type I error under the null hypothesis", x = "Method", y = "Type I error") +
  theme_bw(base_size = base1) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = base1 * 1.0),
        axis.text.y = element_text(size = base1 * 1.0),
        axis.title  = element_text(size = base1 * 1.15),
        legend.position = "top", legend.text = element_text(size = base1),
        plot.title = element_text(face = "bold", size = base1 * 1.4, hjust = 0.5))
ggsave(file.path(OUT, "Fig1_TypeI.png"), p1fig, width = 21, height = 11, dpi = 300)

# ================= Power : 3 x 4, as in the manuscript figure =================
pw <- function(g, ship, lvcol) bind_rows(
  new %>% filter(grid == g) %>% select(method, level, value = Power) %>%
    mutate(level = as.character(level)),
  carry(shipped(ship), "Power") %>% select(method, level = !!lvcol, value) %>%
    mutate(level = as.character(level)))

p1 <- pw("power_size", "d_plot_Power_Sample_Size.rda", "sample_size") %>%
  mutate(cond = "Balanced",
         level = factor(level, c("5","10","30","50"), labels = c("20","40","120","200")),
         row = "Sample size")
p2 <- pw("power_effect", "d_plot_Power_Effect_Size.rda", "effect_size") %>%
  mutate(cond = "Balanced",
         level = factor(level, c("0.3","0.5","0.7","1")), row = "Effect size")
p3 <- pw("power_balance", "d_plot_Power_Imbalance.rda", "balance") %>%
  mutate(cond = COND[level], level = factor(COND[level], levels = COND), row = "Balance")

# R1.3: total number of DE genes, at a fixed 1:1 up:down split. Same design as the
# Balance panel (50 per group, 200 samples). Ten methods only -- COCONUT and Rank-In
# use outcome labels, and Shambhala2 and Meta need an external call per replicate.
p4 <- bind_rows(lapply(list.files(file.path(ROOT, "revision_repo/results_de_sweep"),
                                  "\\.csv$", full.names = TRUE), read.csv)) %>%
  filter(ratio == "1:1") %>%
  transmute(method, value = Power, cond = "Balanced",
            level = factor(total, levels = c(500, 1000, 2000, 4000)), row = "DE count")

# R1.3: up:down ratio at a fixed total of 2,000 DE genes. Same design as the Balance
# panel (50 per group, 200 samples). 1:0 is the extreme all-upregulated case.
p5 <- bind_rows(lapply(list.files(file.path(ROOT, "revision_repo/results_de_ratio"),
                                  "\\.csv$", full.names = TRUE), read.csv)) %>%
  # four levels rather than six: balanced, moderate, strong, and the all-up extreme.
  # 2:1 and 3:1 sit on the same monotone trend as 4:1 and add no shape.
  filter(ratio %in% c("1:1","4:1","9:1","1:0")) %>%
  transmute(method, value = power_raw, cond = "Balanced",
            level = factor(ratio, levels = c("1:1","4:1","9:1","1:0")),
            row = "DE ratio")

# One 3 x 4 figure, as in the manuscript: row 1 varies total sample size, row 2 the
# effect size, row 3 the label-balance design. Facet strips name the varied quantity so
# the rows cannot be read as a single sweep. Fill still encodes condition, which is why
# only row 3 is multi-coloured.
LV <- list(
  paste0("Total n = ", c("20","40","120","200")),
  paste0("Effect size = ", c("0.3","0.5","0.7","1")),
  COND)
fig3 <- bind_rows(
  p1 %>% mutate(level = paste0("Total n = ", as.character(level))),
  p2 %>% mutate(level = paste0("Effect size = ", as.character(level))),
  p3 %>% mutate(level = as.character(level)))
fig3$level <- factor(fig3$level, levels = unlist(LV))
fig3$cond  <- factor(fig3$cond, levels = COND)
stopifnot(!any(is.na(fig3$level)))
# Each row gets its OWN y range. On a shared 0-1 axis the balance row (all values
# 0.5-0.95) wastes half its height, which compresses exactly the differences the row
# exists to show. Rows are stacked with patchwork so within-row comparisons are still
# on one scale; only across-row comparisons need the axis read.
mk_row <- function(d, lv, showx, legend) {
  d <- d[d$level %in% lv, ]
  d$level <- factor(as.character(d$level), levels = lv)
  d$cond  <- factor(d$cond, levels = COND)
  rng <- range(d$value, na.rm = TRUE)
  yb  <- rng[1] - 0.02                       # strictly below the minimum: nothing clipped
  yt  <- rng[2] + 0.16 * diff(rng)           # headroom for the in-panel group labels
  p <- grouped_panel(d, "level", nrow = 1, NULL, NULL, "Power",
                     base = 21, ytop = yt, ybot = yb)
  p <- p + theme(legend.position = if (legend) "bottom" else "none")
  if (!showx) p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  p
}
fig3p <- mk_row(fig3, LV[[1]], FALSE, FALSE) /
         mk_row(fig3, LV[[2]], FALSE, FALSE) /
         mk_row(fig3, LV[[3]], TRUE,  TRUE)
fig3p <- fig3p + patchwork::plot_annotation(
  title = "Power under the alternative hypothesis",
  theme = ggplot2::theme(plot.title = ggplot2::element_text(
    face = "bold", size = 30, hjust = 0.5)))
ggsave(file.path(OUT, "Fig3_Power.png"), fig3p, width = 24, height = 21, dpi = 300)
# ============ Supplementary Figure 1: DE count and DE ratio, one figure ============
# Both vary a property of the DE set at a fixed design (50 per group, 200 samples), so
# they belong together. Row 1 varies HOW MANY genes are DE at a fixed 1:1 split; row 2
# varies the up:down SPLIT at a fixed total of 2,000. Facet labels carry which is which,
# so the two rows cannot be confused.
sup <- bind_rows(
  p4 %>% mutate(level = factor(paste0("Total DE = ", format(as.numeric(as.character(level)),
                                                            big.mark = ",", trim = TRUE)),
                               levels = paste0("Total DE = ", c("500","1,000","2,000","4,000")))),
  p5 %>% mutate(level = factor(paste0("Up:down = ", as.character(level)),
                               levels = paste0("Up:down = ", c("1:1","4:1","9:1","1:0")))))
sup$level <- factor(as.character(sup$level),
                    levels = c(paste0("Total DE = ", c("500","1,000","2,000","4,000")),
                               paste0("Up:down = ",  c("1:1","4:1","9:1","1:0"))))
sup$cond <- factor(sup$cond, levels = COND)
ggsave(file.path(OUT, "Supp1_Power_DE.png"),
  grouped_panel(sup, "level", nrow = 2,
    "Power by the number and directional balance of DE genes",
    NULL, "Power", base = 22, ytop = 1.06),
  width = 24, height = 15, dpi = 300)

cat("written:\n"); print(list.files(OUT))
