# prepare_inputs.R -- how to obtain every input the analysis scripts load.
#
# No primary data is redistributed with this repository. Everything analysed here is
# public: the two METSIM series and the SEQC/MAQC series from GEO, CCLE from the Broad
# portal, and TCGA-LUSC from the GDC. This script records, for each intermediate file the
# analysis scripts `load()`, exactly what it contains and how it is rebuilt from those
# sources, and verifies a rebuilt file against the dimensions used in the manuscript.
#
#   Rscript code/prepare_inputs.R            # report what is present and what is missing
#   Rscript code/prepare_inputs.R --verify   # additionally check contents against the manifest
#
# Four files are needed. Two have builders in this repository; two are assembled from
# downloads and are described below in enough detail to reconstruct them.
#
# ---------------------------------------------------------------------------------------
# 1. real_data_METSIM/METSIM_analysis.rda        builder: code/preprocess_metsim.R
#    GSE70353 (microarray, Affymetrix HG-U219 / GPL13667, RMA log2 from the GEO series
#    matrix, probes mapped with hgu219.db) and GSE135134 (RNA-seq TPM). Probes are
#    collapsed to gene symbols by highest row median, the two platforms are intersected on
#    symbol, and the Matsuda Index is carried across by subject title. Download the two
#    series matrices into real_data_METSIM/raw/ first.
#
# 2. real_data_TCGA/LUSC_analysis.rda            builders: code/preprocess_lusc.R,
#                                                          code/build_lusc.R
#    TCGA-LUSC microarray and RNA-seq from the GDC, same collapse-by-median and
#    intersect-on-symbol treatment, with pathologic stage as the outcome. `paired` records
#    the 129 barcodes present on both platforms.
#
# 3. real_all2.rda                               assembled from three sources
#    The three-dataset object the real-data scripts read. Contents, with the dimensions
#    this manuscript was produced with:
#      d_affy_BMI  16160 x 331   METSIM microarray, log2 RMA        \  the paired subject
#      d_seq_BMI   16160 x 331   METSIM RNA-seq, TPM                /  subset of item 1
#      meta_BMI      331 x 80    GSE70353 sample characteristics (carries `matsuda:ch1`)
#      d_affy_CCLE 16244 x 107   CCLE microarray, log2               \ breast and large
#      d_seq_CCLE  16244 x 105   CCLE RNA-seq                        / intestine lines
#      d_affy_SEQC 14732 x   8   SEQC/MAQC microarray (GSE56457), A1-A4, B1-B4
#      d_seq_SEQC  14732 x  16   SEQC/MAQC RNA-seq (GSE47774), 8 lanes per pool
#    The METSIM block here is the 16,160-gene version behind the published results. It is
#    NOT the output of preprocess_metsim.R, which rebuilds METSIM from the current GEO
#    annotation and yields 18,566 genes on 330 subjects; the two differ in annotation
#    vintage and in the probe-collapsing rule, so the manuscript's METSIM numbers come from
#    this file and not from the rebuild. Use the rebuild to check that a conclusion survives
#    reannotation, not to reproduce a published number. For CCLE, take the microarray
#    and RNA-seq profiles for the breast and large-intestine lines from the CCLE portal,
#    collapse to gene symbols and intersect the two platforms; sample names encode
#    tissue and replicate (`breast_array_1`, `breast_seq_1`, ...), which is what `pkey()`
#    in real_silhouette_all.R matches on. For SEQC, take the two MAQC reference pools
#    (A = UHRR, B = HBRR) from GSE56457 and GSE47774 and intersect on symbol.
#
# 4. sim_BMI_control0.rda                        the plasmode seed for every simulation
#    Two control-only METSIM matrices, NOT a paired subset:
#      GSE70353_control  16160 x 101   log2 RMA microarray, range 2.02 .. 13.99
#      GSE135134_control 16160 x  65   RNA-seq TPM, range 0 .. 42458
#    Control means Matsuda Index >= 4. Columns are named `X<subject>_array` and
#    `X<subject>_seq`. Every simulation resamples these columns with replacement, adds
#    Gaussian noise and injects fold changes into the first 2*n_DE rows, so the simulation
#    cannot be rerun without this file; the aggregated per-replicate results in results/
#    are provided so the figures can be rebuilt without it.
#
# Figures 1 and 3 additionally carry two arms from the original analysis,
# simulation_2024/d_plot_Type_I_error_small.rda and d_plot_Power_Sample_Size.rda, which
# sim_rerun_plots.R reads for the methods that were not rerun.
# ---------------------------------------------------------------------------------------
ROOT <- "/work/users/x/u/xuejun1/Integration_paper_Sim"

MANIFEST <- list(
  list(path = "real_data_METSIM/METSIM_analysis.rda", builder = "code/preprocess_metsim.R",
       objects = c(d_array = "18566 x 330", d_seq = "18566 x 330")),
  list(path = "real_data_TCGA/LUSC_analysis.rda", builder = "code/preprocess_lusc.R + code/build_lusc.R",
       objects = c(d_array = "11894 x 132", d_seq = "11894 x 497", paired = "129")),
  list(path = "real_all2.rda", builder = "see item 3 above",
       objects = c(d_affy_BMI = "16160 x 331", d_seq_BMI = "16160 x 331", meta_BMI = "331 x 80",
                   d_affy_CCLE = "16244 x 107", d_seq_CCLE = "16244 x 105",
                   d_affy_SEQC = "14732 x 8", d_seq_SEQC = "14732 x 16")),
  list(path = "sim_BMI_control0.rda", builder = "see item 4 above",
       objects = c(GSE70353_control = "16160 x 101", GSE135134_control = "16160 x 65")))

verify <- "--verify" %in% commandArgs(TRUE)
cat(sprintf("%-42s %-9s %s\n", "input", "present", "rebuild with"))
cat(strrep("-", 96), "\n")
missing <- character()
for (m in MANIFEST) {
  ok <- file.exists(file.path(ROOT, m$path))
  if (!ok) missing <- c(missing, m$path)
  cat(sprintf("%-42s %-9s %s\n", m$path, if (ok) "yes" else "NO", m$builder))
  if (ok && verify) {
    e <- new.env(); load(file.path(ROOT, m$path), envir = e)
    for (nm in names(m$objects)) {
      if (!exists(nm, e)) { cat("    MISSING OBJECT:", nm, "\n"); next }
      x <- get(nm, e)
      got <- if (is.null(dim(x))) as.character(length(x)) else paste(dim(x), collapse = " x ")
      exp <- m$objects[[nm]]
      cat(sprintf("    %-20s %-14s %s\n", nm, got,
                  if (is.na(exp)) "" else if (identical(got, exp)) "ok" else paste("EXPECTED", exp)))
    }
  }
}
if (length(missing)) {
  cat("\n", length(missing), "input(s) missing. Rebuild them as described in the header,",
      "\nor rebuild only the figures, which read the aggregated results in results/ and",
      "\nneed none of these files:\n", sep = "")
  cat("  Rscript code/sim_rerun_plots.R    # Figures 1 and 3\n",
      "  Rscript code/fig5_real_metrics.R  # Figure 5\n",
      "  Rscript code/fig_hist_plot.R      # Supplementary Figure S4\n",
      "  Rscript code/supp5_cluster_tile.R # Supplementary Figure S5\n", sep = "")
} else cat("\nAll inputs present.\n")
