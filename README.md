# Cross-Platform Integration of Microarray and Bulk RNA-seq Data

Code and aggregated results for a benchmark of thirteen batch-effect correction methods
for integrating bulk microarray and RNA-seq data, evaluated on Type I error, statistical
power, distributional alignment, sample clustering and cross-platform outcome prediction,
across four paired real datasets and a plasmode simulation.

> Manuscript under review. Citation details will be added on acceptance.

```
code/         preprocessing, simulation, per-method application, figures, SLURM scripts
results*/     aggregated per-replicate results — every figure rebuilds from these alone
figures/      the figures as published
retarget.sh   rewrite the absolute paths in code/ to your own location
```

The manuscript source and the response to reviewers are not included here.

## Methods compared

Thirteen correction methods plus a meta-analysis arm. Details of each algorithm are in
Table 1 of the manuscript; this is the roster and where each one is implemented in `code/`.

| Class | Method | Implementation | Original application domain |
|---|---|---|---|
| Unsupervised subject-wise | Quantile normalization (QN) | `preprocessCore` 1.70.0 | within-platform microarray normalization |
| | Angel's method | base R 4.5.0 | cross-platform transcriptome atlas |
| | TDM | `TDM` 0.3 | cross-platform microarray–RNA-seq for ML |
| Unsupervised gene-wise | MatchMixeR (MMR) | `MatchMixeR` 0.1.1 | cross-platform, matched samples |
| | ComBat | `sva` 3.56.0 | within-platform microarray batch correction |
| | ComBat-seq | `sva` 3.56.0 | RNA-seq count batch correction |
| | RNABC | `preprocessCore` + `sva` | cross-platform subtype transfer |
| | Shambhala2 | MATLAB R2024b | cross-platform harmonization to a universal format |
| | limma | `limma` 3.64.1 | differential expression framework |
| | XPN | `MatchMixeR` 0.1.1 | cross-platform merging of two studies |
| | MNN | `batchelor` 1.24.0 | cross-batch integration of single-cell RNA-seq |
| Supervised | COCONUT | `COCONUT` 1.0.2 | multi-study co-normalization using controls |
| | Rank-In | Python 3.9 | cross-platform microarray–RNA-seq for cancer |
| Meta-analysis | Cauchy combination (ACAT) | `ACAT` | combines per-study p-values |

The gene-wise methods are wrapped in [`code/sim_rerun_common.R`](code/sim_rerun_common.R)
(`METHODS`, from line 205), which every analysis script calls, so each is parameterised in
exactly one place. [`code/real_silhouette_all.R`](code/real_silhouette_all.R) is the one
script that applies all thirteen in a single pass, including the two that need external
tools, and is the easiest place to read the whole set side by side.

Two need external tools, and `code/longleaf_external_tools.R` provides drop-in replacements
that run both on Linux, each validated against the original implementation's own reference
output. Rank-In is distributed only as a Windows executable; Shambhala2 requires MATLAB.
Neither tool is redistributed here — see that file's header for how to obtain them.

## Parameters used for each method

Every argument that departs from the package default, taken from the code linked in the
table above. Where a method needs the two platforms on a particular scale, that is given
too, because several of these choices change the result materially.

| Method | Call | Non-default arguments | Input scale |
|---|---|---|---|
| QN | `normalize.quantiles.use.target()` ([code](code/real_silhouette_all.R#L114)) | `target` = first microarray column | array RMA log2; seq log2(TPM+1) |
| Angel | base R ([code](code/real_silhouette_all.R#L117)) | `rank(x)/length(x)` within each sample | array RMA log2; seq log2(TPM+1) |
| TDM | `tdm_transform()` ([code](code/real_silhouette_all.R#L119)) | `ref_data` = microarray, `target_data` = RNA-seq | array RMA log2; seq rounded to integer |
| MMR | `MatchMixeR::MM(a, lg)` ([code](code/sim_rerun_common.R#L93)) | defaults | array RMA log2; seq log(TPM+1) |
| ComBat | `sva::ComBat()` ([code](code/sim_rerun_common.R#L162)) | `mod = NULL` | both quantile-normalized **within** platform first |
| ComBat-seq | `sva::ComBat_seq()` ([code](code/sim_rerun_common.R#L110)) | `group = NULL` | array `2^x` then depth-matched to the seq median library size and rounded; seq rounded |
| RNABC | `normalize.quantiles.use.target()` then `sva::ComBat()` ([code](code/sim_rerun_common.R#L189)) | `target` = `rowMeans` of the microarray; `mod = NULL` | array RMA log2; seq log(TPM+1) |
| Shambhala2 | `Shambhala2(Input, P0, Q0)` ([code](code/real_silhouette_all.R#L130)) | `k = 5`, `delete_buffer_files = TRUE` | `Q0` = microarray, `P0` = `Input` = **raw** seq; output already on the array scale |
| limma | `limma::removeBatchEffect()` ([code](code/sim_rerun_common.R#L172)) | `design = NULL` | both quantile-normalized **within** platform first |
| XPN | `MatchMixeR::xpn()` ([code](code/sim_rerun_common.R#L199)) | defaults; corrects **both** platforms | array RMA log2; seq log(TPM+1) |
| MNN | `batchelor::mnnCorrect()` ([code](code/sim_rerun_common.R#L102)) | `k = 20`, `cos.norm.in = FALSE`, `cos.norm.out = FALSE` | array RMA log2; seq log(TPM+1) |
| COCONUT | `COCONUT()` ([code](code/real_silhouette_all.R#L154)) | `control.0.col = "group"`, `byPlatform = FALSE` | array RMA log2; seq log(TPM+1) |
| Rank-In | vendor program via `RankIn()` ([code](code/real_silhouette_all.R#L138)) | vendor defaults; leading `gene` column, tab-separated, unquoted | array RMA log2; seq **raw** counts |
| Meta (ACAT) | `ACAT()` ([code](code/method_smoke_test.R#L228)) | equal weights, `w_k = 1/K`; NA p-values dropped per gene | per-study p-values |

Shared settings: differential expression is tested with `wilcox.test()` on the pooled
corrected matrix; prediction uses `cv.glmnet(family = "binomial", nfolds = 5)` at
`lambda.min` with `alpha` 1, 0.5 and 0 for lasso, elastic net and ridge; 100 replicates per
simulation cell, seeds 101–200.

Four of these are worth singling out, because the obvious alternative gives a different
answer:

- **ComBat and ComBat-seq take `mod`/`group = NULL`.** Passing the outcome leaks it into the
  corrected matrix. On TCGA-LUSC, passing it produced a perfect in-sample fit and an
  *inverted* test AUC (0.046); with `group = NULL` the same model gives 0.538.
- **ComBat-seq's microarray input must be depth-matched.** `2^x` alone leaves the microarray
  with a median library size 3.9x the RNA-seq, which ComBat-seq reads as real sequencing
  depth and barely corrects.
- **Shambhala2 takes raw counts and its output is used as is.** It maps onto `Q0`'s
  distribution, and `Q0` is the log2 microarray, so the harmonized matrix is already on the
  array scale; logging it again compresses the result.
- **RNABC follows the published recipe**: quantile-normalize the RNA-seq onto the microarray
  row means, then ComBat with `mod = NULL`.

## Datasets

Four paired microarray / RNA-seq datasets. None is redistributed; all are public.

| Dataset | Array | Seq | Genes | Outcome | Pairing | Source |
|---|---|---|---|---|---|---|
| SEQC/MAQC | 8 | 16 | 14,732 | pool A vs B | replicates of the same two reference RNAs | GSE56457, GSE47774 |
| CCLE | 107 | 105 | 16,244 | breast vs large intestine | 101 shared cell lines | CCLE portal |
| METSIM | 331 | 331 | 16,160 | Matsuda index < 4 vs ≥ 4 | same subjects | GSE70353, GSE135134 |
| TCGA-LUSC | 132 | 497 | 11,894 | stage I vs II+ | 129 shared barcodes | GDC |

The simulations are plasmode: every replicate resamples the METSIM **control** samples
(`sim_BMI_control0.rda`) and injects fold changes, so the null retains real covariance
rather than assuming a parametric form.

## What builds each figure and table

Every figure rebuilds from the aggregated results in this repository; none of them needs
the primary data.

| Manuscript item | Built by | Reads |
|---|---|---|
| Figure 1 — Type I error | `code/sim_rerun_plots.R` | `results/` |
| Figure 2 — PCA, supervised methods | `code/fig2_pca.R` | simulates in-script |
| Figure 3 — Power | `code/sim_rerun_plots.R` | `results/` |
| Figure 4 — Prediction (simulation) | `code/fig4_prediction.R` | `results_pred_sweep/` |
| Figure 5 — Real-data metrics | `code/fig5_real_metrics.R` | `results_sil_all/` |
| Table 2 — METSIM prediction | `code/metsim_prediction.R`, `code/metsim_sbl_prediction.R` | → `results/metsim_prediction*.csv` |
| Table 3 — TCGA-LUSC prediction | `code/lusc_prediction.R` | → `results/lusc_prediction.csv` |
| Figure S1 — power vs DE composition | `code/sim_rerun_plots.R` | `results_de_sweep/`, `results_de_ratio/` |
| Figure S2 — quantile normalization on/off | `code/supp2_noqn.R` | `results_noqn_grids/` |
| Figure S3 — one-step vs two-step | `code/supp2_noqn.R` | `results_limma_1step/` |
| Figure S4 — distribution alignment | `code/fig_hist_plot.R` | `results/real_hist_density.rds`, `results_hist_sup/` |
| Figure S5 — clustering tiles | `code/supp5_cluster_tile.R` | `results_sil_all/` |

## Reproducing the analyses

Three stages. Only the third is cheap.

1. **Inputs** — `Rscript code/prepare_inputs.R --verify` (see *Obtaining the inputs*).
2. **Analyses** — SLURM array jobs, one `.sl` per grid in `code/`; each task writes one CSV
   so a partial run is resumable. The simulations are 100 replicates per cell; the
   Shambhala2 and Rank-In arms dominate the wall time because both cost per sample.
   Collate with `code/sim_rerun_collate.R`.
3. **Figures** — the plot scripts above, from the committed results.

Random seeds are fixed in every script (`set.seed`, and seed 101–200 per replicate), so a
rerun reproduces the published numbers rather than merely resembling them.

## Environment

R 4.5.0 with Bioconductor packages installed by `code/install_deps.R` and
`code/install_deps2.R`; MATLAB R2024b for Shambhala2; Python 3.9 for Rank-In. Per-method
versions are in the methods table above and in Table 1 of the manuscript.

Three constraints are load-bearing, and each cost a debugging session:

- **`preprocessCore` must be built without threading.** The stock build dies with
  `pthread_create() is 22` inside `normalize.quantiles.use.target()`, taking out QN and
  RNABC. `OMP_NUM_THREADS=1` does not help; reinstall with
  `configure.args = c(preprocessCore = "--disable-threading")`.
- **Rank-In needs numpy 1.21 / pandas 1.3 / scipy 1.7.** With numpy ≥ 2 the vendor's code
  dies where `np.mat` was removed.
- **ComBat-seq must not be given the outcome labels** in a cross-platform transfer. Doing
  so leaks the outcome into the corrected matrix and inverts the test AUC; `METHODS`
  passes `group = NULL`.

## Obtaining the inputs

No primary data is redistributed here. Everything analysed is public: GSE70353 and
GSE135134 (METSIM), GSE56457 and GSE47774 (SEQC/MAQC), CCLE from the Broad portal, and
TCGA-LUSC from the GDC.

`code/prepare_inputs.R` documents each intermediate file the analysis scripts load — what
it contains, which public sources it comes from, and how it is rebuilt — and verifies a
rebuilt file against the dimensions the manuscript was produced with:

```bash
Rscript code/prepare_inputs.R --verify
```

Two of the four have builders in this repository (`preprocess_metsim.R`,
`preprocess_lusc.R` + `build_lusc.R`); the other two are assembled from downloads as the
header of `prepare_inputs.R` describes. **The figures need none of them** — they read only
the aggregated results in `results*/`.

## Running the code on another machine

The scripts carry absolute paths from the cluster the analysis ran on. Retarget them in
one pass, then read `code/prepare_inputs.R` for the inputs:

```bash
git clone <this repo> /your/path/Integration_paper_Sim/revision_repo
cd /your/path/Integration_paper_Sim/revision_repo
./retarget.sh /your/path/Integration_paper_Sim
Rscript code/prepare_inputs.R            # what is present, what is missing
```

The scripts expect this repository at `<analysis_root>/revision_repo`, with the input data
files directly under `<analysis_root>`. Rebuilding the figures alone needs no inputs and no
retargeting beyond the clone path.

Analyses were run through SLURM; each `code/run_*.sl` submits one grid, and every task
writes its own CSV so a partial run is resumable.
