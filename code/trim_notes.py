#!/usr/bin/env python3
"""Replace each draft note with a trimmed version; drop notes on minor comments."""
import re, sys

P = "revision_repo/response/response_to_reviewers.tex"
s = open(P, encoding="utf8").read()

TRIM = {
"R1.1": r"""\textbf{XS:} easy -- agreed.

Two citations in the Abstract, both line~172 of \texttt{manuscript.tex}:
\verb|\citep{ritchie_limma_2015}| after \texttt{limma} and
\verb|\citep{bolstad_comparison_2003}| after QN. No others.""",

"R1.2": r"""\textbf{XS:} we chose cross-platform normalization methods first, then added other
batch-effect correction methods as comparison. \emph{To discuss:} should we add single-cell
methods such as MNN?

Same request as R2.1, R2.m1 and R3.1 -- write one justification for all four.

Proposed criterion: two tiers (five purpose-built cross-platform methods -- TDM, MMR, RNABC,
Rank-In, Shambhala2; five general batch-correction methods as comparators -- QN, Angel,
ComBat, limma, COCONUT), plus the operational requirement that \textbf{a method must return
a corrected gene $\times$ sample matrix}. Under that rule \texttt{mnnCorrect} qualifies;
\texttt{fastMNN}, Harmony and scVI return embeddings only and cannot enter the DE comparison.""",

"R1.3": r"""\textbf{XS:} should we add number of DE genes and the up:down ratio? \emph{To discuss:}
how to simulate the ratio.

The generator behind the current figures (\texttt{sim\_array\_seq\_balanced}) draws the
down-regulated effect \emph{independently} on rows $(n_{DE}{+}1):(2n_{DE})$, so an unequal
ratio only needs separate \texttt{n\_up} / \texttt{n\_down} arguments. (The older
\texttt{data\_simulation.R} instead \emph{copies} the up block with case/control swapped,
which forces 1:1 -- worth knowing if that arm is revisited.)

The catch: roughly twenty scripts hardcode the truth as ``the first $2\cdot n_{DE}$ genes'',
with \texttt{real\_up} $=$ \texttt{1:n\_DE}. All must be updated consistently.

Expect QN, Angel and Rank-In to degrade as imbalance grows -- they assume the overall
distribution is unchanged -- while gene-wise methods should not. That makes this a
substantive result, not just a robustness check.""",

"R1.4": r"""\textbf{XS:} do we need to add other methods? (easy)

\texttt{prediction\_simulation.R} uses Lasso via \texttt{glmnet}. Ridge
(\texttt{alpha=0}) and elastic net (\texttt{alpha=0.5}) are the same call with one argument
changed; random forest or SVM would test whether the conclusion depends on linearity.
\texttt{caret} and \texttt{pROC} are installed.

Keep separate from R2.6, which asks about DESeq2/edgeR for \emph{differential expression},
not prediction.""",

"R2.1": r"""\textbf{XS:} discuss with professor -- should we add single-cell methods?

All five cited PMIDs are single-cell; none is a bulk cross-platform method:
\textbf{29608177} = MNN; \textbf{33758076} = CellMixS, a batch-effect \emph{metric} not a
correction method (belongs to R2.8); \textbf{37991248} and \textbf{38778573} are scRNA-seq
clustering/deep-learning methods; \textbf{40623818} (\emph{Genome Research} 2025) reports
that scRNA-seq correction methods are \emph{poorly calibrated}, which corroborates our own
Type I error finding and is worth citing.

If we add one, \texttt{batchelor::mnnCorrect} is the candidate -- it returns a gene-level
matrix and meets the R1.2 criterion. Caveat to decide first: MNN assumes many observations
sharing latent structure, so on 20--200 bulk samples it may perform poorly; are we prepared
to report that?""",

"R2.2": r"""\textbf{XS:} discuss with professor. Emphasise that we focus on microarray and RNA-seq,
and that microarray data is still useful in comparative studies.

GEO (Sept 2026): \textbf{69{,}879} array series vs \textbf{144{,}323} sequencing series. So
array is genuinely legacy for \emph{new} data -- do not contest that. The defensible claim is
that $\sim$70{,}000 array series already exist, cannot be regenerated (consumed specimens,
10--20 year follow-up cohorts), and are usable only through cross-platform integration.

R2.2 asks for \emph{framing}, not new methods -- so pair this with the scope declaration
(R3.2a) rather than treating them as alternatives, and put the archive argument in the
\textbf{Introduction}, since the complaint is that the boundary appears only in the Discussion.""",

"R2.4": r"""\textbf{XS:} to discuss. Measured behaviour of the injection
(\texttt{sim\_array\_seq\_balanced}, $m=50$, $n_{DE}=1000$):

\begin{itemize}\itemsep1pt
\item \textbf{Array:} additive shift $+0.181$ / $-0.184$ on a 2--14 log scale, independent
      of baseline.
\item \textbf{RNA-seq: an additive \emph{count} shift, not a fold change.} Induced
      $\log_2$FC is anti-correlated with baseline ($\rho = -0.965$): $+1.18$ for the
      lowest-expressed genes vs $+0.02$ for the highest, so DE signal sits almost entirely
      in low-expression genes. This is plausibly the bias the reviewer suspects.
\item \textbf{No noise is added} -- columns are resampled with replacement and used as-is,
      giving 15/100 duplicate seq columns. The older generator did add $N(0,0.3)$ noise.
\item \texttt{sim\_seq[sim\_seq < 0] <- 0} floors negatives, pushing 145 $\to$ 190
      down-genes to all-zero across every case sample -- trivially detectable.
\end{itemize}

Also: \texttt{sim\_balanced\_type1.rda}, a Figure~1 input, has \textbf{no generating script}
in the package (relevant to R2.9).""",

"R2.5": r"""\textbf{XS:} to discuss. The methods do receive different inputs, so the benchmark
does compare pipelines. Two findings we should raise ourselves:

\begin{enumerate}\itemsep2pt
\item \textbf{The ``ComBat'' arm is not plain ComBat} -- it is within-platform quantile
      normalisation followed by ComBat. Its only difference from RNABC is whether QN is
      applied within platform or jointly, so that contrast isolates \emph{the QN step},
      not two correction algorithms.
\item \textbf{The Methods text does not match the code for limma.} The manuscript says
      RNA-seq for limma is $\log_2$-transformed and quantile-normalised within platform;
      the code applies a natural log with \textbf{no QN}. The QN step described there is
      what the ComBat block does.
\end{enumerate}

Fix: add an explicit preprocessing column to the methods table, and correct the limma
description.""",

"R2.6": r"""\textbf{XS:} asking for limma details, and asking to add DESeq2 and edgeR.

\textbf{limma:} only \texttt{removeBatchEffect} is used, on
$\mathrm{cbind}(\mathrm{array},\log(\mathrm{seq}+1))$. The DE call is a Wilcoxon test on the
corrected matrix, identically for every method -- so \texttt{lmFit}/\texttt{eBayes} play no
part in the reported DE performance. One sentence answers this.

\textbf{DESeq2/edgeR:} they are DE testing frameworks, not harmonisers, so they cannot
produce a corrected array+seq matrix and are excluded by the R1.2 criterion. But note
\textbf{DESeq2 is already in the paper}, inside the Meta arm -- the reviewer appears not to
have noticed. The genuine answer is R2.7's joint-regression baseline, which uses exactly
these frameworks with platform as a covariate.""",

"R2.7": r"""\textbf{XS:} adding joint modelling baseline.

Cheapest high-value addition in the decision, and it answers R2.6 at the same time: fit
$y_g \sim \mathrm{condition} + \mathrm{platform}$ on the \emph{uncorrected} pooled data --
no corrected matrix -- via limma, or DESeq2/edgeR on counts.

One point to settle: every current arm is scored by a common Wilcoxon test on the corrected
matrix, but a joint model supplies its own test (as the Meta arm already does). State that
these two arms are evaluated on their own inferential output; this also connects to R2.m9.

Cost is low -- no new packages, no correction step, reuses the existing simulated data.""",

"R2.8": r"""\textbf{XS:} adding quantitative evaluation.

The key word is \emph{separately}: batch removal and signal preservation must be scored on
\textbf{two axes}, so that a method which removes the platform effect by destroying all
structure scores well on one and badly on the other.

Batch removal: kBET, iLISI, silhouette on the platform label. Signal preservation: ARI/NMI
against known labels, cLISI, silhouette on the biological label. Distribution alignment
(currently histograms): per-gene KS or Wasserstein distance.

\textbf{CellMixS} (PMID 33758076) is one of R2.1's own citations and is a metric package,
not a correction method -- adopting it answers R2.8 and engages that citation. Adding
ARI/silhouette also fixes R3.m3.

\emph{Blocked} by the missing real-data code (same blocker as R2.3).""",

"R2.9": r"""\textbf{XS:} provide code. Repository:
\texttt{github.com/XuejunSun98/\allowbreak Evaluating-Cross-Platform-Batch-Correction-Methods-Revision}

\textbf{In place:} manuscript and response sources, environment setup with exact versions
(also answers R2.m5), SLURM template, working Linux runners for Rank-In and Shambhala2, MMR
verification scripts.

\textbf{Still missing:} the analysis scripts themselves (roughly 150 hardcoded absolute
paths to strip first); the \textbf{real-data pipeline} -- four of seven figures have no
generating script and \texttt{real\_all2.rda}, \texttt{GSE70353.rda},
\texttt{GSE135134.rda} are absent; and \texttt{sim\_balanced\_type1.rda}, also with no
script. This is the largest gap and it also blocks R2.3 and R2.8.

Both external tools now run on Linux and were validated against reference output (Rank-In to
$6.6\times10^{-12}$ with an identical DEG set; Shambhala2 to $9.8\times10^{-15}$).

\emph{The repository must be made public, with an archived DOI given the ``Problem solving
protocol'' designation.}""",

"R2.10": r"""\textbf{XS:} to discuss -- largely settled by accepting R3.2(a).

What remains is wording: replace ``best general-purpose method'' and ``strongest overall''
with a claim bounded by what was tested. The Abstract's final sentence is the specific
target, plus the title, Key Points, Discussion and Conclusions.

Stating the recommendation conditionally -- limma for DE and clustering, QN for prediction,
\emph{in bulk microarray/RNA-seq integration under the designs tested} -- is defensible and
still useful. Revisit once the new results (TCGA, joint baseline, MNN) are in, since they
may support stronger wording again.""",

"R3.1": r"""\textbf{XS:} refer to R1.2.

Two specifics R3 raises that R1 does not: \textbf{ComBat-seq} is already implemented in our
scripts but is not one of the ten -- cheap to promote, and it answers their example
directly. \textbf{XPN} is not cheap: \texttt{xpn()} ships in the MatchMixeR sources we
already source, but that package's \texttt{src/} was never included, so it fails at its
\texttt{.C("XPN\_MLE\_C\_C")} call.

The classic/recent mix they find confusing is explained by the two-tier structure in R1.2.""",

"R3.2": r"""\textbf{XS:} we should accept (a) -- declare the microarray--bulk RNA-seq scope in the
title, abstract and Methods.

The title is already close: adding one word -- ``\dots Microarray and \emph{Bulk} RNA-seq
Data Analysis'' -- makes the scope explicit. This also settles R2.10 and blunts the
single-cell pressure in R2.2 and R3.2, giving one coherent position: narrow the claim rather
than expand the benchmark.

Caveat: narrowing does not excuse omitting methods that are \emph{in} scope. ComBat-seq and
XPN are bulk cross-platform methods (R3.1), so we should still add what we reasonably can,
or the narrowing reads as avoidance.""",

"R3.m1": r"""\textbf{XS:} need discuss.

Verified by fingerprinting against the shipped results:

\begin{itemize}\itemsep1pt
\item \textbf{Figure~1 (Type I) and Figure~3 (power) use \emph{different} MMR
      implementations} -- coefficients applied to raw \texttt{d\_seq} in one and to
      $\log(\mathrm{seq}+1)$ in the other. Each reproduces its own published values exactly,
      so this is not ambiguity. It must be unified and disclosed.
\item \textbf{MMR ranks last or next-to-last under every variant tested}, at every sample
      size. The implementation choice moves it one rank at most. So the limma/MMR reversal
      the reviewer finds counterintuitive is \emph{robust}, which is the answer to give.
\item The old \texttt{integration\_comparison*} no-op bug (only the first $4m$ genes
      transformed) never reached the manuscript.
\end{itemize}

Neither variant is dimensionally coherent -- \texttt{MM(array, log(seq+1))} fits
array $\to$ log-seq, so mapping seq onto the array scale needs the inverse. A
direction-corrected fit collapses under MM's shrinkage (seq SD falls to 6\% of array), so no
well-behaved alternative was found.""",

"R3.m2": r"""\textbf{XS:} modify plot.

The red dashed line is $\alpha = 0.05$ -- confirmed in \texttt{Unbalanced\_Type\_I.R}
(\texttt{geom\_hline(yintercept = 0.05)}). It just needs declaring in the caption.

For the layout, facet by design only with method on the $x$ axis. \textbf{No simulation
rerun is needed} -- Figure~1 is assembled entirely from the small shipped
\texttt{d\_plot\_*.rda} files, so this is a pure plotting change.""",

"R3.m3": r"""\textbf{XS:} edit -- easy.

The caption fix is easy: Figure~2 says ``Sample clustering result'' but shows PCA scatter
plots with no clustering.

Adding COCONUT is tractable -- \texttt{simulation\_clustering\_plot.R} already builds panels
for all eleven methods including \texttt{p\_COCONUT}. But two things must be fixed first: it
loads \texttt{sim\_FDR3.rda}, a stage-1 output that was never shipped and must be
regenerated; and its ComBat block slices the \emph{uncorrected} matrix, so a regenerated
ComBat panel would show uncorrected data labelled as ComBat.""",
}

n_trim = 0
for tag, body in TRIM.items():
    pat = (r"(\\noindent\\textbf\{Comment " + re.escape(tag) +
           r"\}\n\\begin\{comment\}\n.*?\n\\end\{comment\}\n\\begin\{response\}\n)"
           r"\\begin\{note\}.*?\\end\{note\}(\n\\end\{response\})")
    s, k = re.subn(pat, lambda m: m.group(1) + "\\\\begin{note}\n" + body + "\n\\\\end{note}" + m.group(2),
                   s, count=1, flags=re.S)
    if k != 1:
        print("  WARN: could not trim", tag); continue
    n_trim += 1

# minor comments: drop notes entirely, restore the empty TODO block
n_drop = 0
for tag in ["R2.m1"]:
    pat = (r"(\\noindent\\textbf\{Comment " + re.escape(tag) +
           r"\}\n\\begin\{comment\}\n.*?\n\\end\{comment\}\n\\begin\{response\}\n)"
           r"\\begin\{note\}.*?\\end\{note\}(\n\\end\{response\})")
    s, k = re.subn(pat, lambda m: m.group(1) + "%% TODO" + m.group(2), s, count=1, flags=re.S)
    n_drop += k

open(P, "w", encoding="utf8").write(s)
print(f"trimmed {n_trim} notes, dropped {n_drop} minor note(s)")
