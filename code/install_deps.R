.libPaths("~/R/x86_64-pc-linux-gnu-library/4.5")
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 4)

try_install <- function(label, expr) {
  cat("\n========== ", label, " ==========\n")
  r <- try(eval(expr), silent = FALSE)
  cat("---> ", label, if (inherits(r, "try-error")) "FAILED" else "done", "\n")
}

try_install("CRAN: caret pROC", quote(install.packages(c("caret","pROC"))))
try_install("Bioc: sva",        quote(BiocManager::install("sva", ask = FALSE, update = FALSE)))
try_install("CRAN: COCONUT",    quote(install.packages("COCONUT")))
try_install("GH: ACAT",         quote(remotes::install_github("yaowuliu/ACAT", upgrade = "never")))
try_install("GH: TDM",          quote(remotes::install_github("greenelab/TDM", upgrade = "never")))
try_install("GH: MetaDE",       quote(remotes::install_github("metaOmics/MetaDE", upgrade = "never")))

cat("\n=== FINAL STATUS ===\n")
for (p in c("sva","COCONUT","TDM","MetaDE","ACAT","caret","pROC"))
  cat(sprintf("%-10s %s\n", p, if (requireNamespace(p, quietly=TRUE)) "OK" else "MISSING"))
