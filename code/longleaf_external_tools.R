# ============================================================
# Drop-in replacements for the two external tools, for Longleaf.
#
#   source("longleaf_external_tools.R")
#
# Both were validated against reference output on 2026-09-01 (see CLAUDE.md).
# Prerequisites, once per shell / at the top of an sbatch script:
#     source /usr/share/lmod/lmod/init/bash
#     module load r/4.5.0 matlab/2024b
# ============================================================

ROOT <- "/work/users/x/u/xuejun1/Integration_paper_Sim"

# ------------------------------------------------------------
# Rank-In
# ------------------------------------------------------------
# The vendor ships only a Windows .exe, but it is a PyInstaller bundle of a
# pure-Python 3.9 script. run_rankin.py executes the vendor's own code object
# on Linux -- no Wine, no reimplementation. Library versions MUST stay pinned
# to the 2021-era set in rankin_env39 (numpy 1.21 / pandas 1.3 / scipy 1.7):
# numpy >= 2 removes np.mat and the program dies at line 214.
#
# Replaces the processx/"Rank-In.exe" RankIn() helper in
# integration_comparison_Rank_in.R, integration_comparison_FDR_Rank_in.R,
# integration_comparison_effect_size_Rank_in.R, Unbalanced_Type_I_RankIn.R.
# Those scripts write d_all_Rankin.txt + sample.txt into the working directory
# first, then call RankIn(); keeping the same defaults means the call site does
# not have to change.
RankIn <- function(expr_file   = "d_all_Rankin.txt",
                   class_file  = "sample.txt",
                   result_file = "Rank_In_r_result.txt",
                   deg_file    = "deg.txt",
                   quiet       = TRUE) {
  py     <- file.path(ROOT, "rankin_env39/bin/python")
  runner <- file.path(ROOT, "rankin_dl/run_rankin.py")
  if (!file.exists(py))     stop("Rank-In python env missing: ", py)
  if (!file.exists(runner)) stop("Rank-In runner missing: ", runner)

  status <- system2(py, c(runner, expr_file, class_file, result_file, deg_file),
                    stdout = if (quiet) FALSE else "", stderr = if (quiet) FALSE else "")
  if (status != 0) stop("Rank-In failed (exit ", status, ")")
  if (!file.exists(result_file)) stop("Rank-In produced no result file")
  invisible(status)
}

# ------------------------------------------------------------
# Shambhala2
# ------------------------------------------------------------
# The .m files were missing from the package; they are the upstream ones from
# github.com/BorisovNM/Shambhala2, cloned to Shambhala2_upstream/. The only
# change needed versus the scripts' own copy is the system() line: upstream
# calls plain `matlab`, which works once the module is loaded.
#
# MATLAB must find Shambhala2.m / CuBlock.m / readExpressionData.m in the
# working directory, so this copies them in on first use.
shambhala2_setup <- function(dir = getwd()) {
  src <- file.path(ROOT, "Shambhala2_upstream")
  for (f in c("Shambhala2.m", "CuBlock.m", "readExpressionData.m")) {
    dst <- file.path(dir, f)
    if (!file.exists(dst)) file.copy(file.path(src, f), dst)
  }
  invisible(TRUE)
}

# Call MATLAB the Longleaf way. Use this in place of the hardcoded
#   system("/Applications/MATLAB_R2024a.app/bin/maci64/matlab ...")
#   system("\"C:\\Program Files\\MATLAB\\R2024a\\bin\\matlab.exe\" ...")
# lines inside the Shambhala2() helper in integration_comparison_shambhala*.R,
# Unbalanced_Type_I_shambhala.R, prediction_simulation_sbl.R,
# unbalanced_power_*_sbl.R and simulation_clustering_plot.R.
shambhala2_matlab <- function() {
  shambhala2_setup()
  if (Sys.which("matlab") == "") stop("matlab not on PATH -- run: module load matlab/2024b")
  status <- system("matlab -batch \"Shambhala2\"")
  if (status != 0) stop("MATLAB Shambhala2 failed (exit ", status, ")")
  # the .m writes Cu_bis.txt; the R wrapper polls for it
  invisible(status)
}
