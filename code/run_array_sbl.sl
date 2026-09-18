#!/bin/bash
#SBATCH --job-name=sim_sbl
#SBATCH --array=1-800%40
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g
#SBATCH --time=05:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
# One (grid, level, seed) per task. Throttled to 40 concurrent: each task starts
# a MATLAB process and the site licence pool is shared.
set -euo pipefail
source /usr/share/lmod/lmod/init/bash
module load r/4.5.0 matlab/2024b
export R_LIBS_USER=~/R/x86_64-pc-linux-gnu-library/4.5
export OMP_NUM_THREADS=1
export TMPDIR=/work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/scratch_sbl
mkdir -p "$TMPDIR"
cd /work/users/x/u/xuejun1/Integration_paper_Sim
Rscript revision_repo/code/sbl_rerun_task.R
