#!/bin/bash
#SBATCH --job-name=de_xpn
#SBATCH --array=101-200
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
# R1.3 DE sweep: 12 scenarios (total DE x up:down ratio) at m=50, one seed per task.
set -euo pipefail
source /usr/share/lmod/lmod/init/bash
module load r/4.5.0
export R_LIBS_USER=~/R/x86_64-pc-linux-gnu-library/4.5
export OMP_NUM_THREADS=1
cd /work/users/x/u/xuejun1/Integration_paper_Sim
Rscript revision_repo/code/de_sweep_task.R xpn XPN
