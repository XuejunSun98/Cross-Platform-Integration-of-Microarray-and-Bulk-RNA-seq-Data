#!/bin/bash
#SBATCH --job-name=sim_cheap
#SBATCH --array=101-200
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=02:30:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
# One task per seed. 4 grids x 4 levels x 6 methods = 96 method-runs per task.
set -euo pipefail
source /usr/share/lmod/lmod/init/bash
module load r/4.5.0
export R_LIBS_USER=~/R/x86_64-pc-linux-gnu-library/4.5
export OMP_NUM_THREADS=1
cd /work/users/x/u/xuejun1/Integration_paper_Sim
Rscript revision_repo/code/sim_rerun_task.R cheap MMR,MNN,ComBat_seq,ComBat,limma,RNABC
