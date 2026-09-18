#!/bin/bash
#SBATCH --job-name=ratio_cheap
#SBATCH --array=101-200
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
# Up:down ratio sweep, total DE fixed at 2000, one seed per task.
# Evaluated as in Soneson & Delorenzi 2013: the SAME BH-adjusted 0.05 threshold
# for every method, reporting AUC / TPR / true FDR, plus this manuscript's own
# raw-p Type I and power. Includes the all-upregulated extreme (1:0).
set -euo pipefail
source /usr/share/lmod/lmod/init/bash
module load r/4.5.0
export R_LIBS_USER=~/R/x86_64-pc-linux-gnu-library/4.5
export OMP_NUM_THREADS=1
export METHODS="QN,Angel,TDM,MMR,ComBat,ComBat_seq,RNABC,limma,MNN"
cd /work/users/x/u/xuejun1/Integration_paper_Sim
Rscript revision_repo/code/de_ratio_task.R
