#!/bin/bash
#SBATCH --job-name=sim_cbseq
#SBATCH --array=101-200
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=01:30:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
# Re-run ComBat-seq only, with group = NULL. The 100-seed results in
# results/cheap_seed*.csv came from the group= version, which inverts transfer
# when the platforms carry different class proportions (see RUNLIST).
# Right-sized from the measured cost: ComBat-seq was 17-64 s per dataset.
set -euo pipefail
source /usr/share/lmod/lmod/init/bash
module load r/4.5.0
export R_LIBS_USER=~/R/x86_64-pc-linux-gnu-library/4.5
export OMP_NUM_THREADS=1
cd /work/users/x/u/xuejun1/Integration_paper_Sim
Rscript revision_repo/code/sim_rerun_task.R cbseq ComBat_seq
