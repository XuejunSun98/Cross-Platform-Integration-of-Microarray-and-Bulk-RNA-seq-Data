#!/bin/bash
#SBATCH --job-name=histmmr
#SBATCH --time=02:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=2
#SBATCH --output=/work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/logs/hist_mmr.log
source /usr/share/lmod/lmod/init/bash
module load r/4.5.0
export OMP_NUM_THREADS=1
Rscript /work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/code/real_hist_mmr.R
