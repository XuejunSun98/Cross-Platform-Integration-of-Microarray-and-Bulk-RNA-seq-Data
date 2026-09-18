#!/bin/bash
#SBATCH --job-name=luscpred
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=/work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/logs/lusc_prediction.log
source /usr/share/lmod/lmod/init/bash
module load r/4.5.0
export OMP_NUM_THREADS=1
Rscript /work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/code/lusc_prediction.R
