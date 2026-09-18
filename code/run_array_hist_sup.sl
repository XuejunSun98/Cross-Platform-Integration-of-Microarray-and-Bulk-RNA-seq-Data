#!/bin/bash
#SBATCH --job-name=histsup
#SBATCH --array=1-12%6
#SBATCH --time=11:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=2
#SBATCH --output=/work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/logs/histsup_%a.log
source /usr/share/lmod/lmod/init/bash
module load r/4.5.0 matlab/2024b
export OMP_NUM_THREADS=1
Rscript /work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/code/real_hist_supervised.R
