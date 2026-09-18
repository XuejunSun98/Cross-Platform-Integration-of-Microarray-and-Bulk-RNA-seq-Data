#!/bin/bash
#SBATCH --job-name=msplit
#SBATCH --array=1-11
#SBATCH --time=10:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --output=/work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/logs/msplit_%a.log
source /usr/share/lmod/lmod/init/bash
module load r/4.5.0
export OMP_NUM_THREADS=1
Rscript /work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/code/metsim_split_task.R
