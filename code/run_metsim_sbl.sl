#!/bin/bash
#SBATCH --job-name=metsim_sbl
#SBATCH --time=12:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --output=/work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/logs/metsim_sbl_%j.log
source /usr/share/lmod/lmod/init/bash
module load r/4.5.0 matlab/2024b
export OMP_NUM_THREADS=1
Rscript /work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/code/metsim_sbl_prediction.R
