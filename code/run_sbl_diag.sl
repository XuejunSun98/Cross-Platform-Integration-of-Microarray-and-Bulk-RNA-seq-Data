#!/bin/bash
#SBATCH --job-name=sbldiag
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=/work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/logs/sbl_diag.log
source /usr/share/lmod/lmod/init/bash
module load r/4.5.0 matlab/2024b
export OMP_NUM_THREADS=1
export TMPDIR=/work/users/x/u/xuejun1/Integration_paper_Sim/revision_repo/scratch_sbl
mkdir -p "$TMPDIR"
Rscript /tmp/claude-345312/-work-users-x-u-xuejun1-Integration-paper-Sim/aca0bda4-f319-426e-8386-01c564ef3aaf/scratchpad/sbl_diag.R
