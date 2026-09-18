#!/bin/bash
#SBATCH --job-name=integ_sim
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32g
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Usage: sbatch run_r.sl <script.R> [args...]
#   e.g. sbatch run_r.sl simulation_2024/data_simulation.R
# Stage-1 data_simulation*.R scripts write 1-8 GB .rda files; raise --mem/--time for those.

set -euo pipefail
source /usr/share/lmod/lmod/init/bash
module load r/4.5.0

export R_LIBS_USER=~/R/x86_64-pc-linux-gnu-library/4.5
cd /work/users/x/u/xuejun1/Integration_paper_Sim

echo "host=$(hostname) start=$(date)"
Rscript "$@"
echo "end=$(date)"
