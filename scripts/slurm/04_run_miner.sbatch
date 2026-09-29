#!/bin/bash
#SBATCH --job-name=hcc_04
#SBATCH --partition=active
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
# Submit from the project root: sbatch scripts/slurm/04_run_miner.sbatch [--matrix cohort_z] [--steps mechinf,subtypes]
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source "${CONDA_BASE:-/users/sturkars/mambaforge}/etc/profile.d/conda.sh"
export PYTHONUTF8=1
# MINER writes matplotlib figures; keep caches inside the project
export MPLCONFIGDIR="$PWD/envs/cache/matplotlib"
export XDG_CACHE_HOME="$PWD/envs/cache"
conda activate miner3
python scripts/04_run_miner.py "$@"
