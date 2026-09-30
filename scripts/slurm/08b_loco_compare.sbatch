#!/bin/bash
#SBATCH --job-name=hcc_08b_cmp
#SBATCH --partition=active
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x_%j.out
# LOCO comparison (miner3 env: miner + lifelines): sbatch scripts/slurm/08b_loco_compare.sbatch
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source "${CONDA_BASE:-/users/sturkars/mambaforge}/etc/profile.d/conda.sh"
export PYTHONUTF8=1 MPLCONFIGDIR="$PWD/envs/cache/matplotlib" XDG_CACHE_HOME="$PWD/envs/cache"
conda activate miner3
python scripts/08b_loco.py compare "$@"
