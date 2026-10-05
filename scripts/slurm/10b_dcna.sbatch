#!/bin/bash
#SBATCH --job-name=hcc_10b
#SBATCH --partition=active
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
# Drug-constrained network activity (miner3 env): sbatch scripts/slurm/10b_dcna.sbatch
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source "${CONDA_BASE:-/users/sturkars/mambaforge}/etc/profile.d/conda.sh"
export PYTHONUTF8=1 MPLCONFIGDIR="$PWD/envs/cache/matplotlib" XDG_CACHE_HOME="$PWD/envs/cache"
conda activate miner3
python scripts/10b_dcna.py "$@"
