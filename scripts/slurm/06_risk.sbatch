#!/bin/bash
#SBATCH --job-name=hcc_06
#SBATCH --partition=active
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out
# Submit from the project root: sbatch scripts/slurm/06_risk.sbatch [--skip prognostic]
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source "${CONDA_BASE:-/users/sturkars/mambaforge}/etc/profile.d/conda.sh"
export PYTHONUTF8=1
export MPLCONFIGDIR="$PWD/envs/cache/matplotlib"
export XDG_CACHE_HOME="$PWD/envs/cache"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
conda activate miner3
# xgboost 1.7.6 installed in the project (MINER's generatePredictor needs the pre-2.0 fit() API):
#   pip install --no-deps --target envs/pylib-miner xgboost==1.7.6
export PYTHONPATH="$PWD/envs/pylib-miner${PYTHONPATH:+:$PYTHONPATH}"
python scripts/06_risk.py "$@"
