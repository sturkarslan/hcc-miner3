#!/bin/bash
#SBATCH --job-name=hcc_05
#SBATCH --partition=active
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
# Submit from the project root: sbatch scripts/slurm/05_causal_inference.sbatch [--features MUT_CTNNB1,MUT_TP53]
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source "${CONDA_BASE:-/users/sturkars/mambaforge}/etc/profile.d/conda.sh"
export PYTHONUTF8=1
export MPLCONFIGDIR="$PWD/envs/cache/matplotlib"
export XDG_CACHE_HOME="$PWD/envs/cache"
# one BLAS thread per worker process
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
conda activate miner3
python scripts/05_causal_inference.py "$@"
