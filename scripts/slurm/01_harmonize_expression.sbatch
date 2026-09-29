#!/bin/bash
#SBATCH --job-name=hcc_01
#SBATCH --partition=active
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
# Submit from the project root: sbatch scripts/slurm/01_harmonize_expression.sbatch
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source "${CONDA_BASE:-/users/sturkars/mambaforge}/etc/profile.d/conda.sh"
export PYTHONUTF8=1
conda activate "$PWD/envs/hcc-prep"
python scripts/01_harmonize_expression.py "$@"
