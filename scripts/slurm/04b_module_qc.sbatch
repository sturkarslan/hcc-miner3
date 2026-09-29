#!/bin/bash
#SBATCH --job-name=hcc_04b
#SBATCH --partition=active
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
# Submit from the project root: sbatch scripts/slurm/04b_module_qc.sbatch
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source "${CONDA_BASE:-/users/sturkars/mambaforge}/etc/profile.d/conda.sh"
export PYTHONUTF8=1
conda activate "$PWD/envs/hcc-prep"
python scripts/04b_module_qc.py "$@"
