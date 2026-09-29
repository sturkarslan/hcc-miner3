#!/bin/bash
#SBATCH --job-name=hcc_04c
#SBATCH --partition=active
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.out
# Submit from the project root after mechinf: sbatch scripts/slurm/04c_filter_regulons.sbatch
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source "${CONDA_BASE:-/users/sturkars/mambaforge}/etc/profile.d/conda.sh"
export PYTHONUTF8=1
conda activate "$PWD/envs/hcc-prep"
python scripts/04c_filter_regulons.py "$@"
